#' Install ONNX Runtime native library
#'
#' Downloads and installs a prebuilt ONNX Runtime shared library for the
#' current platform. The library is used for ML model inference (CREPE,
#' Brouhaha, Swift-F0) without requiring Python.
#'
#' @param version ONNX Runtime version to install. Default: "1.24.3".
#' @param path Installation directory. Default: platform-appropriate user cache.
#' @param gpu Logical. If TRUE, install GPU-enabled variant (Linux/Windows only).
#'
#' @return Invisible path to installed library directory.
#'
#' @details
#' The ONNX Runtime shared library is downloaded from GitHub releases and
#' extracted to a user-level cache directory. This survives R package reinstalls.
#'
#' Supported platforms:
#' \itemize{
#'   \item macOS arm64 (Apple Silicon)
#'   \item Linux x86_64
#'   \item Linux aarch64
#'   \item Windows x64
#' }
#'
#' macOS x86_64 (Intel) is not supported by ONNX Runtime >= 1.24.1.
#'
#' @examples
#' \dontrun{
#' install_onnxruntime()
#' onnxruntime_available()
#' onnxruntime_info()
#' }
#'
#' @export
install_onnxruntime <- function(version = "1.24.3",
                                 path = NULL,
                                 gpu = FALSE) {

  # Determine platform
  platform <- .ort_detect_platform(gpu = gpu)

  # Determine install path
  if (is.null(path)) {
    path <- file.path(tools::R_user_dir("superassp", "cache"),
                      "onnxruntime")
  }
  lib_dir <- file.path(path, "lib")

  # Construct download URL
  ext <- if (platform$os == "win") "zip" else "tgz"
  archive_name <- sprintf("onnxruntime-%s-%s.%s", platform$slug, version, ext)
  url <- sprintf(
    "https://github.com/microsoft/onnxruntime/releases/download/v%s/%s",
    version, archive_name
  )

  cli::cli_h2("Installing ONNX Runtime v{version}")
  cli::cli_alert_info("Platform: {platform$slug}")
  cli::cli_alert_info("URL: {url}")


  # Create temp dir for download
  tmp_dir <- tempfile("ort_")
  dir.create(tmp_dir, recursive = TRUE)
  on.exit(unlink(tmp_dir, recursive = TRUE), add = TRUE)

  archive_path <- file.path(tmp_dir, archive_name)

  # Download
  cli::cli_alert("Downloading...")
  utils::download.file(url, archive_path, mode = "wb", quiet = TRUE)

  if (!file.exists(archive_path)) {
    cli::cli_abort("Download failed. Check your internet connection and the URL.")
  }

  # Extract
  cli::cli_alert("Extracting...")
  if (ext == "zip") {
    utils::unzip(archive_path, exdir = tmp_dir)
  } else {
    utils::untar(archive_path, exdir = tmp_dir)
  }

  # Find extracted directory (e.g., onnxruntime-osx-arm64-1.24.3/)
  extracted <- list.dirs(tmp_dir, recursive = FALSE, full.names = TRUE)
  # Filter to onnxruntime-* dirs
  ort_dirs <- extracted[grepl("^onnxruntime-", basename(extracted))]
  if (length(ort_dirs) == 0) {
    cli::cli_abort("Could not find extracted onnxruntime directory.")
  }
  ort_dir <- ort_dirs[1]

  # Create target directories
  dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)

  # Copy shared library
  lib_name <- switch(platform$os,
    "osx" = "libonnxruntime.dylib",
    "linux" = "libonnxruntime.so",
    "win" = "onnxruntime.dll"
  )

  src_lib <- file.path(ort_dir, "lib", lib_name)
  if (!file.exists(src_lib)) {
    # Some releases use flat layout
    src_lib <- file.path(ort_dir, lib_name)
  }
  if (!file.exists(src_lib)) {
    cli::cli_abort(
      "Could not find {lib_name} in extracted archive. Contents: {paste(list.files(ort_dir, recursive = TRUE), collapse = ', ')}"
    )
  }

  # Copy library and any versioned symlinks
  lib_files <- list.files(dirname(src_lib), pattern = "onnxruntime",
                          full.names = TRUE)
  file.copy(lib_files, lib_dir, overwrite = TRUE, copy.date = TRUE)

  # Verify
  dest_lib <- file.path(lib_dir, lib_name)
  if (!file.exists(dest_lib)) {
    cli::cli_abort("Installation failed: {dest_lib} not found after copy.")
  }

  # Store path in option for immediate use
  options(superassp.onnxruntime.path = lib_dir)

  # Write version file for later reference
  writeLines(version, file.path(path, "VERSION"))

  cli::cli_alert_success("ONNX Runtime v{version} installed to {lib_dir}")
  cli::cli_alert_info("Library: {dest_lib} ({round(file.info(dest_lib)$size / 1e6, 1)} MB)")

  invisible(lib_dir)
}


#' Check if ONNX Runtime is available
#'
#' @return Logical. TRUE if the ONNX Runtime shared library can be loaded.
#' @export
onnxruntime_available <- function() {
  # Ensure path is set from cached install
  .ort_set_cached_path()
  ort_available_cpp()
}


#' Get ONNX Runtime information
#'
#' @return A list with components:
#' \describe{
#'   \item{available}{Logical. Whether onnxruntime is loadable.}
#'   \item{version}{Character. Version string, or "" if not available.}
#'   \item{lib_path}{Character. Path to loaded library, or "" if not available.}
#'   \item{lib_dir}{Character. Configured library directory.}
#'   \item{platform}{List with os, arch, slug.}
#' }
#' @export
onnxruntime_info <- function() {
  .ort_set_cached_path()
  avail <- ort_available_cpp()
  list(
    available = avail,
    version = if (avail) ort_version_cpp() else "",
    lib_path = if (avail) ort_lib_path_cpp() else "",
    lib_dir = .ort_resolve_lib_dir(),
    platform = tryCatch(.ort_detect_platform(), error = function(e) list())
  )
}


#' Get ONNX Runtime library path
#'
#' @return Character. Path to the onnxruntime library directory.
#' @export
onnxruntime_path <- function() {
  .ort_resolve_lib_dir()
}


# ---------- Internal helpers ----------

#' Detect platform for ONNX Runtime download
#' @param gpu Logical. Request GPU variant.
#' @return List with os, arch, slug.
#' @keywords internal
.ort_detect_platform <- function(gpu = FALSE) {
  si <- Sys.info()
  os_name <- tolower(si[["sysname"]])
  arch <- si[["machine"]]

  os <- switch(os_name,
    "darwin" = "osx",
    "linux" = "linux",
    "windows" = "win",
    cli::cli_abort("Unsupported OS: {os_name}")
  )

  # Normalize architecture
  arch_norm <- switch(arch,
    "arm64" = "arm64",
    "aarch64" = "aarch64",
    "x86_64" = "x64",
    "AMD64" = "x64",
    cli::cli_abort("Unsupported architecture: {arch}")
  )

  # Check macOS x64 (dropped in ORT >= 1.24.1)
  if (os == "osx" && arch_norm == "x64") {
    cli::cli_abort(c(
      "macOS x86_64 (Intel) is not supported by ONNX Runtime >= 1.24.1.",
      "i" = "Use an Apple Silicon (arm64) Mac, or use Python onnxruntime instead."
    ))
  }

  # Build slug
  slug <- if (os == "osx") {
    paste0("osx-", arch_norm)
  } else if (os == "linux") {
    base <- paste0("linux-", arch_norm)
    if (gpu) paste0(base, "-gpu") else base
  } else {
    base <- paste0("win-", arch_norm)
    if (gpu) paste0(base, "-gpu") else base
  }

  list(os = os, arch = arch_norm, slug = slug)
}


#' Resolve ONNX Runtime library directory
#' @return Character path to library directory. May not exist.
#' @keywords internal
.ort_resolve_lib_dir <- function() {
  # 1. R option
  opt <- getOption("superassp.onnxruntime.path")
  if (!is.null(opt) && nzchar(opt)) return(opt)

  # 2. Environment variable
  env <- Sys.getenv("ONNXRUNTIME_LIB_PATH", "")
  if (nzchar(env)) return(env)

  # 3. Default cache location
  file.path(tools::R_user_dir("superassp", "cache"), "onnxruntime", "lib")
}


#' Set cached path if install exists but option is unset
#' @keywords internal
.ort_set_cached_path <- function() {
  # Only set if option not already configured
  if (!is.null(getOption("superassp.onnxruntime.path"))) return(invisible())

  lib_dir <- .ort_resolve_lib_dir()
  if (dir.exists(lib_dir)) {
    # Check that the library file actually exists
    si <- Sys.info()
    lib_name <- switch(tolower(si[["sysname"]]),
      "darwin" = "libonnxruntime.dylib",
      "linux" = "libonnxruntime.so",
      "windows" = "onnxruntime.dll",
      return(invisible())
    )
    if (file.exists(file.path(lib_dir, lib_name))) {
      options(superassp.onnxruntime.path = lib_dir)
      ort_set_lib_dir_cpp(lib_dir)
    }
  }
  invisible()
}

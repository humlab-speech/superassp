# onnxruntime.R — ONNX Runtime installation, availability, and session management
#
# All ONNX-using functions (trk_crepe, etc.) call ensure_onnx() at entry.
# The runtime is downloaded automatically on first use and cached in the
# R user directory (~/.cache/R/superassp/onnxruntime/).


# ── Public-internal entry point ───────────────────────────────────────────────

#' Ensure ONNX Runtime is available, installing automatically if needed
#'
#' Called internally by any function that requires ONNX Runtime inference
#' (e.g. \code{trk_crepe}). On first call the function checks for a cached
#' installation; if none is found it downloads and installs the runtime
#' automatically with an informative message.
#'
#' @param version ONNX Runtime version to install if not already present.
#'   Default: \code{"1.24.3"}.
#' @param gpu Logical. If TRUE, install GPU variant (Linux/Windows only).
#'
#' @return Invisible TRUE.
#' @keywords internal
ensure_onnx <- function(version = "1.24.3", gpu = FALSE) {
  # 1. Try loading from cached path
  .ort_set_cached_path()

  # 2. Already available — done
  if (ort_available_cpp()) return(invisible(TRUE))

  # 3. Not available — auto-install
  cli::cli_inform(c(
    "i" = "ONNX Runtime not found. Installing v{version} automatically...",
    "i" = "This is a one-time download (~30 MB) cached in your R user directory."
  ))

  tryCatch(
    .ort_install(version = version, gpu = gpu),
    error = function(e) {
      cli::cli_abort(c(
        "x" = "ONNX Runtime installation failed: {conditionMessage(e)}",
        "i" = "Check your internet connection and try again."
      ))
    }
  )

  # 4. Re-check after install
  .ort_set_cached_path()
  if (!ort_available_cpp()) {
    cli::cli_abort("ONNX Runtime was installed but could not be loaded.")
  }

  invisible(TRUE)
}


# ── Installation ──────────────────────────────────────────────────────────────

#' Download and install ONNX Runtime native library
#'
#' @param version ONNX Runtime version to install. Default: "1.24.3".
#' @param path Installation directory. Default: platform-appropriate user cache.
#' @param gpu Logical. If TRUE, install GPU-enabled variant (Linux/Windows only).
#'
#' @return Invisible path to installed library directory.
#' @keywords internal
.ort_install <- function(version = "1.24.3", path = NULL, gpu = FALSE) {

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


# ── Path helpers ──────────────────────────────────────────────────────────────

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
      "i" = "Use an Apple Silicon (arm64) Mac."
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


# ── Session helpers ───────────────────────────────────────────────────────────

#' Create an ONNX Runtime inference session
#'
#' @param model_path Path to an ONNX model file (.onnx).
#' @param num_threads Number of intra-op threads. 0 = auto (ORT default).
#'
#' @return An external pointer to the ORT session (class "ort_session").
#'
#' @details
#' The session holds the loaded model and is reusable across multiple
#' inference calls. It is automatically released when garbage collected.
#'
#' @examples
#' \dontrun{
#' sess <- ort_session("model.onnx")
#' result <- ort_run(sess, list(input = rnorm(1024)), list(c(1L, 1024L)))
#' }
#'
#' @keywords internal
ort_session <- function(model_path, num_threads = 0L) {
  ensure_onnx()

  if (!file.exists(model_path)) {
    cli::cli_abort("Model file not found: {model_path}")
  }

  ort_create_session_cpp(normalizePath(model_path, mustWork = TRUE),
                          as.integer(num_threads))
}


#' Run ONNX Runtime inference
#'
#' @param session An ORT session from \code{ort_session()}.
#' @param inputs Named list of input data. Each element is a numeric vector.
#' @param shapes List of integer vectors specifying the shape of each input.
#'   Must be in the same order as \code{inputs}.
#' @param output_names Character vector of output tensor names to fetch.
#'   NULL = fetch all outputs.
#'
#' @return Named list of output tensors (numeric vectors with "shape" attribute).
#'
#' @keywords internal
ort_run <- function(session, inputs, shapes, output_names = NULL) {
  if (!inherits(session, "ort_session")) {
    cli::cli_abort("session must be created by ort_session()")
  }

  input_names <- names(inputs)
  if (is.null(input_names) || any(input_names == "")) {
    cli::cli_abort("All inputs must be named.")
  }

  # Convert shapes to list of integer vectors
  shapes <- lapply(shapes, as.integer)

  ort_run_cpp(session,
              input_names,
              unname(inputs),
              shapes,
              output_names)
}


#' Get ORT session input/output metadata
#'
#' @param session An ORT session from \code{ort_session()}.
#' @return Named list of tensor info (name, shape, type).
#' @keywords internal
ort_input_info <- function(session) {
  if (!inherits(session, "ort_session")) {
    cli::cli_abort("session must be created by ort_session()")
  }
  ort_session_input_info_cpp(session)
}

#' @rdname ort_input_info
#' @keywords internal
ort_output_info <- function(session) {
  if (!inherits(session, "ort_session")) {
    cli::cli_abort("session must be created by ort_session()")
  }
  ort_session_output_info_cpp(session)
}

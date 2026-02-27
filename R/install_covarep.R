#' Install COVAREP Python Module
#'
#' Install the covarep Python module with optional optimizations.
#'
#' The COVAREP module provides speech analysis functions including:
#' \itemize{
#'   \item F0 estimation via SRH algorithm
#'   \item Glottal source analysis via IAIF
#'   \item Voice quality parameter extraction
#' }
#'
#' @param method Installation method:
#'   \describe{
#'     \item{\code{"auto"}}{Automatically detect and install optimizations (default)}
#'     \item{\code{"numba"}}{Install with Numba JIT optimization (5-10x speedup)}
#'     \item{\code{"pure"}}{Pure Python mode (no optimizations, maximum compatibility)}
#'   }
#' @param python Path to Python executable (default: NULL uses reticulate default)
#' @param verbose Logical; show installation progress (default: TRUE)
#' @param force Logical; force reinstallation (default: FALSE)
#'
#' @return Invisible TRUE on success
#'
#' @details
#' **Optimization Levels:**
#'
#' \describe{
#'   \item{\bold{Pure Python (1.0x)}}{
#'     - Requirements: Python 3.8+, NumPy, SciPy
#'     - Performance: Baseline
#'     - Use case: Maximum compatibility
#'   }
#'   \item{\bold{NumPy Vectorization (5-10x)}}{
#'     - Automatic with NumPy installation
#'     - No additional dependencies
#'     - Primary optimization layer
#'   }
#'   \item{\bold{Numba JIT (5-10x for LPC)}}{
#'     - Requirements: + numba >= 0.56.0
#'     - Performance: ~2-3x total speedup
#'     - Use case: Production, no compilation needed
#'     - Just-in-time compilation, automatic fallback
#'   }
#' }
#'
#' **Installation Methods:**
#'
#' \itemize{
#'   \item \code{method="auto"}: Tries Numba, falls back to pure if unavailable
#'   \item \code{method="numba"}: Forces Numba installation (recommended)
#'   \item \code{method="pure"}: Skips all optimizations
#' }
#'
#' **Platform Notes:**
#'
#' \itemize{
#'   \item macOS (Apple Silicon): Numba works with arm64 (set NUMBA_DISABLE_INTEL_SVML=1 if warnings)
#'   \item Linux: Full support, excellent performance
#'   \item Windows: Full support, may need Visual C++ for some dependencies
#' }
#'
#' @seealso
#' \code{\link{covarep_available}} to check installation status,
#' \code{\link{covarep_info}} for detailed optimization info,
#' \code{\link{trk_covarep_srh}} for F0 tracking,
#' \code{\link{trk_covarep_iaif}} for glottal analysis
#'
#' @examples
#' \dontrun{
#' # Automatic installation (recommended)
#' install_covarep()
#'
#' # Force Numba optimization
#' install_covarep(method = "numba")
#'
#' # Pure Python (no optimizations)
#' install_covarep(method = "pure")
#'
#' # Check installation
#' covarep_info()
#' }
#'
install_covarep <- function(method = c("auto", "numba", "pure"),
                            python = NULL,
                            verbose = TRUE,
                            force = FALSE) {

  method <- match.arg(method)

  # Check if already installed
  if (!force && covarep_available()) {
    info <- covarep_info()
    if (verbose) {
      cli::cli_alert_success("COVAREP already installed")
      cli::cli_alert_info("Optimization level: {info$optimization_level}")
      cli::cli_alert_info("Use force=TRUE to reinstall")
    }
    return(invisible(TRUE))
  }

  # Get Python executable
  if (is.null(python)) {
    py_config <- reticulate::py_config()
    python <- py_config$python
  }

  if (verbose) {
    cli::cli_h2("Installing COVAREP Python Module")
    cli::cli_alert_info("Python: {python}")
    cli::cli_alert_info("Method: {method}")
  }

  # Get package path
  pkg_path <- system.file("python", "covarep_python", package = "superassp")

  if (!dir.exists(pkg_path)) {
    stop("COVAREP module not found in package installation", call. = FALSE)
  }

  # Install core dependencies
  if (verbose) cli::cli_alert("Installing core dependencies...")

  core_deps <- c("numpy>=1.20.0", "scipy>=1.7.0", "soundfile>=0.10.0")

  tryCatch({
    reticulate::py_install(core_deps, pip = TRUE)
    if (verbose) cli::cli_alert_success("Core dependencies installed")
  }, error = function(e) {
    stop("Failed to install core dependencies: ", e$message, call. = FALSE)
  })

  # Install optional optimizations
  if (method %in% c("auto", "numba")) {
    if (verbose) cli::cli_alert("Installing Numba optimization...")

    tryCatch({
      reticulate::py_install("numba>=0.56.0", pip = TRUE)
      if (verbose) cli::cli_alert_success("Numba installed (5-10x speedup for LPC)")
    }, error = function(e) {
      if (method == "numba") {
        stop("Failed to install Numba: ", e$message, call. = FALSE)
      } else {
        if (verbose) {
          cli::cli_alert_warning("Numba installation failed, using pure Python")
          cli::cli_alert_info("Error: {e$message}")
        }
      }
    })
  }

  # Install COVAREP package
  if (verbose) cli::cli_alert("Installing COVAREP module...")

  tryCatch({
    # Use pip install in editable mode or direct install
    reticulate::py_install(pkg_path, pip = TRUE)
    if (verbose) cli::cli_alert_success("COVAREP module installed")
  }, error = function(e) {
    stop("Failed to install COVAREP: ", e$message, call. = FALSE)
  })

  # Reload module
  tryCatch({
    covarep_module <<- reticulate::import("covarep", delay_load = TRUE)
  }, error = function(e) {
    stop("COVAREP installed but failed to load: ", e$message, call. = FALSE)
  })

  # Report status
  if (verbose) {
    cli::cli_h2("Installation Complete")
    info <- covarep_info()

    cli::cli_alert_success("COVAREP is ready")
    cli::cli_ul(c(
      "Optimization level: {info$optimization_level}",
      "Numba JIT: {if(info$numba) 'Available' else 'Not available'}",
      "Performance: {if(info$numba) '5-10x speedup' else 'Baseline'}"
    ))

    cli::cli_alert_info("Use trk_covarep_srh() and trk_covarep_iaif() for analysis")
  }

  invisible(TRUE)
}


#' Check COVAREP Module Availability
#'
#' Check if the COVAREP Python module is installed and available.
#'
#' @return Logical; TRUE if covarep module is available, FALSE otherwise
#'
#' @examples
#' \dontrun{
#' if (covarep_available()) {
#'   f0 <- trk_covarep_srh("audio.wav", toFile = FALSE)
#' } else {
#'   message("Install COVAREP with: install_covarep()")
#' }
#' }
#'
#' @seealso \code{\link{install_covarep}}, \code{\link{covarep_info}}
covarep_available <- function() {
  !is.null(covarep_module) && reticulate::py_module_available("covarep")
}


#' Get COVAREP Module Information
#'
#' Get detailed information about COVAREP installation and optimization status.
#'
#' @return List with elements:
#'   \describe{
#'     \item{\code{available}}{Logical; module installed and available}
#'     \item{\code{numba}}{Logical; Numba JIT optimization available}
#'     \item{\code{optimization_level}}{Character; optimization level description}
#'     \item{\code{performance}}{Character; expected performance gain}
#'   }
#'
#' @examples
#' \dontrun{
#' # Get detailed information
#' info <- covarep_info()
#' print(info)
#'
#' # Check specific optimization
#' if (info$numba) {
#'   message("Numba optimization is active")
#' }
#' }
#'
#' @seealso \code{\link{install_covarep}}, \code{\link{covarep_available}}
covarep_info <- function() {
  if (!covarep_available()) {
    return(list(
      available = FALSE,
      numba = FALSE,
      optimization_level = "Not installed",
      performance = "N/A",
      message = "COVAREP module not installed. Use install_covarep()"
    ))
  }

  # Check Numba availability
  numba_available <- tryCatch({
    py_has_numba <- reticulate::py_eval("
try:
    import numba
    True
except ImportError:
    False
", convert = TRUE)
    py_has_numba
  }, error = function(e) FALSE)

  # Determine optimization level
  if (numba_available) {
    opt_level <- "NumPy vectorization + Numba JIT"
    performance <- "5-10x speedup (F0), 2-3x speedup (IAIF)"
  } else {
    opt_level <- "NumPy vectorization only"
    performance <- "5-10x speedup (F0 only)"
  }

  list(
    available = TRUE,
    numba = numba_available,
    optimization_level = opt_level,
    performance = performance
  )
}


# Module cache (set in .onLoad)
covarep_module <- NULL

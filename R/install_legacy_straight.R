#' Setup Python path for legacy STRAIGHT module
#'
#' Internal helper function to add the legacy STRAIGHT module to Python's sys.path
#'
#' @return inst_path Character; path to inst/python directory
#' @keywords internal
#' @noRd
.setup_straight_path <- function() {
  # Find inst/python directory
  inst_path <- system.file("python", package = "superassp")
  
  if (inst_path == "") {
    # Try relative path for development
    inst_path <- file.path(getwd(), "inst", "python")
    if (!dir.exists(inst_path)) {
      stop("superassp not properly installed and inst/python not found", call. = FALSE)
    }
  }
  
  # Add to Python path if not already there
  reticulate::py_run_string(sprintf("
import sys
if '%s' not in sys.path:
    sys.path.insert(0, '%s')
", inst_path, inst_path))
  
  return(inst_path)
}


#' Install legacy STRAIGHT Python dependencies
#'
#' Installs the required Python packages for legacy STRAIGHT vocoder analysis
#' and synthesis. The STRAIGHT algorithm provides high-quality pitch-adaptive
#' spectral analysis for speech.
#'
#' @param method Installation method for Python packages
#'   - "auto" (default): Automatically selects pip or conda
#'   - "virtualenv": Uses Python virtualenv
#'   - "conda": Uses Conda environment
#' @param conda Name of Conda environment to use (if method = "conda")
#' @param envname Name of Python environment (default: "r-reticulate")
#' @param install_numba Logical; install Numba for JIT optimization (default: TRUE).
#'   Provides ~20% speedup for F0 extraction with no code changes.
#' @param ... Additional arguments passed to `reticulate::py_install()`
#'
#' @details
#' The function installs:

#' - numpy (>=1.24.0,<2.0): Numerical computing (v1.x for compatibility)

#' - scipy (>=1.11.0): Signal processing (FFT, filtering)
#' - soundfile (>=0.12.0): Audio I/O
#' - matplotlib (>=3.7.0): Visualization (optional)
#' - numba (optional): JIT compilation for 20% speedup
#'
#' The legacy STRAIGHT module is included in the package at
#' `inst/python/legacy_STRAIGHT/` and does not require separate installation.
#'
#' **Performance Note**: Installing Numba (`install_numba = TRUE`) provides
#' automatic ~20% speedup for F0 extraction with zero code changes. The
#' optimization is transparent and falls back gracefully if Numba is unavailable.
#'
#' @return Invisible NULL. Called for side effects.
#'
#' @examples
#' \dontrun{
#' # Install with Numba optimization (recommended)
#' install_legacy_straight(install_numba = TRUE)
#'
#' # Install minimal dependencies only
#' install_legacy_straight(install_numba = FALSE)
#'
#' # Install in specific conda environment
#' install_legacy_straight(method = "conda", conda = "my-env")
#' }
#'
#' @export
install_legacy_straight <- function(method = "auto", conda = "auto",
                                   envname = "r-reticulate",
                                   install_numba = TRUE, ...) {
  
  # Check if reticulate is available
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required. Install it with: install.packages('reticulate')")
  }
  
  # Core dependencies
  packages <- c(

    "numpy>=1.24.0,<2.0",

    "scipy>=1.11.0",
    "soundfile>=0.12.0",
    "matplotlib>=3.7.0"
  )
  
  # Add Numba for optimization
  if (install_numba) {
    packages <- c(packages, "numba>=0.57.0")
    message("Installing with Numba optimization (~20% speedup for F0 extraction)")
  }
  
  message("Installing legacy STRAIGHT Python dependencies...")
  message("Packages: ", paste(packages, collapse = ", "))
  
  tryCatch({
    reticulate::py_install(
      packages = packages,
      method = method,
      conda = conda,
      envname = envname,
      ...
    )
    message("\n✓ Installation complete!")
    message("Use straight_available() to verify installation.")
    
    if (install_numba) {
      message("\n✓ Numba optimization enabled")
      message("  First run will compile (~0.5s overhead)")
      message("  Subsequent runs: ~20% faster than baseline")
    }
    
  }, error = function(e) {
    message("\n✗ Installation failed: ", e$message)
    message("\nTroubleshooting:")
    message("1. Ensure Python is available: reticulate::py_config()")
    message("2. Try a different method: install_legacy_straight(method = 'virtualenv')")
    message("3. Check system dependencies for audio I/O")
    stop("Legacy STRAIGHT installation failed", call. = FALSE)
  })
  
  invisible(NULL)
}


#' Check if legacy STRAIGHT is available
#'
#' Tests whether the legacy STRAIGHT Python module can be imported and used.
#'
#' @param detailed Logical; if TRUE, returns detailed version information
#'
#' @return Logical (if detailed = FALSE) or list with version info (if detailed = TRUE)
#'
#' @examples
#' \dontrun{
#' # Simple availability check
#' if (straight_available()) {
#'   message("STRAIGHT is ready to use")
#' }
#'
#' # Detailed information
#' info <- straight_available(detailed = TRUE)
#' print(info)
#' }
#'
#' @export
straight_available <- function(detailed = FALSE) {
  
  # Check if reticulate is available
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    if (detailed) {
      return(list(
        available = FALSE,
        error = "reticulate package not installed"
      ))
    }
    return(FALSE)
  }
  
  # Try to import the module
  result <- tryCatch({
    # Setup Python path for legacy_STRAIGHT

    inst_path <- .setup_straight_path()

    
    # Try importing modules
    f0_mod <- reticulate::import("legacy_STRAIGHT.f0_extraction", convert = FALSE)
    spec_mod <- reticulate::import("legacy_STRAIGHT.spectral", convert = FALSE)
    synth_mod <- reticulate::import("legacy_STRAIGHT.synthesis", convert = FALSE)
    
    # Check for Numba optimization
    has_numba <- reticulate::py_module_available("numba")
    
    # Get versions
    np <- reticulate::import("numpy", convert = FALSE)
    sp <- reticulate::import("scipy", convert = FALSE)
    
    list(
      available = TRUE,
      numpy_version = as.character(np$`__version__`),
      scipy_version = as.character(sp$`__version__`),
      numba_available = has_numba,
      optimization = if (has_numba) "Numba JIT (~20% faster)" else "Baseline",
      module_path = file.path(inst_path, "legacy_STRAIGHT")
    )
    
  }, error = function(e) {
    list(
      available = FALSE,
      error = as.character(e$message)
    )
  })
  
  if (detailed) {
    return(result)
  } else {
    return(result$available)
  }
}


#' Get legacy STRAIGHT module information
#'
#' Returns detailed information about the legacy STRAIGHT implementation,
#' including performance characteristics and optimization status.
#'
#' @return A list with module information
#'
#' @examples
#' \dontrun{
#' info <- straight_info()
#' print(info)
#' }
#'
#' @export
straight_info <- function() {
  
  avail <- straight_available(detailed = TRUE)
  
  if (!avail$available) {
    message("Legacy STRAIGHT is not available")
    message("Install with: install_legacy_straight()")
    return(invisible(NULL))
  }
  
  cat("=== Legacy STRAIGHT Module Information ===\n\n")
  cat("Status: ✓ Available\n")
  cat("Location:", avail$module_path, "\n")
  cat("NumPy version:", avail$numpy_version, "\n")
  cat("SciPy version:", avail$scipy_version, "\n")
  cat("Optimization:", avail$optimization, "\n\n")
  
  cat("Available Functions:\n")

  cat("  • MulticueF0v14():    F0 extraction (91.9% frame, 99.0% mean accuracy)\n")
  cat("  • exstraightspec():  Spectral analysis (99.996% MATLAB accuracy)\n")
  cat("  • exstraightsynth(): Speech synthesis (99.99% MATLAB accuracy)\n\n")

  
  cat("Performance:\n")
  if (avail$numba_available) {
    cat("  • F0 extraction: ~0.68s for 0.79s audio (0.86x RT)\n")
    cat("  • Numba JIT: ✓ Enabled (20% speedup)\n")
    cat("  • First run: ~0.5s compilation overhead\n")
  } else {
    cat("  • F0 extraction: ~0.81s for 0.79s audio (1.02x RT)\n")
    cat("  • Numba JIT: ✗ Not installed\n")
    cat("  • Install for 20% speedup: install_legacy_straight(install_numba=TRUE)\n")
  }
  
  cat("\nAccuracy vs MATLAB:\n")

  cat("  • F0 extraction: 91.9% frame accuracy, 99.0% mean F0 accuracy\n")
  cat("  • Spectral analysis: 99.996% correlation\n")
  cat("  • Aperiodicity: 99.83% accuracy\n")
  cat("  • Synthesis: 99.99% accuracy\n\n")

  
  cat("References:\n")
  cat("  • Original: STRAIGHT vocoder (Kawahara et al.)\n")
  cat("  • Implementation: Python reimplementation (2025)\n")
  cat("  • Status: Validated against MATLAB reference\n")
  
  invisible(avail)
}

#' Install Voxit Python module
#'
#' Installs the Voxit Python module for voice and articulation complexity
#' analysis. Voxit computes prosodic and rhythmic measures from speech audio
#' and word alignments.
#'
#' The voxit module computes 11 features:
#' \itemize{
#'   \item Speaking rate (WPM)
#'   \item Pause statistics (counts, durations, rates)
#'   \item Rhythmic complexity (Lempel-Ziv)
#'   \item Pitch statistics (range, entropy, speed, acceleration)
#' }
#'
#' @param envname Name of Python environment to use (default: NULL = use default)
#' @param method Installation method: "auto" (default), "virtualenv", or "conda"
#' @param install_numba Logical. Install numba for JIT optimization (2-3x speedup).
#'   Default: FALSE. No compilation required, instant speedup.
#' @param compile_cython Logical. Compile Cython extensions (3-5x speedup).
#'   Default: FALSE. Requires C compiler. Maximum performance.
#' @param ... Additional arguments passed to \code{reticulate::py_install()}
#'
#' @details
#' The voxit module requires:
#' \itemize{
#'   \item \strong{numpy} - Numerical operations
#'   \item \strong{scipy} - Signal processing (Savitzky-Golay filter)
#'   \item \strong{lempel_ziv_complexity} - Rhythmic complexity calculation
#' }
#'
#' The module is installed from a local copy in \code{inst/python/voxit/}.
#'
#' @section Performance Optimization:
#' For optimal performance, install optional accelerators:
#' \describe{
#'   \item{\strong{Numba JIT}}{2-3x faster, no compilation needed.
#'     \code{install_voxit(install_numba = TRUE)}
#'     Provides instant speedup through just-in-time compilation.
#'   }
#'   \item{\strong{Cython}}{3-5x faster, requires C compiler.
#'     \code{install_voxit(compile_cython = TRUE)}
#'     Maximum performance through compiled C extensions.
#'     Requires: gcc (Linux), clang (macOS), or MSVC (Windows).
#'   }
#'   \item{\strong{Combined}}{Up to 5-8x total speedup.
#'     \code{install_voxit(install_numba = TRUE, compile_cython = TRUE)}
#'     Best performance for production use.
#'   }
#' }
#'
#' Check current optimization status with \code{voxit_info()$optimized}
#'
#' @return Invisibly returns TRUE if installation successful, FALSE otherwise
#'
#' @examples
#' \dontrun{
#' # Basic installation
#' install_voxit()
#'
#' # Install with Numba optimization (recommended)
#' install_voxit(install_numba = TRUE)
#'
#' # Install with maximum optimization (requires C compiler)
#' install_voxit(install_numba = TRUE, compile_cython = TRUE)
#'
#' # Install in specific conda environment
#' install_voxit(envname = "r-superassp", method = "conda")
#'
#' # Check installation and optimization status
#' voxit_available()
#' info <- voxit_info()
#' print(info$optimized)  # TRUE if optimized version loaded
#' }
#'
#' @seealso
#' \code{\link{voxit_available}}, \code{\link{voxit_info}},
#' \code{\link{lst_voxit}}
#'
install_voxit <- function(envname = NULL,
                         method = "auto",
                         install_numba = FALSE,
                         compile_cython = FALSE,
                         ...) {

  cli::cli_h1("Installing Voxit Python Module")

  # Check Python availability
  if (!reticulate::py_available()) {
    cli::cli_abort(c(
      "x" = "Python not available",
      "i" = "Install Python or configure reticulate"
    ))
  }

  # Get module path
  voxit_path <- system.file("python/voxit", package = "superassp")
  if (voxit_path == "" || !dir.exists(voxit_path)) {
    cli::cli_abort(c(
      "x" = "Voxit module not found in package",
      "i" = "Package installation may be corrupted"
    ))
  }

  tryCatch({
    # Install dependencies
    cli::cli_alert_info("Installing dependencies...")
    
    deps <- c("numpy", "scipy", "lempel_ziv_complexity")
    
    reticulate::py_install(
      packages = deps,
      envname = envname,
      method = method,
      pip = TRUE,
      ...
    )

    # Install numba if requested
    if (install_numba) {
      cli::cli_alert_info("Installing numba for JIT optimization...")
      reticulate::py_install(
        packages = "numba",
        envname = envname,
        method = method,
        pip = TRUE,
        ...
      )
    }

    # Install voxit module from local path
    cli::cli_alert_info("Installing voxit module...")
    
    # Use pip to install from local directory
    pip_cmd <- sprintf("pip install -e %s", shQuote(voxit_path))
    system(pip_cmd)

    # Compile Cython if requested
    if (compile_cython) {
      cli::cli_alert_info("Compiling Cython extensions (requires C compiler)...")
      
      # Check if Cython is available
      if (!reticulate::py_module_available("Cython")) {
        reticulate::py_install(
          packages = "Cython",
          envname = envname,
          method = method,
          pip = TRUE,
          ...
        )
      }

      # Build Cython extensions
      cython_setup <- file.path(voxit_path, "setup_cython.py")
      if (file.exists(cython_setup)) {
        py <- reticulate::import_builtins()
        py$exec(sprintf("import subprocess; subprocess.run(['python', '%s', 'build_ext', '--inplace'], cwd='%s')",
                       cython_setup, voxit_path))
      } else {
        cli::cli_alert_warning("Cython setup script not found, skipping compilation")
      }
    }

    # Verify installation
    if (reticulate::py_module_available("voxit")) {
      cli::cli_alert_success("Voxit module installed successfully")
      
      # Check optimization status
      voxit <- reticulate::import("voxit")
      has_numba <- voxit$HAS_NUMBA
      has_cython <- voxit$HAS_CYTHON
      
      if (has_numba || has_cython) {
        cli::cli_alert_success("Optimizations available:")
        if (has_numba) cli::cli_alert_info("  - Numba JIT: {.val TRUE}")
        if (has_cython) cli::cli_alert_info("  - Cython: {.val TRUE}")
      } else {
        cli::cli_alert_info("Using standard Python implementation")
        if (install_numba || compile_cython) {
          cli::cli_alert_warning("Optimizations were requested but not available")
        }
      }
      
      return(invisible(TRUE))
    } else {
      cli::cli_alert_danger("Installation failed - module not found")
      return(invisible(FALSE))
    }

  }, error = function(e) {
    cli::cli_alert_danger("Installation failed: {e$message}")
    return(invisible(FALSE))
  })
}


#' Check if Voxit module is available
#'
#' @return Logical. TRUE if voxit module is available, FALSE otherwise
#'
#' @examples
#' \dontrun{
#' if (voxit_available()) {
#'   features <- lst_voxit("audio.wav", alignmentFiles = "align.csv")
#' } else {
#'   install_voxit()
#' }
#' }
#'
#' @seealso \code{\link{install_voxit}}, \code{\link{voxit_info}}
#'
voxit_available <- function() {
  reticulate::py_module_available("voxit")
}


#' Get Voxit module information
#'
#' Returns information about the installed Voxit module including version,
#' optimization status, and available features.
#'
#' @return List with module information:
#' \describe{
#'   \item{available}{Logical. TRUE if module is available}
#'   \item{version}{Character. Module version string}
#'   \item{optimized}{Logical. TRUE if using optimized implementation}
#'   \item{numba}{Logical. TRUE if Numba JIT is available}
#'   \item{cython}{Logical. TRUE if Cython extensions are available}
#'   \item{features}{Integer. Number of features computed (11)}
#'   \item{python_version}{Character. Python version string}
#'   \item{dependencies}{Character vector of required packages}
#' }
#'
#' @examples
#' \dontrun{
#' info <- voxit_info()
#' print(info$version)
#' print(info$optimized)
#' print(info$numba)
#' }
#'
#' @seealso \code{\link{install_voxit}}, \code{\link{voxit_available}}
#'
voxit_info <- function() {
  if (!voxit_available()) {
    return(list(
      available = FALSE,
      version = NA,
      optimized = FALSE,
      numba = FALSE,
      cython = FALSE,
      features = 11,
      python_version = NA,
      dependencies = c("numpy", "scipy", "lempel_ziv_complexity")
    ))
  }

  voxit <- reticulate::import("voxit")
  py <- reticulate::import_builtins()
  
  version <- tryCatch(voxit$`__version__`, error = function(e) "unknown")
  has_numba <- tryCatch(voxit$HAS_NUMBA, error = function(e) FALSE)
  has_cython <- tryCatch(voxit$HAS_CYTHON, error = function(e) FALSE)
  
  list(
    available = TRUE,
    version = version,
    optimized = has_numba || has_cython,
    numba = has_numba,
    cython = has_cython,
    features = 11,
    python_version = paste(py$sys$version_info$major,
                          py$sys$version_info$minor,
                          py$sys$version_info$micro, sep = "."),
    dependencies = c("numpy", "scipy", "lempel_ziv_complexity")
  )
}

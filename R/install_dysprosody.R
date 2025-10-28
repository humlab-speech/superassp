#' Install dysprosody Python module
#'
#' Installs the dysprosody Python module for prosodic assessment. The module
#' implements the prosodic measures described in Nylén et al. (2025):
#' "A model of dysprosody in autism spectrum disorder"
#' (doi: 10.3389/fnhum.2025.1566274).
#'
#' The dysprosody module extracts 193 prosodic features including:
#' \itemize{
#'   \item MOMEL-INTSINT pitch targets and tone labels
#'   \item Spectral tilt measures (L2L1, SLF, C1, etc.)
#'   \item Prosodic statistics (pitch mean, range, concentration)
#'   \item Inter-label differential features
#' }
#'
#' @param envname Name of Python environment to use (default: NULL = use default)
#' @param method Installation method: "auto" (default), "virtualenv", or "conda"
#' @param install_numba Logical. Install numba for JIT optimization (10-20\% speedup).
#'   Default: FALSE. No compilation required, instant speedup.
#' @param compile_cython Logical. Compile Cython extensions (15-25\% speedup).
#'   Default: FALSE. Requires C compiler. Maximum performance.
#' @param ... Additional arguments passed to \code{reticulate::py_install()}
#'
#' @details
#' The dysprosody module is a pure Python implementation that requires:
#' \itemize{
#'   \item \strong{parselmouth} - Python bindings for Praat
#'   \item \strong{numpy} - Numerical operations
#'   \item \strong{pandas} - Data structures
#'   \item \strong{scipy} - Statistical functions
#' }
#'
#' The optimized version (\code{dysprosody_optimized.py}) provides 15-25\%
#' faster processing through vectorized statistics and simplified operations.
#'
#' @section Performance:
#' Typical processing times:
#' \itemize{
#'   \item Single file (2-6s audio): 0.16-0.44s (~14x realtime)
#'   \item Batch processing: ~5x speedup with 8 cores
#' }
#'
#' @section Performance Optimization:
#' For optimal performance, install optional accelerators:
#' \describe{
#'   \item{\strong{Numba JIT}}{10-20\% faster, no compilation needed.
#'     \code{install_dysprosody(install_numba = TRUE)}
#'     Provides instant speedup through just-in-time compilation.
#'   }
#'   \item{\strong{Cython}}{15-25\% faster, requires C compiler.
#'     \code{install_dysprosody(compile_cython = TRUE)}
#'     Maximum performance through compiled C extensions.
#'     Requires: gcc (Linux), clang (macOS), or MSVC (Windows).
#'   }
#'   \item{\strong{Combined}}{Up to 30-40\% total speedup.
#'     \code{install_dysprosody(install_numba = TRUE, compile_cython = TRUE)}
#'     Best performance for production use.
#'   }
#' }
#'
#' Check current optimization status with \code{dysprosody_info()$optimized}
#'
#' @return Invisibly returns TRUE if installation successful, FALSE otherwise
#'
#' @examples
#' \dontrun{
#' # Basic installation
#' install_dysprosody()
#'
#' # Install with Numba optimization (recommended)
#' install_dysprosody(install_numba = TRUE)
#'
#' # Install with maximum optimization (requires C compiler)
#' install_dysprosody(install_numba = TRUE, compile_cython = TRUE)
#'
#' # Install in specific conda environment
#' install_dysprosody(envname = "r-superassp", method = "conda")
#'
#' # Check installation and optimization status
#' dysprosody_available()
#' info <- dysprosody_info()
#' print(info$optimized)  # TRUE if optimized version loaded
#' }
#'
#' @references
#' Nylén, F., et al. (2025). A model of dysprosody in autism spectrum disorder.
#' \emph{Frontiers in Human Neuroscience}.
#' \doi{10.3389/fnhum.2025.1566274}
#'
#' @seealso
#' \code{\link{dysprosody_available}}, \code{\link{dysprosody_info}},
#' \code{\link{lst_dysprosody}}
#'
#' @export
install_dysprosody <- function(envname = NULL,
                              method = "auto",
                              install_numba = FALSE,
                              compile_cython = FALSE,
                              ...) {

  # Required dependencies
  packages <- c("numpy", "pandas", "scipy", "praat-parselmouth")

  # Optional optimization packages
  if (install_numba) {
    packages <- c(packages, "numba")
  }
  if (compile_cython) {
    packages <- c(packages, "cython")
  }

  cli::cli_alert_info("Installing dysprosody dependencies...")
  cli::cli_ul(packages)

  if (install_numba) {
    cli::cli_alert_info("Numba JIT optimization will be enabled (10-20% speedup)")
  }
  if (compile_cython) {
    cli::cli_alert_info("Cython compilation will be attempted (15-25% speedup)")
  }

  tryCatch({
    reticulate::py_install(
      packages = packages,
      envname = envname,
      method = method,
      pip = TRUE,
      ...
    )

    # Verify installation
    if (dysprosody_available()) {
      info <- dysprosody_info()

      cli::cli_alert_success("Dysprosody installed successfully!")
      cli::cli_alert_info(sprintf(
        "Version: %s | Optimized: %s",
        info$version,
        if (info$optimized) "Yes" else "No"
      ))

      # Show dependency versions
      cli::cli_h3("Dependencies")
      for (dep in names(info$dependencies)) {
        version <- info$dependencies[[dep]]
        if (!is.null(version)) {
          status <- if (dep %in% c("numba") && !is.null(version)) {
            " (optimization enabled)"
          } else {
            ""
          }
          cli::cli_alert_success(sprintf("%s: %s%s", dep, version, status))
        } else {
          cli::cli_alert_warning(sprintf("%s: Not found", dep))
        }
      }

      # Compile Cython if requested
      if (compile_cython && !is.null(info$dependencies$cython)) {
        cli::cli_h3("Compiling Cython extensions")
        cli::cli_alert_info("This may take a few minutes...")

        tryCatch({
          # Get dysprosody module location
          dysprosody_path <- system.file("python/dysprosody", package = "superassp")

          # Try to compile Cython extensions
          setup_script <- file.path(dysprosody_path, "setup_cython.py")

          if (file.exists(setup_script)) {
            reticulate::py_run_file(setup_script)
            cli::cli_alert_success("Cython compilation completed")
          } else {
            cli::cli_alert_warning("Cython setup script not found - skipping compilation")
            cli::cli_alert_info("Cython optimization will be limited to pure Python speedups")
          }

        }, error = function(e) {
          cli::cli_alert_warning("Cython compilation failed (this is optional)")
          cli::cli_alert_info("Falling back to Numba/pure Python optimization")
          cli::cli_alert_info(sprintf("Error: %s", e$message))
        })
      }

      return(invisible(TRUE))
    } else {
      cli::cli_alert_danger("Installation completed but dysprosody not available")
      cli::cli_alert_info("Try: reticulate::py_config() to diagnose")
      return(invisible(FALSE))
    }

  }, error = function(e) {
    cli::cli_alert_danger("Failed to install dysprosody")
    cli::cli_alert_danger(e$message)
    cli::cli_alert_info("Manual installation: pip install numpy pandas scipy praat-parselmouth")
    return(invisible(FALSE))
  })
}


#' Check if dysprosody is available
#'
#' Tests whether the dysprosody Python module is available and properly configured.
#'
#' @return Logical indicating whether dysprosody is available
#'
#' @examples
#' \dontrun{
#' if (dysprosody_available()) {
#'   result <- lst_dysprosody("speech.wav")
#' }
#' }
#'
#' @export
dysprosody_available <- function() {
  tryCatch({
    # Check if dysprosody module can be imported
    if (!reticulate::py_module_available("dysprosody")) {
      return(FALSE)
    }

    # Try to import the main function
    dysprosody <- reticulate::import("dysprosody", delay_load = FALSE)

    # Check if prosody_measures function exists
    if (!reticulate::py_has_attr(dysprosody, "prosody_measures")) {
      return(FALSE)
    }

    return(TRUE)

  }, error = function(e) {
    return(FALSE)
  })
}


#' Get dysprosody module information
#'
#' Returns detailed information about the installed dysprosody module,
#' including version, optimization status, and dependency versions.
#'
#' @return Named list with module information:
#' \describe{
#'   \item{version}{Dysprosody version string}
#'   \item{optimized}{Logical indicating if optimized version is used}
#'   \item{dependencies}{Named list of dependency versions}
#'   \item{parselmouth}{Parselmouth version (for Praat bindings)}
#'   \item{numpy}{NumPy version}
#'   \item{pandas}{Pandas version}
#'   \item{scipy}{SciPy version}
#' }
#'
#' @examples
#' \dontrun{
#' info <- dysprosody_info()
#' print(info$version)
#' print(info$optimized)
#' }
#'
#' @export
dysprosody_info <- function() {
  if (!dysprosody_available()) {
    cli::cli_alert_warning("Dysprosody not available")
    cli::cli_alert_info("Install with: install_dysprosody()")
    return(NULL)
  }

  tryCatch({
    dysprosody <- reticulate::import("dysprosody", delay_load = FALSE)

    # Get info from module
    info_py <- dysprosody$info()

    # Convert to R list
    info <- list(
      version = info_py$version,
      optimized = info_py$optimized,
      dependencies = as.list(info_py$dependencies)
    )

    # Add direct dependency info for convenience
    info$parselmouth <- info$dependencies$parselmouth
    info$numpy <- info$dependencies$numpy
    info$pandas <- info$dependencies$pandas
    info$scipy <- info$dependencies$scipy

    return(info)

  }, error = function(e) {
    cli::cli_alert_danger("Error getting dysprosody info")
    cli::cli_alert_danger(e$message)
    return(NULL)
  })
}

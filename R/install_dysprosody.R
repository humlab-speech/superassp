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
#' @return Invisibly returns TRUE if installation successful, FALSE otherwise
#'
#' @examples
#' \dontrun{
#' # Install with auto-detection
#' install_dysprosody()
#'
#' # Install in specific conda environment
#' install_dysprosody(envname = "r-superassp", method = "conda")
#'
#' # Check installation
#' dysprosody_available()
#' dysprosody_info()
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
install_dysprosody <- function(envname = NULL, method = "auto", ...) {

  # Required dependencies
  packages <- c("numpy", "pandas", "scipy", "praat-parselmouth")

  cli::cli_alert_info("Installing dysprosody dependencies...")
  cli::cli_ul(packages)

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
          cli::cli_alert_success(sprintf("%s: %s", dep, version))
        } else {
          cli::cli_alert_warning(sprintf("%s: Not found", dep))
        }
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

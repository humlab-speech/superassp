#' Install SAcC (Subband Autocorrelation Classification) pitch tracker
#'
#' This function installs the required Python dependencies for the SAcC pitch tracker.
#' SAcC is a robust pitch tracking algorithm based on subband autocorrelation analysis
#' and neural network classification.
#'
#' @param envname The name of the Python environment to use. If NULL (default),
#'   uses the default reticulate environment.
#' @param method Installation method to use. Options are "auto" (default),
#'   "virtualenv", or "conda".
#' @param ... Additional arguments passed to \code{reticulate::py_install()}.
#'
#' @details
#' The SAcC pitch tracker requires the following Python packages:
#' \itemize{
#'   \item numpy - Numerical computing
#'   \item scipy - Scientific computing (signal processing, filters)
#'   \item soundfile - Audio file I/O (for SPH format support)
#' }
#'
#' SAcC uses pre-trained neural network weights and PCA mappings that are
#' included with the superassp package in \code{inst/python/calc_sbpca/python/aux/}.
#'
#' The algorithm processes audio at 8kHz (automatic downsampling) and returns
#' pitch tracks with 10ms frame shifts.
#'
#' @return Invisible NULL. Called for side effects.
#'
#' @examples
#' \dontrun{
#' # Install SAcC dependencies
#' install_sacc()
#'
#' # Install to a specific conda environment
#' install_sacc(envname = "r-superassp", method = "conda")
#'
#' # Check if SAcC is available after installation
#' sacc_available()
#' }
#'
#' @seealso
#' \code{\link{sacc_available}} to check if SAcC is available
#' \code{\link{sacc_info}} to get SAcC configuration information
#' \code{\link{trk_sacc}} to run SAcC pitch tracking
#'
#' @export
install_sacc <- function(envname = NULL, method = "auto", ...) {
  # Required packages for SAcC
  packages <- c("numpy", "scipy", "soundfile")

  message("Installing Python dependencies for SAcC pitch tracker...")
  message("  - numpy (numerical computing)")
  message("  - scipy (signal processing)")
  message("  - soundfile (audio I/O)")

  reticulate::py_install(packages, envname = envname, method = method, ...)

  message("\nSAcC installation complete!")
  message("Test with: sacc_available()")

  invisible(NULL)
}

#' Check if SAcC pitch tracker is available
#'
#' Checks whether the required Python modules for SAcC are available.
#'
#' @return Logical. TRUE if all required modules are available, FALSE otherwise.
#'
#' @examples
#' \dontrun{
#' if (sacc_available()) {
#'   result <- trk_sacc("speech.wav", toFile = FALSE)
#' } else {
#'   install_sacc()
#' }
#' }
#'
#' @seealso
#' \code{\link{install_sacc}} to install SAcC dependencies
#' \code{\link{sacc_info}} to get SAcC configuration information
#'
#' @export
sacc_available <- function() {
  numpy_available <- reticulate::py_module_available("numpy")
  scipy_available <- reticulate::py_module_available("scipy")
  soundfile_available <- reticulate::py_module_available("soundfile")

  all(numpy_available, scipy_available, soundfile_available)
}

#' Get SAcC pitch tracker information
#'
#' Returns information about the SAcC installation and configuration.
#'
#' @return A list with information about:
#' \itemize{
#'   \item Python path and version
#'   \item NumPy version
#'   \item SciPy version
#'   \item soundfile availability
#'   \item SAcC module path
#'   \item Available configuration files
#' }
#'
#' @examples
#' \dontrun{
#' if (sacc_available()) {
#'   info <- sacc_info()
#'   print(info)
#' }
#' }
#'
#' @seealso
#' \code{\link{install_sacc}} to install SAcC dependencies
#' \code{\link{sacc_available}} to check availability
#'
#' @export
sacc_info <- function() {
  if (!sacc_available()) {
    message("SAcC is not available. Install with: install_sacc()")
    return(invisible(NULL))
  }

  # Get Python configuration
  py_config <- reticulate::py_config()

  # Get module versions
  np <- reticulate::import("numpy")
  scipy <- reticulate::import("scipy")

  # Get SAcC module path
  sacc_path <- system.file("python", "calc_sbpca", "python", package = "superassp")
  aux_path <- file.path(sacc_path, "aux")

  # List available configuration files
  config_files <- list.files(aux_path, pattern = "\\.(mat|wgt|norms|txt)$")

  info <- list(
    python_path = py_config$python,
    python_version = py_config$version,
    numpy_version = np$`__version__`,
    scipy_version = scipy$`__version__`,
    soundfile_available = reticulate::py_module_available("soundfile"),
    sacc_module_path = sacc_path,
    aux_files_path = aux_path,
    config_files = config_files,
    description = "SAcC (Subband Autocorrelation Classification) pitch tracker by Dan Ellis"
  )

  class(info) <- c("sacc_info", "list")
  return(info)
}

#' @export
print.sacc_info <- function(x, ...) {
  cat("SAcC Pitch Tracker Configuration\n")
  cat("=================================\n\n")
  cat("Python Environment:\n")
  cat("  Path:", x$python_path, "\n")
  cat("  Version:", x$python_version, "\n\n")
  cat("Required Modules:\n")
  cat("  NumPy:", x$numpy_version, "\n")
  cat("  SciPy:", x$scipy_version, "\n")
  cat("  soundfile:", ifelse(x$soundfile_available, "available", "NOT AVAILABLE"), "\n\n")
  cat("SAcC Module:\n")
  cat("  Location:", x$sacc_module_path, "\n")
  cat("  Aux files:", x$aux_files_path, "\n")
  cat("  Config files:", length(x$config_files), "files\n\n")
  cat("Description:\n")
  cat(" ", x$description, "\n")

  invisible(x)
}

#' Install DeepFormants dependencies
#'
#' This function installs the required Python dependencies for DeepFormants formant
#' tracking and estimation. DeepFormants uses deep neural networks trained on labeled
#' formant data to provide accurate F1-F4 estimates.
#'
#' @param envname The name of the Python environment to use. If NULL (default),
#'   uses the default reticulate environment.
#' @param method Installation method to use. Options are "auto" (default),
#'   "virtualenv", or "conda".
#' @param ... Additional arguments passed to \code{reticulate::py_install()}.
#'
#' @details
#' DeepFormants requires the following Python packages:
#' \itemize{
#'   \item torch - PyTorch deep learning framework
#'   \item numpy - Numerical computing
#'   \item scipy - Scientific computing (LPC analysis, signal processing)
#'   \item pandas - Data manipulation
#'   \item numba - JIT compilation for 2-3x performance boost
#' }
#'
#' The package includes pre-trained models:
#' \itemize{
#'   \item Estimation model (estimation_model.dat) - Single time-window estimates
#'   \item Tracking model (LPC_NN.pt) - Continuous 10ms tracking
#' }
#'
#' @return Invisible NULL. Called for side effects.
#'
#' @examples
#' \dontrun{
#' # Install DeepFormants dependencies
#' install_deepformants()
#'
#' # Install to a specific conda environment
#' install_deepformants(envname = "r-superassp", method = "conda")
#'
#' # Check if DeepFormants is available after installation
#' deepformants_available()
#' }
#'
#' @seealso
#' \code{\link{deepformants_available}} to check if DeepFormants is available
#' \code{\link{deepformants_info}} to get DeepFormants configuration information
#' \code{\link{trk_deepformants}} to track formants
#' \code{\link{lst_deepformants}} to estimate formants
#'
#' @export
install_deepformants <- function(envname = NULL, method = "auto", ...) {
  # Required packages for DeepFormants
  packages <- c("torch", "numpy", "scipy", "pandas", "numba")

  message("Installing Python dependencies for DeepFormants...")
  message("  - torch (PyTorch deep learning framework)")
  message("  - numpy (numerical computing)")
  message("  - scipy (LPC analysis, signal processing)")
  message("  - pandas (data manipulation)")
  message("  - numba (JIT compilation for performance)")

  reticulate::py_install(packages, envname = envname, method = method, ...)

  message("\nDeepFormants installation complete!")
  message("Test with: deepformants_available()")
  message("\nPre-trained models included:")
  message("  - Estimation: estimation_model.dat")
  message("  - Tracking: LPC_NN.pt")

  invisible(NULL)
}

#' Check if DeepFormants is available
#'
#' Checks whether the required Python modules for DeepFormants are available.
#'
#' @return Logical. TRUE if all required modules are available, FALSE otherwise.
#'
#' @examples
#' \dontrun{
#' if (deepformants_available()) {
#'   result <- trk_deepformants("speech.wav", toFile = FALSE)
#' } else {
#'   install_deepformants()
#' }
#' }
#'
#' @seealso
#' \code{\link{install_deepformants}} to install DeepFormants dependencies
#' \code{\link{deepformants_info}} to get DeepFormants configuration information
#'
#' @export
deepformants_available <- function() {
  torch_available <- reticulate::py_module_available("torch")
  numpy_available <- reticulate::py_module_available("numpy")
  scipy_available <- reticulate::py_module_available("scipy")
  pandas_available <- reticulate::py_module_available("pandas")
  numba_available <- reticulate::py_module_available("numba")

  all(torch_available, numpy_available, scipy_available, pandas_available, numba_available)
}

#' Get DeepFormants information
#'
#' Returns information about the DeepFormants installation and configuration.
#'
#' @return A list with information about:
#' \itemize{
#'   \item Python path and version
#'   \item PyTorch version
#'   \item NumPy version
#'   \item SciPy version
#'   \item Pandas version
#'   \item Numba version
#'   \item DeepFormants module path
#'   \item Available model files
#' }
#'
#' @examples
#' \dontrun{
#' if (deepformants_available()) {
#'   info <- deepformants_info()
#'   print(info)
#' }
#' }
#'
#' @seealso
#' \code{\link{install_deepformants}} to install DeepFormants dependencies
#' \code{\link{deepformants_available}} to check availability
#'
#' @export
deepformants_info <- function() {
  if (!deepformants_available()) {
    message("DeepFormants is not available. Install with: install_deepformants()")
    return(invisible(NULL))
  }

  # Get Python configuration
  py_config <- reticulate::py_config()

  # Get module versions
  torch <- reticulate::import("torch")
  np <- reticulate::import("numpy")
  scipy <- reticulate::import("scipy")
  pandas <- reticulate::import("pandas")
  numba <- reticulate::import("numba")

  # Get DeepFormants module path
  df_path <- system.file("python", "DeepFormants", package = "superassp")

  # List available model files
  model_files <- list.files(df_path, pattern = "\\.(dat|pt|pth)$", recursive = TRUE)

  info <- list(
    python_path = py_config$python,
    python_version = py_config$version,
    torch_version = torch$`__version__`,
    numpy_version = np$`__version__`,
    scipy_version = scipy$`__version__`,
    pandas_version = pandas$`__version__`,
    numba_version = numba$`__version__`,
    deepformants_path = df_path,
    model_files = model_files,
    description = "DeepFormants: Deep learning formant tracking and estimation by Dissen & Keshet"
  )

  class(info) <- c("deepformants_info", "list")
  return(info)
}

#' @export
print.deepformants_info <- function(x, ...) {
  cat("DeepFormants Configuration\n")
  cat("==========================\n\n")
  cat("Python Environment:\n")
  cat("  Path:", x$python_path, "\n")
  cat("  Version:", x$python_version, "\n\n")
  cat("Required Modules:\n")
  cat("  PyTorch:", x$torch_version, "\n")
  cat("  NumPy:", x$numpy_version, "\n")
  cat("  SciPy:", x$scipy_version, "\n")
  cat("  Pandas:", x$pandas_version, "\n")
  cat("  Numba:", x$numba_version, "\n\n")
  cat("DeepFormants Module:\n")
  cat("  Location:", x$deepformants_path, "\n")
  cat("  Model files:", length(x$model_files), "files\n")
  if (length(x$model_files) > 0) {
    cat("    -", paste(x$model_files, collapse = "\n    - "), "\n")
  }
  cat("\nDescription:\n")
  cat(" ", x$description, "\n")

  invisible(x)
}

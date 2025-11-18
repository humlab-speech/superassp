#' Install Phonet dependencies
#'
#' This function installs the required Python dependencies for Phonet phonological
#' posterior extraction. Phonet \insertCite{vasquez2019phonet}{superassp} uses
#' Bidirectional Gated Recurrent Unit (BGRU) neural networks to compute posterior
#' probabilities of phonological classes from speech audio.
#'
#' @param envname The name of the Python environment to use. If NULL (default),
#'   uses the default reticulate environment.
#' @param method Installation method to use. Options are "auto" (default),
#'   "virtualenv", or "conda".
#' @param ... Additional arguments passed to \code{reticulate::py_install()}.
#'
#' @details
#' Phonet requires the following Python packages:
#' \itemize{
#'   \item tensorflow - Deep learning framework for BGRU models (TensorFlow 2.x)
#'   \item keras - High-level neural network API (2.3.1)
#'   \item pandas - Data manipulation for posteriors
#'   \item pysptk - Speech signal processing toolkit
#'   \item python_speech_features - Mel-filterbank feature extraction
#'   \item numpy - Numerical computing
#'   \item scipy - Signal processing (resampling, filtering)
#'   \item scikit-learn - Machine learning utilities
#'   \item numba - JIT compilation for performance
#'   \item tqdm - Progress bars for batch processing
#' }
#'
#' The package includes pre-trained models for:
#' \itemize{
#'   \item Multi-task phonological classifier (17 binary classifiers)
#'   \item Phoneme recognition model
#'   \item Feature normalization parameters (mu, sigma)
#' }
#'
#' \strong{Phonological Classes Available:}
#' vocalic, consonantal, back, anterior, open, close, nasal, stop, continuant,
#' lateral, flap, trill, voice, strident, labial, dental, velar, pause
#'
#' @return Invisible NULL. Called for side effects.
#'
#' @examples
#' \dontrun{
#' # Install Phonet dependencies
#' install_phonet()
#'
#' # Install to a specific conda environment
#' install_phonet(envname = "r-superassp", method = "conda")
#'
#' # Check if Phonet is available after installation
#' phonet_available()
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso
#' \code{\link{phonet_available}} to check if Phonet is available
#' \code{\link{phonet_info}} to get Phonet configuration information
#' \code{\link{lst_phonet}} to extract phonological posteriors
#'
#' @export
install_phonet <- function(envname = NULL, method = "auto", ...) {
  # Required packages for Phonet
  packages <- c(
    "tensorflow",
    "tf-keras",  # Keras 2.x compatibility for Python 3.12+
    "pandas",
    "pysptk",
    "python_speech_features",
    "numpy",
    "scipy",
    "scikit-learn",
    "numba",
    "tqdm",
    "six"
  )

  message("Installing Python dependencies for Phonet...")
  message("  - tensorflow (TensorFlow 2.x deep learning framework)")
  message("  - tf-keras (Keras 2.x API compatibility)")
  message("  - pandas (data manipulation)")
  message("  - pysptk (speech signal processing)")
  message("  - python_speech_features (Mel-filterbank extraction)")
  message("  - numpy (numerical computing)")
  message("  - scipy (signal processing)")
  message("  - scikit-learn (machine learning utilities)")
  message("  - numba (JIT compilation for performance)")
  message("  - tqdm (progress bars)")

  reticulate::py_install(packages, envname = envname, method = method, ...)

  # Install phonet from GitHub
  message("\nInstalling Phonet from GitHub...")
  reticulate::py_install(
    "git+https://github.com/jcvasquezc/phonet.git",
    envname = envname,
    method = method,
    ...
  )

  message("\nPhonet installation complete!")
  message("Test with: phonet_available()")
  message("\nPre-trained models included:")
  message("  - Multi-task phonological classifier (model.h5)")
  message("  - Phoneme recognition model (phonemes.hdf5)")
  message("  - Feature normalization (mu.npy, std.npy)")
  message("\n17 phonological classes available:")
  message("  vocalic, consonantal, back, anterior, open, close, nasal,")
  message("  stop, continuant, lateral, flap, trill, voice, strident,")
  message("  labial, dental, velar, pause")

  invisible(NULL)
}

#' Check if Phonet is available
#'
#' Checks whether the required Python modules for Phonet are available.
#'
#' @return Logical. TRUE if all required modules are available, FALSE otherwise.
#'
#' @examples
#' \dontrun{
#' if (phonet_available()) {
#'   result <- lst_phonet("speech.wav", classes = c("nasal", "stop"))
#' } else {
#'   install_phonet()
#' }
#' }
#'
#' @seealso
#' \code{\link{install_phonet}} to install Phonet dependencies
#' \code{\link{phonet_info}} to get Phonet configuration information
#'
#' @export
phonet_available <- function() {
  phonet_available <- reticulate::py_module_available("phonet")
  tensorflow_available <- reticulate::py_module_available("tensorflow")
  keras_available <- reticulate::py_module_available("keras")
  pysptk_available <- reticulate::py_module_available("pysptk")
  python_speech_features_available <- reticulate::py_module_available("python_speech_features")
  numba_available <- reticulate::py_module_available("numba")

  all(
    phonet_available,
    tensorflow_available,
    keras_available,
    pysptk_available,
    python_speech_features_available,
    numba_available
  )
}

#' Get Phonet information
#'
#' Returns information about the Phonet installation and configuration.
#'
#' @return A list with information about:
#' \itemize{
#'   \item Python path and version
#'   \item TensorFlow version
#'   \item Keras version
#'   \item NumPy version
#'   \item SciPy version
#'   \item Numba version
#'   \item Phonet version and location
#'   \item Available phonological classes
#' }
#'
#' @examples
#' \dontrun{
#' if (phonet_available()) {
#'   info <- phonet_info()
#'   print(info)
#' }
#' }
#'
#' @seealso
#' \code{\link{install_phonet}} to install Phonet dependencies
#' \code{\link{phonet_available}} to check availability
#'
#' @export
phonet_info <- function() {
  if (!phonet_available()) {
    message("Phonet is not available. Install with: install_phonet()")
    return(invisible(NULL))
  }

  # Get Python configuration
  py_config <- reticulate::py_config()

  # Get module versions
  tryCatch({
    tensorflow <- reticulate::import("tensorflow")
    keras <- reticulate::import("keras")
    np <- reticulate::import("numpy")
    scipy <- reticulate::import("scipy")
    numba <- reticulate::import("numba")
    phonet <- reticulate::import("phonet")

    # Get phonological classes
    Phon <- phonet$Phonological()
    phon_classes <- Phon$get_list_phonological_keys()

    info <- list(
      python_path = py_config$python,
      python_version = py_config$version,
      tensorflow_version = tensorflow$`__version__`,
      keras_version = keras$`__version__`,
      numpy_version = np$`__version__`,
      scipy_version = scipy$`__version__`,
      numba_version = numba$`__version__`,
      phonet_version = if (!is.null(phonet$`__version__`)) phonet$`__version__` else "0.3.7",
      phonological_classes = phon_classes,
      num_classes = length(phon_classes),
      description = "Phonet: BGRU-based phonological posterior extraction by Vásquez-Correa et al."
    )

    class(info) <- c("phonet_info", "list")
    return(info)
  }, error = function(e) {
    warning("Error retrieving Phonet info: ", e$message)
    return(invisible(NULL))
  })
}

#' @export
print.phonet_info <- function(x, ...) {
  cat("Phonet Configuration\n")
  cat("====================\n\n")
  cat("Python Environment:\n")
  cat("  Path:", x$python_path, "\n")
  cat("  Version:", x$python_version, "\n\n")
  cat("Required Modules:\n")
  cat("  TensorFlow:", x$tensorflow_version, "\n")
  cat("  Keras:", x$keras_version, "\n")
  cat("  NumPy:", x$numpy_version, "\n")
  cat("  SciPy:", x$scipy_version, "\n")
  cat("  Numba:", x$numba_version, "\n\n")
  cat("Phonet Module:\n")
  cat("  Version:", x$phonet_version, "\n")
  cat("  Phonological classes:", x$num_classes, "\n")
  if (length(x$phonological_classes) > 0) {
    cat("    ", paste(x$phonological_classes, collapse = ", "), "\n")
  }
  cat("\nDescription:\n")
  cat(" ", x$description, "\n")

  invisible(x)
}

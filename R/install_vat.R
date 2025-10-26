# ==============================================================================
# Voice Analysis Toolkit - Installation and Availability Functions
# ==============================================================================

#' Install Voice Analysis Toolkit Python Dependencies
#'
#' Installs required Python packages for the Voice Analysis Toolkit.
#'
#' @param method Installation method: "auto", "virtualenv", or "conda"
#' @param conda Conda environment specification
#' @param pip Logical; use pip for installation (default: TRUE)
#'
#' @return Invisible TRUE on success
#'
#' @details
#' The Voice Analysis Toolkit requires:
#' \itemize{
#'   \item numpy >= 1.21.0
#'   \item scipy >= 1.7.0
#'   \item soundfile >= 0.11.0
#'   \item pywavelets >= 1.3.0
#' }
#'
#' Optional packages for enhanced functionality:
#' \itemize{
#'   \item matplotlib - For visualization
#'   \item librosa - Additional audio processing
#' }
#'
#' @examples
#' \dontrun{
#' # Install with default settings
#' install_vat()
#'
#' # Install in specific conda environment
#' install_vat(method = "conda", conda = "r-reticulate")
#' }
#'
#' @export
install_vat <- function(method = "auto", conda = "auto", pip = TRUE) {

  cli::cli_h1("Installing Voice Analysis Toolkit Dependencies")

  # Required packages
  required_pkgs <- c(
    "numpy>=1.21.0",
    "scipy>=1.7.0",
    "soundfile>=0.11.0",
    "pywavelets>=1.3.0"
  )

  cli::cli_alert_info("Installing required packages: {.pkg {required_pkgs}}")

  tryCatch({
    reticulate::py_install(
      packages = required_pkgs,
      method = method,
      conda = conda,
      pip = pip
    )

    cli::cli_alert_success("Voice Analysis Toolkit dependencies installed successfully")

    # Check if installation worked
    if (vat_available()) {
      cli::cli_alert_success("Voice Analysis Toolkit is now available")

      # Show info
      info <- vat_info()
      cli::cli_alert_info("Python location: {.path {info$python_path}}")

      return(invisible(TRUE))
    } else {
      cli::cli_alert_warning("Installation completed but VAT not detected")
      cli::cli_alert_info("Try restarting R session")
      return(invisible(FALSE))
    }

  }, error = function(e) {
    cli::cli_alert_danger("Installation failed: {e$message}")
    cli::cli_alert_info("Try manual installation:")
    cli::cli_code("pip install numpy scipy soundfile pywavelets")
    return(invisible(FALSE))
  })
}

#' Check Voice Analysis Toolkit Availability
#'
#' Check if the Voice Analysis Toolkit Python module is available.
#'
#' @return Logical; TRUE if VAT is available, FALSE otherwise
#'
#' @examples
#' if (vat_available()) {
#'   # Use VAT functions
#' } else {
#'   install_vat()
#' }
#'
#' @export
vat_available <- function() {
  tryCatch({
    # Check Python modules
    has_numpy <- reticulate::py_module_available("numpy")
    has_scipy <- reticulate::py_module_available("scipy")

    # Check VAT package location
    vat_path <- system.file("python", "voice_analysis_toolkit",
                           package = "superassp")
    has_vat <- file.exists(vat_path) && file.exists(file.path(vat_path, "__init__.py"))

    all(has_numpy, has_scipy, has_vat)

  }, error = function(e) {
    FALSE
  })
}

#' Get Voice Analysis Toolkit Information
#'
#' Retrieve information about the installed Voice Analysis Toolkit.
#'
#' @return List with VAT installation details:
#' \itemize{
#'   \item \code{available}: Logical; is VAT available
#'   \item \code{python_path}: Path to Python interpreter
#'   \item \code{vat_path}: Path to VAT Python package
#'   \item \code{numpy_version}: NumPy version
#'   \item \code{scipy_version}: SciPy version
#'   \item \code{modules}: Available VAT modules
#' }
#'
#' @examples
#' \dontrun{
#' info <- vat_info()
#' print(info)
#' }
#'
#' @export
vat_info <- function() {
  if (!vat_available()) {
    return(list(
      available = FALSE,
      message = "Voice Analysis Toolkit not available. Install with install_vat()"
    ))
  }

  # Get Python info
  py_config <- reticulate::py_config()

  # Get module versions
  numpy_version <- tryCatch({
    np <- reticulate::import("numpy")
    np$`__version__`
  }, error = function(e) "not available")

  scipy_version <- tryCatch({
    sp <- reticulate::import("scipy")
    sp$`__version__`
  }, error = function(e) "not available")

  # Get VAT location
  vat_path <- system.file("python", "voice_analysis_toolkit",
                         package = "superassp")

  # Check available modules
  modules <- c(
    "general" = file.exists(file.path(vat_path, "general", "__init__.py")),
    "se_vq" = file.exists(file.path(vat_path, "se_vq", "__init__.py")),
    "creak" = file.exists(file.path(vat_path, "creak", "__init__.py")),
    "utils" = file.exists(file.path(vat_path, "utils", "__init__.py"))
  )

  list(
    available = TRUE,
    python_path = py_config$python,
    python_version = py_config$version,
    vat_path = vat_path,
    numpy_version = numpy_version,
    scipy_version = scipy_version,
    modules = modules,
    all_modules_available = all(modules)
  )
}

# Internal function to get VAT module
get_vat_module <- function() {
  if (!vat_available()) {
    stop("Voice Analysis Toolkit not available.\n",
         "Install with: install_vat()\n",
         "Check status with: vat_info()",
         call. = FALSE)
  }

  vat_path <- system.file("python", package = "superassp")

  tryCatch({
    reticulate::import_from_path(
      "voice_analysis_toolkit",
      path = vat_path
    )
  }, error = function(e) {
    stop("Failed to import Voice Analysis Toolkit: ", e$message,
         call. = FALSE)
  })
}

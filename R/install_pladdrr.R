#' Install pladdrr R Package
#'
#' Helper function to install the pladdrr package, which provides R bindings
#' to Praat's phonetic analysis library. This is required for all Praat-based
#' DSP functions in superassp (pitch, formants, intensity, voice quality, etc.).
#'
#' @param method Installation method: "auto", "cran", or "github"
#' @param version Specific version to install (optional)
#' @param ... Additional arguments passed to install.packages() or remotes::install_github()
#'
#' @details
#' pladdrr (Praat Library Bindings for R) provides direct R6 object-oriented
#' access to Praat's C phonetic analysis library. It replaces the previous
#' Python-based parselmouth integration, offering:
#' \itemize{
#'   \item Native R integration (no Python dependency)
#'   \item R6 object-oriented API
#'   \item Direct C library access for performance
#'   \item Full Praat feature set
#' }
#'
#' Minimum version required: 4.8.16
#'
#' @references
#' \itemize{
#'   \item pladdrr GitHub: \url{https://github.com/tjmahr/pladdrr}
#'   \item Praat: \url{http://www.praat.org}
#' }
#'
#' @examples
#' \dontrun{
#' # Install from CRAN (when available)
#' install_pladdrr()
#'
#' # Install from GitHub
#' install_pladdrr(method = "github")
#'
#' # Check if installed
#' if (pladdrr_available()) {
#'   message("pladdrr is ready!")
#' }
#' }
#'
install_pladdrr <- function(method = "auto", version = NULL, ...) {
  
  if (method == "auto") {
    # Try CRAN first, fall back to GitHub
    tryCatch({
      if (is.null(version)) {
        utils::install.packages("pladdrr", ...)
      } else {
        remotes::install_version("pladdrr", version = version, ...)
      }
      message("✓ pladdrr installed from CRAN")
    }, error = function(e) {
      message("CRAN installation failed, trying GitHub...")
      install_pladdrr(method = "github", version = version, ...)
    })
    
  } else if (method == "cran") {
    if (is.null(version)) {
      utils::install.packages("pladdrr", ...)
    } else {
      if (!requireNamespace("remotes", quietly = TRUE)) {
        stop("Package 'remotes' required to install specific version. Install it first.")
      }
      remotes::install_version("pladdrr", version = version, ...)
    }
    message("✓ pladdrr installed from CRAN")
    
  } else if (method == "github") {
    if (!requireNamespace("remotes", quietly = TRUE)) {
      stop("Package 'remotes' required for GitHub installation. Install it first.")
    }
    
    repo <- "tjmahr/pladdrr"
    if (!is.null(version)) {
      repo <- paste0(repo, "@v", version)
    }
    
    remotes::install_github(repo, ...)
    message("✓ pladdrr installed from GitHub")
    
  } else {
    stop("Invalid method. Choose 'auto', 'cran', or 'github'")
  }
  
  # Verify installation
  if (!pladdrr_available()) {
    stop("pladdrr installation completed but package not loadable. Try restarting R.")
  }
  
  invisible(TRUE)
}


#' Check if pladdrr is Available
#'
#' Tests whether the pladdrr package is installed and loadable.
#'
#' @return Logical: TRUE if pladdrr is available, FALSE otherwise
#'
#' @examples
#' if (pladdrr_available()) {
#'   message("pladdrr is ready for Praat-based analysis")
#' } else {
#'   message("Install pladdrr with: install_pladdrr()")
#' }
#'
pladdrr_available <- function() {
  requireNamespace("pladdrr", quietly = TRUE)
}


#' Get pladdrr Package Information
#'
#' Returns version and configuration information about the installed
#' pladdrr package.
#'
#' @return List with package information:
#' \itemize{
#'   \item installed: Whether pladdrr is installed
#'   \item version: Package version (if installed)
#'   \item path: Installation path
#'   \item minimum_version: Minimum version required by superassp
#'   \item status: "ok", "outdated", or "not_installed"
#' }
#'
#' @examples
#' pladdrr_info()
#'
pladdrr_info <- function() {
  
  minimum_version <- "4.8.16"
  
  if (!pladdrr_available()) {
    return(list(
      installed = FALSE,
      version = NA,
      path = NA,
      minimum_version = minimum_version,
      status = "not_installed",
      message = "pladdrr not installed. Run install_pladdrr()"
    ))
  }
  
  version <- as.character(utils::packageVersion("pladdrr"))
  path <- system.file(package = "pladdrr")
  
  # Compare versions
  installed_ver <- package_version(version)
  required_ver <- package_version(minimum_version)
  
  status <- if (installed_ver >= required_ver) {
    "ok"
  } else {
    "outdated"
  }
  
  result <- list(
    installed = TRUE,
    version = version,
    path = path,
    minimum_version = minimum_version,
    status = status
  )
  
  if (status == "outdated") {
    result$message <- sprintf(
      "pladdrr %s installed but %s required. Update with install_pladdrr()",
      version, minimum_version
    )
  } else {
    result$message <- sprintf("pladdrr %s is installed and ready", version)
  }
  
  return(result)
}


#' Get pladdrr Specifications
#'
#' Returns detailed technical specifications about pladdrr capabilities.
#' Used internally for validation and feature detection.
#'
#' @return List with specifications:
#' \itemize{
#'   \item api_version: API version string
#'   \item features: Available feature set
#'   \item direct_api: Whether direct C API is available
#'   \item ultra_api: Whether Ultra API optimizations are available
#' }
#'
#' @keywords internal
pladdrr_specs <- function() {
  
  if (!pladdrr_available()) {
    stop("pladdrr not available. Install with install_pladdrr()")
  }
  
  version <- utils::packageVersion("pladdrr")
  
  # Feature detection based on version
  # v4.8.16 includes formant bug fix and Ultra API
  has_ultra_api <- version >= "4.6.4"
  has_formant_fix <- version >= "4.8.16"
  
  list(
    version = as.character(version),
    api_version = "4.x",
    direct_api = TRUE,  # All recent versions have direct API
    ultra_api = has_ultra_api,
    formant_bug_fixed = has_formant_fix,
    features = list(
      pitch = c("cc", "ac"),  # Cross-correlation, autocorrelation
      formants = "burg",       # Burg's method
      formant_tracking = "hmm", # HMM tracking (may have limitations)
      intensity = TRUE,
      spectrogram = TRUE,
      power_cepstrogram = TRUE,
      cpps = has_ultra_api,   # CPPS via Ultra API
      voice_quality = has_ultra_api, # Combined HNR/shimmer
      voiced_extraction = has_ultra_api
    ),
    notes = c(
      if (!has_formant_fix) "Formant values may be systematically low (35-55%). Update to 4.8.16+",
      if (has_ultra_api) "Ultra API available for performance optimization"
    )
  )
}

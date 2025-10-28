# DisVoice Python Module Initialization
#
# Provides lazy loading and initialization of DisVoice Python modules
# for optimized speech analysis using Parselmouth.

.disvoice_env <- NULL

#' Initialize DisVoice Python Environment
#'
#' Lazily loads the DisVoice Python modules using reticulate.
#' The environment is cached after first load for performance.
#'
#' @return DisVoice Python module environment (invisible)
#' @keywords internal
init_disvoice <- function() {
  if (is.null(.disvoice_env)) {
    # Try installed package path first
    python_path <- system.file("python", package = "superassp", mustWork = FALSE)

    # Check if DisVoice exists in installed location
    if (python_path != "") {
      disvoice_path <- file.path(python_path, "DisVoice")
      if (!dir.exists(disvoice_path)) {
        # Installed package doesn't have DisVoice, try development path
        python_path <- ""
      }
    }

    # If not found in installed location, try development mode (inst/python)
    if (python_path == "") {
      python_path <- file.path(getwd(), "inst", "python")
      disvoice_path <- file.path(python_path, "DisVoice")

      if (!dir.exists(disvoice_path)) {
        stop("DisVoice Python modules not found.\n",
             "  Tried: ", disvoice_path, "\n",
             "  Expected structure: inst/python/DisVoice/")
      }
    }

    .disvoice_env <<- reticulate::import_from_path(
      "DisVoice",
      path = python_path,
      convert = FALSE
    )
  }
  invisible(.disvoice_env)
}

#' Check if DisVoice Python Support is Available
#'
#' Tests whether DisVoice Python modules can be loaded successfully.
#' This checks for:
#' - reticulate availability
#' - Python installation
#' - parselmouth package
#' - DisVoice modules
#'
#' @return Logical: TRUE if DisVoice is available, FALSE otherwise
#'
#' @references
#' \insertCite{Jadoul2018}{superassp}
#'
#' \insertCite{OrozcoArroyave2018}{superassp}
#'
#' @export
#' @examples
#' if (has_disvoice_support()) {
#'   # Use DisVoice functions
#' } else {
#'   message("DisVoice not available, using fallback methods")
#' }
has_disvoice_support <- function() {
  tryCatch({
    init_disvoice()
    TRUE
  }, error = function(e) {
    FALSE
  })
}

#' Get DisVoice Environment
#'
#' Returns the cached DisVoice Python environment.
#' Initializes if not already loaded.
#'
#' @return DisVoice Python module environment
#' @keywords internal
get_disvoice_env <- function() {
  init_disvoice()
  .disvoice_env
}

#' Install DisVoice Python Dependencies
#'
#' Installs required Python packages (parselmouth, numpy) for DisVoice.
#' Uses reticulate's py_install() function.
#'
#' @param method Installation method passed to reticulate::py_install()
#' @param conda Path to conda executable (if using conda)
#' @param ... Additional arguments passed to reticulate::py_install()
#'
#' @return NULL (invisible), called for side effects
#'
#' @references
#' \insertCite{Jadoul2018}{superassp}
#'
#' \insertCite{OrozcoArroyave2018}{superassp}
#'
#' @export
#' @examples
#' \dontrun{
#' install_disvoice_python()
#' }
install_disvoice_python <- function(method = "auto", conda = "auto", ...) {
  requirements_file <- system.file(
    "python", "DisVoice", "requirements.txt",
    package = "superassp",
    mustWork = TRUE
  )

  # Read requirements
  requirements <- readLines(requirements_file)

  message("Installing DisVoice Python dependencies...")
  reticulate::py_install(
    packages = requirements,
    method = method,
    conda = conda,
    pip = TRUE,
    ...
  )

  message("DisVoice Python dependencies installed successfully.")
  invisible(NULL)
}

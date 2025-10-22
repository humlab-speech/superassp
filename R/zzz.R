##' @keywords internal
.onAttach <- function(libname, pkgname) {
  # Check voice_analysis module status on package load
  if (interactive()) {
    check_voice_analysis_status()
    check_covarep_status()
    check_voice_sauce_status()
  }
}

##' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Try to import covarep module (silent - availability checked when needed)
  covarep_module <<- NULL

  tryCatch({
    if (reticulate::py_module_available("covarep")) {
      covarep_module <<- reticulate::import("covarep", delay_load = TRUE)
    }
  }, error = function(e) {
    # Silent - covarep not required for other package functions
  })

  # Try to import voicesauce module (silent - availability checked when needed)
  voicesauce_module <<- NULL

  tryCatch({
    if (reticulate::py_module_available("voicesauce")) {
      voicesauce_module <<- reticulate::import("voicesauce", delay_load = TRUE)
    }
  }, error = function(e) {
    # Silent - voicesauce not required for other package functions
  })

  # Setup S7 method dispatch for DSP functions (lst_*, trk_*)
  # This enables AVAudio object support while maintaining backward compatibility
  tryCatch({
    .setup_s7_methods()
  }, error = function(e) {
    warning("Failed to setup S7 method dispatch: ", e$message, call. = FALSE)
  })

  # Register psychoacoustic units (Bark scale, etc.) with units package
  tryCatch({
    .onLoad_psychoacoustic_units()
  }, error = function(e) {
    # Silent - units package may not be available
  })
}

##' Check and report voice_analysis module status
##'
##' Internal function called on package attach to inform users about
##' voice_analysis module availability and optimization status.
##'
##' @keywords internal
check_voice_analysis_status <- function() {
  # Only check if module directory exists
  module_dir <- system.file("python", "voice_analysis_python", package = "superassp")

  if (!dir.exists(module_dir)) {
    return(invisible(NULL))  # Module not included in package
  }

  # Check if module is installed
  if (!voice_analysis_available()) {
    return(invisible(NULL))  # Don't spam on every load if not installed
  }

  # Module is installed - check optimization status
  tryCatch({
    info <- voice_analysis_info()

    # Build status message
    status_parts <- character(0)

    if (info$cython_available) {
      status_parts <- c(status_parts, "Cython")
    }

    if (info$numba_available) {
      status_parts <- c(status_parts, "Numba")
    }

    if (length(status_parts) > 0) {
      optimizations <- paste(status_parts, collapse = " + ")
      # Only show message occasionally (e.g., once per session)
      if (!exists(".superassp_vat_msg_shown", envir = .GlobalEnv)) {
        packageStartupMessage(
          sprintf("voice_analysis: %s optimizations active", optimizations)
        )
        assign(".superassp_vat_msg_shown", TRUE, envir = .GlobalEnv)
      }
    } else {
      # No optimizations - suggest installation with Cython
      if (!exists(".superassp_vat_warning_shown", envir = .GlobalEnv)) {
        packageStartupMessage(
          "voice_analysis: Running in pure Python mode.\n",
          "  For 2-3x speedup, install with: install_voice_analysis(method='cython')"
        )
        assign(".superassp_vat_warning_shown", TRUE, envir = .GlobalEnv)
      }
    }

  }, error = function(e) {
    # Silently fail - don't spam users with errors on load
    invisible(NULL)
  })

  invisible(NULL)
}

##' Check and report COVAREP module status
##'
##' Internal function called on package attach to inform users about
##' COVAREP module availability and optimization status.
##'
##' @keywords internal
check_covarep_status <- function() {
  # Only check if module directory exists
  module_dir <- system.file("python", "covarep_python", package = "superassp")

  if (!dir.exists(module_dir)) {
    return(invisible(NULL))  # Module not included in package
  }

  # Check if module is installed
  if (!covarep_available()) {
    return(invisible(NULL))  # Don't spam on every load if not installed
  }

  # Module is installed - check optimization status
  tryCatch({
    info <- covarep_info()

    # Build status message
    status_parts <- character(0)

    if (info$numba) {
      status_parts <- c(status_parts, "Numba")
    }

    if (length(status_parts) > 0) {
      optimizations <- paste(status_parts, collapse = " + ")
      # Only show message occasionally (e.g., once per session)
      if (!exists(".superassp_covarep_msg_shown", envir = .GlobalEnv)) {
        packageStartupMessage(
          sprintf("covarep: %s optimizations active", optimizations)
        )
        assign(".superassp_covarep_msg_shown", TRUE, envir = .GlobalEnv)
      }
    } else {
      # No optimizations - suggest installation with Numba
      if (!exists(".superassp_covarep_warning_shown", envir = .GlobalEnv)) {
        packageStartupMessage(
          "covarep: Running in pure Python mode.\n",
          "  For 5-10x speedup, install with: install_covarep(method='numba')"
        )
        assign(".superassp_covarep_warning_shown", TRUE, envir = .GlobalEnv)
      }
    }

  }, error = function(e) {
    # Silently fail - don't spam users with errors on load
    invisible(NULL)
  })

  invisible(NULL)
}

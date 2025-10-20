#' S7 Method System for DSP Functions
#'
#' This file sets up S7 method dispatch for all lst_* and trk_* functions,
#' allowing them to accept both character vectors (file paths) and AVAudio objects.
#'
#' The system works by:
#' 1. Converting each existing function to an S7 generic
#' 2. Registering the original implementation as the character method
#' 3. Adding an AVAudio method that converts to temp file and calls original
#'
#' This preserves full backward compatibility while adding AVAudio support.
#'
#' Note: listOfFiles is now a mandatory parameter in all DSP functions (no default value).
#'
#' @name s7-methods
NULL

#' Setup S7 Method Dispatch for DSP Functions
#'
#' Internal function called during .onLoad() to set up S7 dispatch for all
#' lst_* and trk_* functions. This enables them to work with AVAudio objects
#' while maintaining backward compatibility with file paths.
#'
#' @return NULL (called for side effects)
#' @keywords internal
.setup_s7_methods <- function() {

  # List of all DSP functions (will be populated dynamically)
  # Get all exported functions starting with lst_ or trk_
  ns <- getNamespace("superassp")
  all_fns <- ls(ns, all.names = FALSE)

  dsp_fns <- all_fns[grepl("^(lst_|trk_)", all_fns)]

  converted_count <- 0
  failed_count <- 0

  for (fn_name in dsp_fns) {
    result <- tryCatch({
      .convert_to_s7_generic(fn_name)
      converted_count <- converted_count + 1
      TRUE
    }, error = function(e) {
      # Only warn in interactive sessions or if verbose
      if (getOption("superassp.s7.verbose", FALSE)) {
        warning("Could not convert ", fn_name, " to S7 generic: ", e$message,
                call. = FALSE)
      }
      failed_count <- failed_count + 1
      FALSE
    })
  }

  # Report summary if verbose
  if (getOption("superassp.s7.verbose", FALSE)) {
    message("S7 dispatch setup: ", converted_count, " functions converted, ",
            failed_count, " skipped")
  }

  invisible(NULL)
}

#' Convert a DSP Function to S7 Generic
#'
#' Internal function to convert an existing DSP function to an S7 generic
#' with methods for character (file paths) and AVAudio objects.
#'
#' @param fn_name Character; name of the function
#' @return NULL (called for side effects)
#' @keywords internal
.convert_to_s7_generic <- function(fn_name) {

  ns <- getNamespace("superassp")

  # Get the original function
  if (!exists(fn_name, envir = ns, inherits = FALSE)) {
    return(invisible(NULL))
  }

  original_fn <- get(fn_name, envir = ns)

  # Skip if already an S7 generic (check for S7_generic class)
  if (inherits(original_fn, "S7_generic")) {
    return(invisible(NULL))
  }

  # Create S7 generic
  generic_fn <- S7::new_generic(
    name = fn_name,
    dispatch_args = "listOfFiles"
  )

  # Register character method (original implementation)
  S7::method(generic_fn, S7::class_character) <- original_fn

  # Register AVAudio method
  avaudio_method <- function(listOfFiles, ...) {
    # Convert AVAudio to temporary file
    temp_file <- avaudio_to_tempfile(listOfFiles, verbose = FALSE)

    # Ensure cleanup
    on.exit(unlink(temp_file), add = TRUE)

    # Call original function with temp file (as character vector)
    result <- original_fn(as.character(temp_file), ...)

    result
  }
  S7::method(generic_fn, AVAudio) <- avaudio_method

  # Replace function in namespace
  unlockBinding(fn_name, ns)
  assign(fn_name, generic_fn, envir = ns)
  lockBinding(fn_name, ns)

  invisible(NULL)
}

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
        cli::cli_warn("Could not convert {.fn {fn_name}} to S7 generic: {e$message}")
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


  # Save attributes from original function (ext, tracks, outputType, etc.)
  original_attrs <- attributes(original_fn)

  # Create S7 generic.
  #
  # The generic is built with the *original* function's formals rather than the
  # `(listOfFiles, ...)` shape S7 would synthesise. R CMD check compares the
  # installed signature against the documented \usage{} in each Rd file, so a
  # synthesised `(listOfFiles, ...)` reported a "code/documentation mismatch"
  # for every converted function. Dispatch still keys on `listOfFiles`.
  #
  # S7 requires dispatch_args to be a prefix of the generic's formals, so
  # functions whose first argument is not `listOfFiles` (lst_avqi, lst_dsi,
  # lst_vowel_space, ...) still fail conversion and are left untouched, exactly
  # as they were before this change.
  generic_formals <- formals(original_fn)

  generic_fun <- function() S7::S7_dispatch()
  formals(generic_fun) <- generic_formals

  generic_fn <- S7::new_generic(
    name = fn_name,
    dispatch_args = "listOfFiles",
    fun = generic_fun
  )

  # Register character method (original implementation)
  S7::method(generic_fn, S7::class_character) <- original_fn

  # Register AVAudio method.
  #
  # The helper methods below reuse the generic's formals rather than the
  # `(listOfFiles, ...)` shorthand: S7 rejects a method whose formals do not
  # match the generic exactly whenever the generic has no `...`, which is the
  # case for every lst_*/trk_* function that takes no free parameters.
  avaudio_method <- function(listOfFiles, ...) {
    # The `(listOfFiles, ...)` signature written here is only a placeholder:
    # `formals()` is replaced with the generic's below. Writing the dispatch
    # argument explicitly keeps it visible to R CMD check's static analysis,
    # which cannot see formals that are assigned at run time.
    forwarded <- as.list(match.call())[-1]
    forwarded$listOfFiles <- NULL

    # Convert AVAudio to temporary file
    temp_file <- avaudio_to_tempfile(listOfFiles, verbose = FALSE)

    # Ensure cleanup
    on.exit(unlink(temp_file), add = TRUE)

    # Call original function with temp file (as character vector)
    do.call(original_fn, c(list(as.character(temp_file)), forwarded))
  }
  formals(avaudio_method) <- generic_formals
  S7::method(generic_fn, AVAudio) <- avaudio_method

  # Fallback for unsupported input (NULL, numeric, logical, …): emit a clear
  # validation error instead of S7's cryptic "Can't find method" message.
  # Character and AVAudio are more specific, so they still take precedence.
  unsupported_input_method <- function(listOfFiles, ...) {
    cli::cli_abort(
      "No input files specified: {.arg listOfFiles} must be a character vector of file paths or an AVAudio object."
    )
  }
  formals(unsupported_input_method) <- generic_formals
  S7::method(generic_fn, S7::class_any) <- unsupported_input_method



  # Restore custom attributes to S7 generic
  for (attr_name in c("ext", "tracks", "outputType", "nativeFiletypes", "suggestCaching")) {
    if (!is.null(original_attrs[[attr_name]])) {
      attr(generic_fn, attr_name) <- original_attrs[[attr_name]]
    }
  }


  # Replace function in namespace
  unlockBinding(fn_name, ns)
  assign(fn_name, generic_fn, envir = ns)
  lockBinding(fn_name, ns)

  invisible(NULL)
}

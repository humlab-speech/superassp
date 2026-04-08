#' Error Message Formatting Helpers
#'
#' @description
#' This file provides helper functions for consistent error and warning
#' message formatting throughout superassp.
#'
#' @name error_helpers
#' @keywords internal
NULL

#' Format processing error message
#'
#' Creates a consistently formatted error message for file processing errors.
#' Includes the file name, error message, and optional context about what
#' operation was being performed.
#'
#' @param file_path Character string with the file path that caused the error
#' @param error_msg Character string with the error message
#' @param context Character string with optional context (e.g., "COVAREP computation", "pitch tracking")
#' @param include_path Logical indicating whether to include full path (default: FALSE, shows basename only)
#'
#' @return Character string with formatted error message
#' @keywords internal
#' @noRd
#'
#' @examples
#' \dontrun{
#' tryCatch({
#'   result <- compute_something(file)
#' }, error = function(e) {
#'   warning(format_processing_error(file, e$message, "COVAREP computation"))
#' })
#' }
format_processing_error <- function(file_path, error_msg, context = NULL,
                                    include_path = FALSE) {

  # Use basename by default for cleaner messages
  file_display <- if (include_path) file_path else basename(file_path)

  # Build base message
  base_msg <- sprintf("Error processing '%s': %s", file_display, error_msg)

  # Add context if provided
  if (!is.null(context) && nchar(context) > 0) {
    base_msg <- sprintf("%s\n  Context: %s", base_msg, context)
  }

  base_msg
}

#' Format processing warning message
#'
#' Creates a consistently formatted warning message for file processing issues.
#' Similar to format_processing_error but for non-fatal issues.
#'
#' @param file_path Character string with the file path
#' @param warning_msg Character string with the warning message
#' @param context Character string with optional context
#' @param include_path Logical indicating whether to include full path (default: FALSE)
#'
#' @return Character string with formatted warning message
#' @keywords internal
#' @noRd
format_processing_warning <- function(file_path, warning_msg, context = NULL,
                                     include_path = FALSE) {

  file_display <- if (include_path) file_path else basename(file_path)

  base_msg <- sprintf("Warning for '%s': %s", file_display, warning_msg)

  if (!is.null(context) && nchar(context) > 0) {
    base_msg <- sprintf("%s\n  Context: %s", base_msg, context)
  }

  base_msg
}

#' Format batch processing summary
#'
#' Creates a summary message for batch processing operations showing
#' success/failure counts.
#'
#' @param n_total Integer total number of files processed
#' @param n_success Integer number of successful files
#' @param n_failed Integer number of failed files
#' @param operation Character string describing the operation (e.g., "pitch tracking", "JSTF writing")
#'
#' @return Character string with formatted summary
#' @keywords internal
#' @noRd
#'
#' @examples
#' \dontrun{
#' message(format_batch_summary(10, 8, 2, "pitch tracking"))
#' # "Batch pitch tracking complete: 8/10 succeeded, 2 failed"
#' }
format_batch_summary <- function(n_total, n_success, n_failed, operation = "processing") {

  if (n_failed == 0) {
    sprintf("Batch %s complete: %d/%d succeeded", operation, n_success, n_total)
  } else {
    sprintf("Batch %s complete: %d/%d succeeded, %d failed",
            operation, n_success, n_total, n_failed)
  }
}

#' Format parameter validation error
#'
#' Creates a consistently formatted error message for parameter validation failures.
#'
#' @param param_name Character string with the parameter name
#' @param expected Character string describing expected value/type
#' @param actual Character string describing actual value/type
#' @param function_name Character string with the calling function name
#'
#' @return Character string with formatted error message
#' @keywords internal
#' @noRd
format_validation_error <- function(param_name, expected, actual,
                                   function_name = "function") {
  sprintf("%s(): Parameter '%s' validation failed\n  Expected: %s\n  Got: %s",
          function_name, param_name, expected, actual)
}

#' Safe error extraction
#'
#' Safely extracts error message from an error object, handling various error types.
#'
#' @param error Error object from tryCatch
#'
#' @return Character string with error message
#' @keywords internal
#' @noRd
safe_error_message <- function(error) {
  if (inherits(error, "error")) {
    return(conditionMessage(error))
  } else if (inherits(error, "simpleError")) {
    return(error$message)
  } else if (is.character(error)) {
    return(error)
  } else {
    return(as.character(error))
  }
}


#' Format file I/O error
#'
#' Creates a formatted error message for file I/O operations (read/write).
#'
#' @param file_path Character string with the file path
#' @param operation Character string describing the operation ("read", "write", "create", "delete")
#' @param error_msg Character string with the error message
#'
#' @return Character string with formatted error message
#' @keywords internal
#' @noRd
format_io_error <- function(file_path, operation = "access", error_msg = NULL) {

  base_msg <- sprintf("Failed to %s file: %s", operation, file_path)

  if (!is.null(error_msg) && nchar(error_msg) > 0) {
    base_msg <- sprintf("%s\n  Error: %s", base_msg, error_msg)
  }

  # Add helpful hints based on operation
  if (operation == "read") {
    base_msg <- sprintf("%s\n  Check: File exists and is readable", base_msg)
  } else if (operation == "write" || operation == "create") {
    base_msg <- sprintf("%s\n  Check: Directory exists and is writable", base_msg)
  }

  base_msg
}

#' Create error handler for batch processing
#'
#' Creates a standardized error handler function for use in batch processing loops.
#' Returns a function that can be used in tryCatch error handlers.
#'
#' @param operation Character string describing the operation
#' @param verbose Logical indicating whether to print errors
#'
#' @return Function that handles errors consistently
#' @keywords internal
#' @noRd
#'
#' @examples
#' \dontrun{
#' error_handler <- create_error_handler("pitch tracking", verbose = TRUE)
#'
#' result <- tryCatch({
#'   compute_pitch(file)
#' }, error = error_handler(file))
#' }
create_error_handler <- function(operation, verbose = TRUE) {
  function(file_path) {
    function(error) {
      msg <- format_processing_error(
        file_path,
        safe_error_message(error),
        operation
      )

      if (verbose) {
        warning(msg, call. = FALSE)
      }

      # Return NULL or error indicator
      return(NULL)
    }
  }
}


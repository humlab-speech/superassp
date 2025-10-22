# Track Attribute Helpers
#
# Helper functions for setting track attributes on returned AsspDataObj objects.

#' Set tracks attribute on AsspDataObj result
#'
#' Internal helper to propagate the tracks attribute from a function to its
#' returned AsspDataObj object(s). This enables proper template expansion
#' in `as.data.frame.AsspDataObj()`.
#'
#' @param result AsspDataObj or list of AsspDataObj objects
#' @param func Function object whose tracks attribute should be copied
#' @param n_files Number of files processed (for list handling)
#'
#' @return The result with tracks attribute set
#'
#' @keywords internal
.set_tracks_attribute <- function(result, func, n_files = 1) {
  tracks_attr <- attr(func, "tracks")

  if (is.null(tracks_attr)) {
    # No tracks attribute on function, nothing to do
    return(result)
  }

  # Handle function-based tracks (MFCC, snack formants)
  if (is.function(tracks_attr)) {
    # These use dynamic naming, can't set static attribute
    # The function will be called at as.data.frame() time
    return(result)
  }

  # Set attribute on result object(s)
  if (n_files == 1) {
    attr(result, "tracks") <- tracks_attr
  } else {
    # List of objects
    for (i in seq_along(result)) {
      attr(result[[i]], "tracks") <- tracks_attr
    }
  }

  result
}

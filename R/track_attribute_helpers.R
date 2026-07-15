# Track Attribute Helpers
#
# Helper functions for setting track attributes on returned AsspDataObj objects.

#' Construct an AsspDataObj from track data + metadata
#'
#' Single source of truth for the hand-built `list(...) + attr() + class()`
#' pattern repeated across the `create_*_asspobj` helpers. Sets exactly the
#' attributes it is given (optional ones are omitted when `NULL`), so it
#' reproduces each caller's previous attribute set byte-for-byte — attribute
#' order does not affect `write.AsspDataObj()` output, which keys by name.
#'
#' @param tracks Named list of track matrices/vectors (the SSFF tracks).
#' @param sampleRate Frame rate of the tracks (frames per second).
#' @param trackFormats Character vector of per-track storage formats
#'   (e.g. \code{"REAL32"}), one per track.
#' @param startTime Start time of the first frame in seconds (default 0).
#' @param startRecord First record index (default 1).
#' @param endRecord Last record index. Defaults to the row count of the first
#'   track.
#' @param origFreq Original audio sample rate in Hz. Omitted if \code{NULL}.
#' @param fileInfo Integer vector \code{c(fileFormat, nTracks)}. Omitted if
#'   \code{NULL}.
#' @param filePath Source file path. Omitted if \code{NULL}.
#'
#' @return An object of class \code{AsspDataObj}.
#' @keywords internal
#' @noRd
new_asspdataobj <- function(tracks, sampleRate, trackFormats,
                            startTime = 0.0, startRecord = 1L,
                            endRecord = NULL, origFreq = NULL,
                            fileInfo = NULL, filePath = NULL) {
  obj <- tracks
  if (is.null(endRecord)) endRecord <- NROW(tracks[[1]])

  attr(obj, "trackFormats") <- trackFormats
  attr(obj, "sampleRate")   <- sampleRate
  if (!is.null(origFreq)) attr(obj, "origFreq") <- origFreq
  attr(obj, "startTime")    <- startTime
  attr(obj, "startRecord")  <- as.integer(startRecord)
  attr(obj, "endRecord")    <- as.integer(endRecord)
  if (!is.null(fileInfo)) attr(obj, "fileInfo") <- fileInfo
  if (!is.null(filePath)) attr(obj, "filePath") <- filePath
  class(obj) <- "AsspDataObj"
  obj
}

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

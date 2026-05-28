#' Accessor methods for AsspDataObj and JsonTrackObj
#'
#' @description
#' Read-only accessors that work on both \link{AsspDataObj} (SSFF signal data)
#' and \link{JsonTrackObj} (JSTF summary objects).
#'
#' | Function | Returns |
#' | :--- | :--- |
#' | `sample_rate(x)` | Sample rate in Hz |
#' | `n_records(x)` | Number of records / slices |
#' | `signal_duration(x)` | Duration in seconds |
#' | `start_time(x)` | Start time of first record/slice |
#' | `track_names(x)` | Track (field) names |
#' | `file_path(x)` | Source audio file path |
#' | `track_formats(x)` | Data-type string per track |
#'
#' @section Deprecated aliases (since v2.8.0):
#' | Old | New |
#' | :--- | :--- |
#' | `rate(x)` | `sample_rate(x)` |
#' | `numRecs(x)` | `n_records(x)` |
#' | `dur(x)` | `signal_duration(x)` |
#' | `startTime(x)` | `start_time(x)` |
#' | `tracks(x)` | `track_names(x)` |
#'
#' @param x An \link{AsspDataObj} or \link{JsonTrackObj}.
#' @param ... Currently unused; reserved for future extensions.
#' @name assp_accessors
#' @aliases assp_accessors
NULL

#' @rdname assp_accessors
#' @export
sample_rate <- function(x, ...) UseMethod("sample_rate")

#' @rdname assp_accessors
#' @exportS3Method
sample_rate.AsspDataObj <- function(x, ...) attr(x, "sampleRate")

#' @rdname assp_accessors
#' @exportS3Method
sample_rate.JsonTrackObj <- function(x, ...) x$sample_rate


#' @rdname assp_accessors
#' @export
n_records <- function(x, ...) UseMethod("n_records")

#' @rdname assp_accessors
#' @exportS3Method
n_records.AsspDataObj <- function(x, ...) {
  attr(x, "endRecord") - attr(x, "startRecord") + 1L
}

#' @rdname assp_accessors
#' @exportS3Method
n_records.JsonTrackObj <- function(x, ...) length(x$slices)


#' @rdname assp_accessors
#' @export
signal_duration <- function(x, ...) UseMethod("signal_duration")

#' @rdname assp_accessors
#' @exportS3Method
signal_duration.AsspDataObj <- function(x, ...) {
  n_records.AsspDataObj(x) / attr(x, "sampleRate")
}

#' @rdname assp_accessors
#' @exportS3Method
signal_duration.JsonTrackObj <- function(x, ...) x$audio_duration


#' @rdname assp_accessors
#' @export
start_time <- function(x, ...) UseMethod("start_time")

#' @rdname assp_accessors
#' @exportS3Method
start_time.AsspDataObj <- function(x, ...) attr(x, "startTime")

#' @rdname assp_accessors
#' @exportS3Method
start_time.JsonTrackObj <- function(x, ...) {
  if (length(x$slices) > 0L) x$slices[[1L]]$begin_time else 0.0
}


#' @rdname assp_accessors
#' @export
track_names <- function(x, ...) UseMethod("track_names")

#' @rdname assp_accessors
#' @exportS3Method
track_names.AsspDataObj <- function(x, ...) names(x)

#' @rdname assp_accessors
#' @exportS3Method
track_names.JsonTrackObj <- function(x, ...) names(x$field_schema)


#' @rdname assp_accessors
#' @export
file_path <- function(x, ...) UseMethod("file_path")

#' @rdname assp_accessors
#' @exportS3Method
file_path.AsspDataObj <- function(x, ...) attr(x, "filePath")

#' @rdname assp_accessors
#' @exportS3Method
file_path.JsonTrackObj <- function(x, ...) x$file_path


#' @rdname assp_accessors
#' @export
track_formats <- function(x, ...) UseMethod("track_formats")

#' @rdname assp_accessors
#' @exportS3Method
track_formats.AsspDataObj <- function(x, ...) attr(x, "trackFormats")

#' @rdname assp_accessors
#' @exportS3Method
track_formats.JsonTrackObj <- function(x, ...) unlist(x$field_schema)

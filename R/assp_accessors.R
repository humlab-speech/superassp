#' @describeIn AsspDataObj Sample rate in Hz.
#' @param x AsspDataObj
#' @param ... Ignored.
#' @export
sample_rate <- function(x, ...) UseMethod("sample_rate")

#' @exportS3Method
sample_rate.AsspDataObj <- function(x, ...) attr(x, "sampleRate")

#' @exportS3Method
sample_rate.JsonTrackObj <- function(x, ...) x$sample_rate


#' @describeIn AsspDataObj Number of records (frames).
#' @param x AsspDataObj
#' @param ... Ignored.
#' @export
n_records <- function(x, ...) UseMethod("n_records")

#' @exportS3Method
n_records.AsspDataObj <- function(x, ...) {
  attr(x, "endRecord") - attr(x, "startRecord") + 1L
}


#' @describeIn AsspDataObj Duration in seconds.
#' @param x AsspDataObj
#' @param ... Ignored.
#' @export
signal_duration <- function(x, ...) UseMethod("signal_duration")

#' @exportS3Method
signal_duration.AsspDataObj <- function(x, ...) {
  n_records.AsspDataObj(x) / attr(x, "sampleRate")
}

#' @exportS3Method
signal_duration.JsonTrackObj <- function(x, ...) x$audio_duration


#' @describeIn AsspDataObj Start time of the first record in seconds.
#' @param x AsspDataObj
#' @param ... Ignored.
#' @export
start_time <- function(x, ...) UseMethod("start_time")

#' @exportS3Method
start_time.AsspDataObj <- function(x, ...) attr(x, "startTime")


#' @describeIn AsspDataObj Track names (equivalent to names()).
#' @param x AsspDataObj
#' @param ... Ignored.
#' @export
track_names <- function(x, ...) UseMethod("track_names")

#' @exportS3Method
track_names.AsspDataObj <- function(x, ...) names(x)

#' @exportS3Method
track_names.JsonTrackObj <- function(x, ...) names(x$field_schema)


#' @describeIn AsspDataObj Source file path, or NULL for in-memory objects.
#' @param x AsspDataObj
#' @param ... Ignored.
#' @export
file_path <- function(x, ...) UseMethod("file_path")

#' @exportS3Method
file_path.AsspDataObj <- function(x, ...) attr(x, "filePath")

#' @exportS3Method
file_path.JsonTrackObj <- function(x, ...) x$file_path


#' @describeIn AsspDataObj SSFF data-type string per track.
#' @param x AsspDataObj
#' @param ... Ignored.
#' @export
track_formats <- function(x, ...) UseMethod("track_formats")

#' @exportS3Method
track_formats.AsspDataObj <- function(x, ...) attr(x, "trackFormats")

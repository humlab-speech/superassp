#' @describeIn AsspDataObj Return the duration in seconds.
#' @param x AsspDataObj
#' @param ... additional arguments (ignored)
#' @export
dur <- function(x, ...) {
  lifecycle::deprecate_warn("2.8.0", "dur()", "signal_duration()")
  UseMethod("dur")
}

#' @describeIn AsspDataObj Return the number of records.
#' @param x AsspDataObj
#' @param ... Additional arguments (ignored; present for S3 generic compatibility).
#' @export
numRecs <- function(x, ...) {
  lifecycle::deprecate_warn("2.8.0", "numRecs()", "n_records()")
  UseMethod("numRecs")
}

#' @describeIn AsspDataObj Return the sample rate in Hz.
#' @param x AsspDataObj
#' @param ... Additional arguments (ignored; present for S3 generic compatibility).
#' @export
rate <- function(x, ...) {
  lifecycle::deprecate_warn("2.8.0", "rate()", "sample_rate()")
  UseMethod("rate")
}

#' @describeIn AsspDataObj Return the start time of the first sample.
#' @param x AsspDataObj
#' @param ... Additional arguments (ignored; present for S3 generic compatibility).
#' @export
startTime <- function(x, ...) {
  lifecycle::deprecate_warn("2.8.0", "startTime()", "start_time()")
  UseMethod("startTime")
}

#' @describeIn AsspDataObj Return track names (equivalent to names()).
#' @param x AsspDataObj
#' @param ... Additional arguments (ignored; present for S3 generic compatibility).
#' @export
tracks <- function(x, ...) {
  lifecycle::deprecate_warn("2.8.0", "tracks()", "track_names()")
  UseMethod("tracks")
}

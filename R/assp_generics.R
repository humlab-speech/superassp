#' @rdname assp_accessors
#' @export
dur <- function(x, ...) {
  lifecycle::deprecate_warn("2.8.0", "dur()", "signal_duration()")
  UseMethod("dur")
}

#' @rdname assp_accessors
#' @export
numRecs <- function(x, ...) {
  lifecycle::deprecate_warn("2.8.0", "numRecs()", "n_records()")
  UseMethod("numRecs")
}

#' @rdname assp_accessors
#' @export
rate <- function(x, ...) {
  lifecycle::deprecate_warn("2.8.0", "rate()", "sample_rate()")
  UseMethod("rate")
}

#' @rdname assp_accessors
#' @export
startTime <- function(x, ...) {
  lifecycle::deprecate_warn("2.8.0", "startTime()", "start_time()")
  UseMethod("startTime")
}

#' @rdname assp_accessors
#' @export
tracks <- function(x, ...) {
  lifecycle::deprecate_warn("2.8.0", "tracks()", "track_names()")
  UseMethod("tracks")
}

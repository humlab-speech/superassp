#' @describeIn AsspDataObj Return the duration in seconds.
#' @param x AsspDataObj
#' @param ... additional arguments (ignored)
#' @export
dur <- function(x, ...) UseMethod("dur")

#' @describeIn AsspDataObj Return the number of records.
#' @export
numRecs <- function(x, ...) UseMethod("numRecs")

#' @describeIn AsspDataObj Return the sample rate in Hz.
#' @export
rate <- function(x, ...) UseMethod("rate")

#' @describeIn AsspDataObj Return the start time of the first sample.
#' @export
startTime <- function(x, ...) UseMethod("startTime")

#' @describeIn AsspDataObj Return track names (equivalent to names()).
#' @export
tracks <- function(x, ...) UseMethod("tracks")

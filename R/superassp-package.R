#' superassp: Praat, MATLAB(TM) and wrassp speech signal processing using a wrassp-like interface
#'
#' This package bundles together routines that utilizes Praat, matlab or python
#' to do the signal processing. The output of the original routines are then
#' wrapped into a common interface that outputs an SSFF file that is compatible
#' with the output of the wrassp package functions. This package also re-exports
#' all of wrassp's functions so it is usually enough to just load superassp.
#' 
#'
#' @docType package
#' @name superassp-package
#' @importFrom Rdpack reprompt
#' 
NULL
#> NULL

PRAAT_DEVEL = FALSE
logger::log_threshold(logger::WARN)
#' superassp: speech signal processing using various framworks using a wrassp-like interface
#'
#' This package bundles together routines that utilizes Praat, matlab or python
#' to do the signal processing. The output of the original routines are then
#' wrapped into a common interface that outputs an SSFF file or returns a wrassp compatible AsspDataObj object.
#' 
#'
#' @docType package
#' @importFrom Rdpack reprompt
#' @name superassp-package
NULL
#> NULL

PRAAT_DEVEL = FALSE
logger::log_threshold(logger::WARN)



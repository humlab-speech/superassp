##' returns all valid AsspWindowTypes according to the assp library
##'
##' wrapper function for AsspWindowTypes of wrassp
##' @title AsspWindowTypes
##' @return vector containing window types
##' @author Raphael Winkelmann
##' @useDynLib superassp, .registration = TRUE
##' @export
'AsspWindowTypes' <- function(){
	
	return(.Call("AsspWindowTypes_", PACKAGE = "superassp"))
  
}

##' returns all valid AsspLpTypes according to the assp library
##'
##' wrapper function for AsspLpTypes of wrassp
##' @title AsspLpTypes
##' @return vector containing lp types
##' @author Raphael Winkelmann
##' @useDynLib superassp, .registration = TRUE
##' @export
'AsspLpTypes' <- function(){
	
	return(.Call("AsspLpTypes_", PACKAGE = "superassp"))

}

##' returns all valid AsspSpectTypes according to the assp library
##'
##' wrapper function for AsspSpectTypes of wrassp
##' @title AsspSpectTypes
##' @return vector containing spectrogram types
##' @author Raphael Winkelmann
##' @useDynLib superassp, .registration = TRUE
##' @export
'AsspSpectTypes' <- function(){
	
	return(.Call("AsspSpectTypes_", PACKAGE = "superassp"))

}




#' Procedure reporting lossless file formats
#'
#' This procedure reports file extensions that are known to contain
#' losslessly encoded sound data. These formats preserve the original
#' audio signal without quality loss, which is essential for accurate
#' DSP (Digital Signal Processing) analysis.
#'
#' The list includes formats supported by the av package and native
#' ASSP/wrassp formats:
#'
#' @details
#' **Container formats (lossless codecs):**
#' - wav: WAV / WAVE (Waveform Audio) - most common
#' - flac: Free Lossless Audio Codec
#' - aiff: Audio Interchange File Format (Apple)
#' - wv: WavPack
#' - ape: Monkey's Audio
#' - tta: True Audio
#' - caf: Apple Core Audio Format
#' - au: Sun/NeXT Audio
#' - w64: Sony Wave64 (64-bit WAV variant)
#'
#' **High-resolution audio:**
#' - dsf: DSD Stream File (Direct Stream Digital)
#' - dff: DSD Interchange File Format
#'
#' **Professional/scientific formats:**
#' - kay: Kay Elemetrics CSL files
#' - nist: NIST SPHERE
#' - nsp: NSP (used in speech research)
#'
#' @return Character vector of file extensions (without dots)
#' @export
#' @author Fredrik Nylén
#' @examples
#' # Get list of lossless formats
#' knownLossless()
#'
#' # Check if a file extension is lossless
#' "flac" %in% knownLossless()  # TRUE
#' "mp3" %in% knownLossless()   # FALSE

knownLossless <- function() {
  c(
    # Common lossless formats (av + wrassp supported)
    "wav", "flac", "aiff", "wv", "ape", "tta", "caf", "au", "w64",
    # High-resolution audio
    "dsf", "dff",
    # Professional/scientific formats
    "kay", "nist", "nsp"
  )
}
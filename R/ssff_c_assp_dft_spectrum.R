##' Track short-term DFT power spectrum
##'
##' Computes a short-term power spectrum via the Fast Fourier Transform using
##' the *libassp* C library \insertCite{s5h}{superassp}. Produces an
##' unsmoothed narrow-band spectrum from 0 Hz to the Nyquist rate. Prefer this
##' function when raw spectral detail is needed; use \code{trk_css_spectrum} or
##' \code{trk_lps_spectrum} for smoothed spectral envelopes.
##'
##' @inheritParams trk_lps_spectrum
##' @param bandwidth Numeric. Effective analysis bandwidth in Hz. Default 0
##'   yields the minimum bandwidth determined by the FFT length.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
##'   \describe{
##'     \item{\code{DFT[dB]}}{REAL32, dB power, n_frames x (FFT_length/2 + 1) columns.
##'       Power spectral amplitude from 0 Hz to the Nyquist rate.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 200 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @details
##' The FFT length is determined by \code{resolution} unless overridden by
##' \code{fftLength}. \code{bandwidth} widens the effective analysis window,
##' trading spectral resolution for reduced side-lobe leakage.
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##'
##' @seealso [wrassp::dftSpectrum]
##' @seealso [superassp::AsspWindowTypes]
##' @seealso [av::av_audio_convert]
##'
##' @useDynLib superassp, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("extdata", package = "wrassp"),
##'                        pattern = glob2rx("*.wav"),
##'                        full.names = TRUE)[1]
##'
##' # calculate dft spectrum
##' res <- trk_dft_spectrum(path2wav, toFile=FALSE)
##'
##' # plot spectral values at midpoint of signal
##' plot(res$dft[dim(res$dft)[1]/2,],
##'      type='l',
##'      xlab='spectral value index',
##'      ylab='spectral value')
##'
##' @export
##'
##' @references
##'   \insertAllCited{}
##'
'trk_dft_spectrum' <- function(listOfFiles = NULL,
                          beginTime = 0.0,
                          centerTime = FALSE,
                          endTime = 0.0,
                          resolution = 40.0,
                          fftLength = 0,
                          windowShift = 5.0,
                          window = 'BLACKMAN',
                          bandwidth = 0.0,
                          toFile = TRUE,
                          explicitExt = "dft",
                          outputDirectory = NULL,
                          assertLossless = NULL,
                          logToFile = FALSE,
                          keepConverted=FALSE,
                          convertOverwrites=FALSE,
                          verbose = TRUE){

  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"dft",explicitExt)
  newTracknames <- "DFT[dB]"
  nativeFiletypes <- c("wav","au","kay","nist","nsp")

  if(!isAsspWindowType(toupper(window))){
    cli::cli_abort(c("WindowFunction of type {.val {window}} is not supported!",
                     "i"="Accepted window types for routines implemented in *libassp* are {.field {AsspWindowTypes()}}.")
    )
  }

  ## Initial constants -- generics
  currCall <- rlang::current_call()
  funName <- rlang::call_name(currCall)
  preferedFiletype <- nativeFiletypes[[1]]

  knownLossless <- c(assertLossless,knownLossless())

  # Normalize time parameters
  beginTime <- if(is.null(beginTime)) 0.0 else beginTime
  endTime <- if(is.null(endTime)) 0.0 else endTime

  n_files <- length(listOfFiles)

  # Validate time parameter lengths
  if(length(beginTime) > 1 && length(beginTime) != n_files) {
    cli::cli_abort("The {.par beginTime} must be length 1 or match {.par listOfFiles} length.")
  }
  if(length(endTime) > 1 && length(endTime) != n_files) {
    cli::cli_abort("The {.par endTime} must be length 1 or match {.par listOfFiles} length.")
  }

  # Use Rcpp for efficient time parameter recycling
  beginTime <- fast_recycle_times(beginTime, n_files)
  endTime <- fast_recycle_times(endTime, n_files)

  #### Setup logging ####
  makeOutputDirectory(outputDirectory, logToFile, funName)

  #### Use unified memory-based processing for all files ####
  if(verbose) {
    format_apply_msg(funName, n_files, beginTime, endTime)
  }

  # Use unified load-and-process helper (works for all file formats)
  result <- processMediaFiles_LoadAndProcess(
    listOfFiles = listOfFiles,
    beginTime = beginTime,
    endTime = endTime,
    nativeFiletypes = nativeFiletypes,
    fname = "spectrum",
    toFile = toFile,
    verbose = verbose,
    centerTime = centerTime,
    resolution = resolution,
    fftLength = as.integer(fftLength),
    windowShift = windowShift,
    window = window,
    bandwidth = bandwidth,
    explicitExt = explicitExt,
    outputDirectory = outputDirectory
  )

  externalRes <- result$externalRes
  listOfFilesDF <- result$listOfFilesDF
  toClear <- character(0)  # No files to clean up with load-and-process

  # Use Rcpp for fast track renaming (only when data is returned, not written to file)
  if(!toFile && !is.null(newTracknames)) {
    n_tracks <- length(names(externalRes[[1]]))
    if(n_tracks != length(newTracknames)) {
      cli::cli_abort(c(
        "Wrong number of track names supplied:",
        "i" = "Track{?s} named: {.field {names(externalRes[[1]])}}"
      ))
    }
    externalRes <- fast_rename_tracks(externalRes, newTracknames)
  }

  # Note: When toFile=TRUE, the C code writes files directly and returns 0
  # No need to write files again here

  # Simplify output for single file
  if(n_files == 1) externalRes <- externalRes[[1]]

  #### Cleanup ####
  cleanupConvertedInputMediaFiles(toClear, keepConverted, verbose)

  return(externalRes)
}
attr(trk_dft_spectrum,"ext") <-  "dft"
attr(trk_dft_spectrum,"tracks") <-  c("DFT[dB]")
attr(trk_dft_spectrum,"outputType") <-  "SSFF"
attr(trk_dft_spectrum,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(trk_dft_spectrum,"suggestCaching") <-  FALSE

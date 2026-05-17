##' Track cepstrally-smoothed spectrum
##'
##' Computes a short-term power spectrum smoothed by liftering in the cepstral
##' domain, using the *libassp* C library \insertCite{s5h}{superassp}. Cepstral
##' smoothing suppresses spectral fine structure (harmonics), leaving the vocal
##' tract envelope; prefer \code{trk_css_spectrum} over \code{trk_dft_spectrum}
##' when a smooth spectral envelope is needed without LP assumptions.
##'
##' @inheritParams trk_cepstrum
##' @param numCeps Integer. Number of cepstral coefficients used for liftering.
##'   Default 0 sets this to sample rate in kHz + 1 (minimum 2).
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
##'   \describe{
##'     \item{\code{CSS[dB]}}{REAL32, dB, n_frames x (FFT_length/2 + 1) columns.
##'       Cepstrally-smoothed power spectral amplitude from 0 Hz to the Nyquist
##'       rate.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 200 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @details
##' \code{numCeps} controls the number of cepstral terms retained before
##' back-transforming; lower values give a smoother envelope. The FFT size is
##' governed by \code{resolution} (or overridden by \code{fftLength}).
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##'
##' @seealso [wrassp::cssSpectrum]
##' @seealso [superassp::AsspWindowTypes]
##' @seealso [av::av_audio_convert]
##'
##' @useDynLib superassp, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
##'
##' # calculate cepstrally smoothed spectrum
##' res <- trk_css_spectrum(path2wav, toFile=FALSE)
##' resolution <- attr(res,"origFreq") / ncol(res[[1]])
##'
##' # plot spectral values at midpoint of signal
##' plot(y=res[["CSS[dB]"]][400,],
##'     x=seq(1,ncol(res[[1]]),1)* resolution,
##'     type='l',
##'     xlab='Frequency (Hz)',
##'     ylab='Amplitude (dB)')
##'
##' @export
##'

trk_css_spectrum <- function(listOfFiles,
                          beginTime = 0.0,
                          centerTime = FALSE,
                          endTime = 0.0,
                          resolution = 40.0,
                          fftLength = 0,
                          windowShift = 5.0,
                          numCeps = 0,
                          window = 'BLACKMAN',
                          toFile = TRUE,
                          explicitExt = "css",
                          outputDirectory = NULL,
                          assertLossless = NULL,
                          logToFile = FALSE,
                          keepConverted=FALSE,
                          convertOverwrites=FALSE,
                          verbose = TRUE){

  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"css",explicitExt)
  newTracknames <- "CSS[dB]"
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
    spectrumType = 'CSS',
    resolution = resolution,
    fftLength = as.integer(fftLength),
    windowShift = windowShift,
    window = window,
    numCeps = as.integer(numCeps),
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
attr(trk_css_spectrum,"ext") <-  "css"
attr(trk_css_spectrum,"tracks") <-  c("CSS[dB]")
attr(trk_css_spectrum,"outputType") <-  "SSFF"
attr(trk_css_spectrum,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(trk_css_spectrum,"suggestCaching") <-  FALSE

##' Track LP-smoothed spectrum
##'
##' Computes a short-term spectral envelope smoothed by linear prediction using
##' the *libassp* C library \insertCite{s5h}{superassp}. LP smoothing produces a
##' parametric model-based envelope; prefer \code{trk_lps_spectrum} over
##' \code{trk_css_spectrum} when the LP order is known and explicit control over
##' the number of resonances is desired.
##'
##' @param listOfFiles Character vector of audio file paths. Any format supported by
##'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
##' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
##' @param centerTime Numeric or logical. Single-frame analysis time point in seconds;
##'   overrides \code{beginTime}, \code{endTime}, and \code{windowShift}. Default \code{FALSE}.
##' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
##' @param resolution Numeric. Target FFT frequency resolution in Hz. Default 40.0.
##' @param fftLength Integer. Explicit FFT length in points; overrides \code{resolution}.
##'   Default 0 (use \code{resolution}).
##' @param windowSize Numeric. Analysis window size in milliseconds. Default 20 ms.
##' @param windowShift Numeric. Frame shift in milliseconds; sets output frame rate
##'   (1000 / windowShift Hz). Default 5 ms.
##' @param window Character. Analysis window function type. Default \code{"BLACKMAN"}.
##'   See [superassp::AsspWindowTypes].
##' @param order Integer. LP prediction order; 0 defaults to sample rate in kHz + 3.
##'   Default 0.
##' @param preemphasis Numeric. Pre-emphasis factor applied before LP analysis.
##'   Default -0.95.
##' @param deemphasize Logical. Undo the spectral tilt introduced by pre-emphasis
##'   when computing the output spectrum. Default \code{TRUE}.
##' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
##'   count written (invisibly). If \code{FALSE}, return an \code{AsspDataObj}.
##'   Default \code{TRUE}.
##' @param explicitExt Character. Output file extension. Default \code{"lps"}.
##' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
##'   writes alongside the input file.
##' @param assertLossless Character vector of additional file extensions to treat as
##'   losslessly encoded.
##' @param logToFile Logical. Write processing log to a file in \code{outputDirectory}
##'   rather than the console. Default \code{FALSE}.
##' @param keepConverted Logical. Retain intermediate transcoded files. Default \code{FALSE}.
##' @param convertOverwrites Logical. Allow transcoding to overwrite existing files.
##'   Default \code{FALSE}.
##' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
##'   \describe{
##'     \item{\code{LPS[dB]}}{REAL32, dB, n_frames x (FFT_length/2 + 1) columns.
##'       LP-smoothed spectral amplitude from 0 Hz to the Nyquist rate.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 200 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @details
##' \code{order = 0} sets the LP order to sample rate in kHz + 3. When
##' \code{deemphasize = TRUE} (default), the output spectrum is corrected for
##' the pre-emphasis spectral tilt applied before LP analysis.
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##'
##' @seealso [wrassp::lpsSpectrum]
##' @seealso [superassp::AsspWindowTypes]
##' @seealso [av::av_audio_convert]
##'
##' @useDynLib superassp, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
##'
##' # calculate linear prediction smoothed spectrum
##' res <- trk_lps_spectrum(path2wav, toFile=FALSE)
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

trk_lps_spectrum <- function(listOfFiles,
                          beginTime = 0.0,
                          centerTime = FALSE,
                          endTime = 0.0,
                          resolution = 40.0,
                          fftLength = 0,
                          windowSize = 20.0,
                          windowShift = 5.0,
                          window = 'BLACKMAN',
                          order = 0,
                          preemphasis = -0.95,
                          deemphasize = TRUE,
                        toFile = TRUE,
                        explicitExt = "lps",
                        outputDirectory = NULL,
                        assertLossless = NULL,
                        logToFile = FALSE,
                        keepConverted=FALSE,
                        convertOverwrites=FALSE,
                        verbose = TRUE){

  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"lps",explicitExt)
  newTracknames <- "LPS[dB]"
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
    spectrumType = 'LPS',
    resolution = resolution,
    fftLength = as.integer(fftLength),
    windowSize = windowSize,
    windowShift = windowShift,
    window = window,
    effectiveLength = TRUE,
    order = as.integer(order),
    preemphasis = preemphasis,
    deemphasize = deemphasize,
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
attr(trk_lps_spectrum,"ext") <-  "lps"
attr(trk_lps_spectrum,"tracks") <-  c("LPS[dB]")
attr(trk_lps_spectrum,"outputType") <-  "SSFF"
attr(trk_lps_spectrum,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(trk_lps_spectrum,"suggestCaching") <-  FALSE

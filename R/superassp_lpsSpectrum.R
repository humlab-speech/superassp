##' Calculate Linear Prediction smoothed spectrum (Rcpp-optimized)
##'
##' Short-term spectral analysis of the signal in `listOfFiles`
##' using the Fast Fourier Transform and linear predictive smoothing.
##'
##' The results will be will be written to an SSFF formated file with the base
##' name of the input file and extension *.lps* in a track *LPS[dB]* which contains amplitudes (on a dB scale) of
##' all frequencies in the computed spectrum.
##'
##' @details The function is a re-write of the [wrassp::lpsSpectrum] function, but
##' with media pre-conversion, better checking of preconditions such as the
##' input file existence, structured logging, and the use of a more modern
##' framework for user feedback. This version includes Rcpp optimizations
##' for improved performance on large batches of files.
##'
##' The native file type of this function is "wav" files (in "pcm_s16le"
##' format), SUNs "au", NIST, or CSL formats (kay or NSP extension). Input
##' signal conversion, when needed, is done by
##' [libavcodec](https://ffmpeg.org/libavcodec.html) and the excellent [av::av_audio_convert]
##' wrapper function
##'
##' @note
##' This function takes some time to apply but also result in data in a relatively large matrix.
##' It is therefore not usually efficient to store intermediate results in a cache.
##' However, if the number of signals it will be applied to
##' is *very* large, then caching of results may be warranted.
##'
##' @inheritParams cssSpectrum
##' @param resolution = <freq>: set FFT length to the smallest value which
##' results in a frequency resolution of <freq> Hz or better (default: 40.0)
##' @param fftLength = <num>: set FFT length to <num> points (overrules default
##' and 'resolution' option)
##' @param order = <num>: set prediction order to <num> (default: sampling
##' rate in kHz + 3)
##' @param preemphasis = <val>: set pre-emphasis factor to <val> (default:
##' -0.95)
##' @param deemphasize (default: undo spectral tilt due to
##' pre-emphasis used in LP analysis, i.e. TRUE)
##' @param windowSize window size in milliseconds
##'
##' @return The number of successfully written files (if `toFile=TRUE`), or a vector of `AsspDataObj` objects (if `toFile=FALSE`).
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##'
##' @seealso \code{\link{dftSpectrum}}, \code{\link{cssSpectrum}}, \code{\link{cepstrum}};
##' all derived from *libassp* \insertCite{s5h}{superassp} spectrum function.
##'
##' @useDynLib superassp, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
##'
##' # calculate linear prediction smoothed spectrum
##' res <- lpsSpectrum(path2wav, toFile=FALSE)
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

lpsSpectrum <- function(listOfFiles = NULL,
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

  #### Fast-path: check if all files are native and need no conversion ####
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles))

  file_exts <- fast_file_ext(listOfFiles)
  is_native <- fast_is_native(file_exts, nativeFiletypes)
  needs_timewindow <- (beginTime != 0.0 | endTime != 0.0)

  # Fast path: all files native, no time windows
  if(all(is_native) && !any(needs_timewindow)) {
    if(verbose) {
      cli::cli_inform("Applying {.fun {funName}} to {cli::no(n_files)} recording{?s}")
    }

    # Direct call - no conversion needed
    externalRes <- Map(
      function(x, bt, et) {
        .External("performAssp", x,
                  fname = "spectrum",
                  beginTime = bt,
                  centerTime = centerTime,
                  endTime = et,
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
                  toFile = toFile,
                  explicitExt = explicitExt,
                  progressBar = NULL,
                  outputDirectory = outputDirectory,
                  PACKAGE = "superassp")
      },
      listOfFiles,
      beginTime,
      endTime
    )

    listOfFilesDF <- data.frame(
      audio = listOfFiles,
      dsp_input = listOfFiles,
      beginTime = beginTime,
      endTime = endTime,
      stringsAsFactors = FALSE
    )
    toClear <- character(0)

  } else {
    # Slow path: needs conversion or time windowing
    #### Input file conversion (Rcpp-optimized) ####
    listOfFiles_toClear <- convertInputMediaFiles(
      listOfFiles, beginTime, endTime, windowShift,
      nativeFiletypes, preferedFiletype, knownLossless,
      funName, keepConverted, verbose
    )

    listOfFilesDF <- listOfFiles_toClear[[1]]
    toClear <- listOfFiles_toClear[[3]]

    # Verify all files are in native format (using Rcpp)
    file_exts <- fast_file_ext(listOfFilesDF$dsp_input)
    if(!all(fast_is_native(file_exts, nativeFiletypes))) {
      cli::cli_abort("File conversion failed - non-native formats remain")
    }

    #### Application of DSP C function  ####

    if(verbose) {
      cli::cli_inform("Applying {.fun {funName}} to {cli::no(n_files)} recording{?s}")
    }

    # Process files with vectorized approach
    externalRes <- Map(
      function(x, bt, et) {
        .External("performAssp", x,
                  fname = "spectrum",
                  beginTime = bt,
                  centerTime = centerTime,
                  endTime = et,
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
                  toFile = toFile,
                  explicitExt = explicitExt,
                  progressBar = NULL,
                  outputDirectory = outputDirectory,
                  PACKAGE = "superassp")
      },
      listOfFilesDF$dsp_input,
      listOfFilesDF$beginTime,
      listOfFilesDF$endTime
    )
  }

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
attr(lpsSpectrum,"ext") <-  "lps"
attr(lpsSpectrum,"tracks") <-  c("LPS[dB]")
attr(lpsSpectrum,"outputType") <-  "SSFF"
attr(lpsSpectrum,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(lpsSpectrum,"suggestCaching") <-  FALSE

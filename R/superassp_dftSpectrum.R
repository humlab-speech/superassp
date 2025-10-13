##' Computes the Discrete Fourier Transform spectrum (Rcpp-optimized)
##'
##' Short-term spectral analysis of the signal in <listOfFiles>
##' using the Fast Fourier Transform. The default is to
##' calculate an unsmoothed narrow-band spectrum with the
##' size of the analysis window equal to the length of the
##' FFT. The output from the FFT will be converted to a
##' power spectrum in dB from 0 Hz up to and including the
##' Nyquist rate.
##'
##'
##' The results will be will be written to an SSFF formated file with the base
##' name of the input file and extension *.dft* in a track *DFT[dB]* which contains amplitudes (on a dB scale) of
##' all frequencies in the computed spectrum.
##'
##' @details The function is a re-write of the [wrassp::dftSpectrum] function, but
##' with media pre-conversion, better checking of preconditions such as the
##' input file existence, structured logging, and the use of a more modern
##' framework for user feedback. The *libassp* \insertCite{s5h}{superassp} C library code is used for
##' DSP. This version includes Rcpp optimizations for improved performance on large batches of files.
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
##' @inheritParams lpsSpectrum
##' @param bandwidth = <freq>: set the effective analysis bandwidth to <freq>
##' Hz (default: 0, yielding the smallest possible value given the length of
##' the FFT)
##'
##' @return The number of successfully written files (if `toFile=TRUE`), or a vector of `AsspDataObj` objects (if `toFile=FALSE`).
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##'
##' @seealso \code{\link{cssSpectrum}}, \code{\link{lpsSpectrum}}, \code{\link{cepstrum}};
##' all derived from *libassp* \insertCite{s5h}{superassp} spectrum function.
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
##' res <- dftSpectrum(path2wav, toFile=FALSE)
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
'dftSpectrum' <- function(listOfFiles = NULL,
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
                  resolution = resolution,
                  fftLength = as.integer(fftLength),
                  windowShift = windowShift,
                  window = window,
                  bandwidth = bandwidth,
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
                  resolution = resolution,
                  fftLength = as.integer(fftLength),
                  windowShift = windowShift,
                  window = window,
                  bandwidth = bandwidth,
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
attr(dftSpectrum,"ext") <-  "dft"
attr(dftSpectrum,"tracks") <-  c("DFT[dB]")
attr(dftSpectrum,"outputType") <-  "SSFF"
attr(dftSpectrum,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(dftSpectrum,"suggestCaching") <-  FALSE

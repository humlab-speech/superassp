##' Performs a short-term cepstral analysis of the signal (Rcpp-optimized)
##'
##' @description
##' This function performs a cepstral analysis \insertCite{Oppenheim.2004.10.1109/msp.2004.1328092,Childers.1977.10.1109/proc.1977.10747}{superassp} on the signals in  `listOfFiles`
##' using the Fast Fourier Transform. The number of
##' coefficients per output record will also equal the
##' FFT length / 2 + 1 (which means that it will not be mirrored).
##'
##'
##' The results will be will be written to an SSFF formated file with the base
##' name of the input file and extension *.cep* in a track *C[dB]* which contains amplitudes of
##' at each ½ Quefrency (in ms).
##'
##' @details The function is a re-write of the [wrassp::cepstrum] function, but
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
##' @inheritParams acfana
##' @param resolution = <freq>: set FFT length to the smallest value which
##' results in a frequency resolution of <freq> Hz or better (default: 40.0)
##' @param fftLength = <num>: set FFT length to <num> points (overrules default
##' and 'resolution' option)
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##'
##' @seealso \code{\link{dftSpectrum}}, \code{\link{cssSpectrum}}, \code{\link{lpsSpectrum}};
##' all derived from libassp's spectrum function
##'
##' @return The number of successfully written files (if `toFile=TRUE`), or a vector of `AsspDataObj` objects (if `toFile=FALSE`).
##'
##' @seealso [wrassp::cepstrum]
##' @seealso [superassp::AsspWindowTypes]
##' @seealso [av::av_audio_convert]
##'
##' @useDynLib superassp, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
##'
##' # calulate cepstrum
##' res <- cepstrum(path2wav, toFile=FALSE)
##'
##' # plot cepstral values at midpoint of signal
##' plot(y=res[["C[dB]"]][400,],
##'     x=seq(1,ncol(res[["C[dB]"]])),
##'     type='l',
##'     xlab='Quefrency (ms)',
##'     ylab='Amplitude (dB)')
##'
##' @export
##'
##' @references
##'   \insertAllCited{}
##'
cepstrum<- function(listOfFiles = NULL,
                       beginTime = 0.0,
                       centerTime = FALSE,
                       endTime = 0.0,
                       resolution = 40.0,
                       fftLength = 0,
                       windowShift = 5.0,
                       window = 'BLACKMAN',
                    toFile = TRUE,
                    explicitExt = "cep",
                    outputDirectory = NULL,
                    assertLossless = NULL,
                    logToFile = FALSE,
                    keepConverted=FALSE,
                    convertOverwrites=FALSE,
                    verbose = TRUE){

  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"cep",explicitExt)
  newTracknames <- "C[dB]"
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
                  spectrumType = 'CEP',
                  resolution = resolution,
                  fftLength = as.integer(fftLength),
                  windowShift = windowShift,
                  window = window,
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
    # New path: load-and-process pattern using av package
    if(verbose) {
      cli::cli_inform("Applying {.fun {funName}} to {cli::no(n_files)} recording{?s}")
    }

    # Use new load-and-process helper
    result <- processMediaFiles_LoadAndProcess(
      listOfFiles = listOfFiles,
      beginTime = beginTime,
      endTime = endTime,
      nativeFiletypes = nativeFiletypes,
      fname = "spectrum",
      toFile = toFile,
      verbose = verbose,
      centerTime = centerTime,
      spectrumType = 'CEP',
      resolution = resolution,
      fftLength = as.integer(fftLength),
      windowShift = windowShift,
      window = window,
      explicitExt = explicitExt,
      outputDirectory = outputDirectory
    )

    externalRes <- result$externalRes
    listOfFilesDF <- result$listOfFilesDF
    toClear <- character(0)  # No files to clean up with load-and-process
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
attr(cepstrum,"ext") <-  "cep"
attr(cepstrum,"tracks") <-  c("C[dB]")
attr(cepstrum,"outputType") <-  "SSFF"
attr(cepstrum,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(cepstrum,"suggestCaching") <-  FALSE

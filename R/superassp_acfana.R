##' Analysis of short-term autocorrelation function of signals (Rcpp-optimized)
##'
##' @description Applies the autocorrelation function to windows of the input
##' signals listed in `listOfFiles`. Input signals not in a file format natively
##' supported will be converted before the autocorrelation functions are
##' computed. The conversion process will display warnings about input files
##' that are not in known losslessly encoded formats.
##'
##' The results will be will be written to an SSFF formated file with the base
##' name of the input file and extension *.acf* in a track *ACF*.
##'
##' @details The function is a re-write of the [wrassp::acfana] function, but
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
##' This function is not considered computationally expensive enough to require caching of
##' results if applied to many signals. However, if the number of signals it will be applied to
##' is *very* large, then caching of results may be warranted.
##'
##' The autocorrelation function is reported as a dimensionless quantity, as it represents
##' the product of signal amplitudes at different time lags.
##'
##' @param listOfFiles vector of file paths to be processed by function
##' @param beginTime the time point (in seconds) of the start of the analysed
##'   interval. A NULL or 0 is interpreted as the start of the signal file.
##'   If a vector of time points is supplied, the length of that vector needs
##'   to correspond with the length of `listOfFiles`.
##' @param centerTime sets a single-frame analysis time point (in seconds).
##'   Overrides `beginTime`, `endTime` and `windowShift` parameters.
##' @param endTime the time point (in seconds) of the end of the analysed
##'   interval. A NULL or 0 is interpreted as the end of the signal file.
##'   If a vector of time points is supplied, the length of that vector needs
##'   to correspond with the length of `listOfFiles`.
##' @param windowShift the amount of time (in ms) that the analysis window will
##'   be shifted between analysis frames
##' @param windowSize the analysis window size (in ms); overrides the effect of
##'   the `effectiveLength` parameter
##' @param effectiveLength make window size effective rather than exact
##' @param window = the analysis window function type ("BLACKMAN" by default).
##'   See [superassp::AsspWindowTypes] for a list of supported window types.
##' @param analysisOrder the analysis order. The `NULL` or `0` sets the analysis
##'   order to the sample rate (in kHz) + 3, so that a signal with a 16000 Hz
##'   sampling rate will be analysed using an `analysisOrder` of 19.
##' @param energyNormalization calculate energy-normalized autocorrelation
##' @param lengthNormalization calculate length-normalized autocorrelation
##' @param toFile Should the function write the results to a file, with the
##'   (default) file extension (`TRUE`) or returned as a list of
##'   [AsspDataObj] objects (`FALSE`)?
##' @param explicitExt the file extension will be used when
##'   result files are written (`toFile=TRUE`), but the file extension can be
##'   set to something else using this function argument.
##' @param outputDirectory directory in which output files are stored. Defaults
##'   to NULL which means that the result file will be stored in the same
##'   directory as the input file.
##' @param verbose display verbose information about processing steps taken, as
##'   well as progress bars.
##' @param assertLossless an optional list of file extensions that the user wants to assert
##'   contains losslessly encoded signals data.
##' @param logToFile whether to log commands to a separate logfile in the
##'   `outputDirectory`. Logging will otherwise be in the function-specific logging
##'   namespace of [logger] and will be put wherever this namespace is defined to place its output.
##'   See [logger::log_appender] for details.
##' @param keepConverted whether to keep converted files
##' @param convertOverwrites whether conversion should overwrite existing files
##'
##' @return The number of successfully written files (if `toFile=TRUE`), or a vector of `AsspDataObj` objects (if `toFile=FALSE`).
##'
##' @seealso [wrassp::acfana]
##' @seealso [superassp::AsspWindowTypes]
##' @seealso [av::av_audio_convert]
##'
##' @useDynLib superassp, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @examples
##' # get path to audio file
#' path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
#'
#' # calculate short-term autocorrelation
#' res <- acfana(path2wav, toFile=FALSE)
#'
#' # plot short-term autocorrelation values
#' matplot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
#'         attr(res, 'startTime'),
#'         res$acf,
#'         type='l',
#'         xlab='time (s)',
#'         ylab='Short-term autocorrelation values')
##'
##' @export
##' @references
##'   \insertAllCited{}
##'

acfana <- function(listOfFiles = NULL,
                   beginTime = 0,
                   centerTime = FALSE,
                   endTime = 0,
                   windowShift = 5,
                   windowSize = 20,
                   effectiveLength = TRUE,
                   window = "BLACKMAN",
                   analysisOrder = 0,
                   energyNormalization = FALSE,
                   lengthNormalization = FALSE,
                   toFile = TRUE,
                   explicitExt = "acf",
                   outputDirectory = NULL,
                   assertLossless = NULL,
                   logToFile = FALSE,
                   keepConverted=FALSE,
                   convertOverwrites=FALSE,
                   verbose = TRUE) {



  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"acf",explicitExt)
  newTracknames <- "ACF"  ## Only used if SSFF tracks needs to be renamed from the called function (in C) before returning the SSFF track obj
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

  knownLossless <- c(assertLossless,knownLossless()) #Use the user asserted information about lossless encoding, in addition to what is already known by superassp

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
                  fname = "acfana",
                  beginTime = bt,
                  centerTime = centerTime,
                  endTime = et,
                  windowShift = windowShift,
                  windowSize = windowSize,
                  effectiveLength = effectiveLength,
                  window = window,
                  analysisOrder = as.integer(analysisOrder),
                  energyNormalization = energyNormalization,
                  lengthNormalization = lengthNormalization,
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
                  fname = "acfana",
                  beginTime = bt,
                  centerTime = centerTime,
                  endTime = et,
                  windowShift = windowShift,
                  windowSize = windowSize,
                  effectiveLength = effectiveLength,
                  window = window,
                  analysisOrder = as.integer(analysisOrder),
                  energyNormalization = energyNormalization,
                  lengthNormalization = lengthNormalization,
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

# Function attributes
attr(acfana,"ext") <-  "acf"
attr(acfana,"tracks") <-  c("ACF")
attr(acfana,"outputType") <-  "SSFF"
attr(acfana,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(acfana,"suggestCaching") <-  FALSE

##' Analysis of the averages of the short-term positive and negative zero-crossing rates (Rcpp-optimized)
##'
##' Analysis of the averages of the short-term positive and
##' negative zero-crossing rates of the signal in `listOfFiles` using the *libassp* C library
##'  \insertCite{s5h}{superassp} function. If `toFile` is `TRUE`, the results will be written to an output fil in the SSFF binary format, with the
##' same name as the input file, but with an extension *.zcr* and with a track named 'ZCR[Hz]'.
##'
##' Input signals not in a natively supported file format will be converted
##' before the autocorrelation functions are computed. The conversion process
##' will display warnings about input files that are not in known losslessly
##' encoded formats.
##'
##'
##' @details The function is a re-write of the [wrassp::zcrana] function, but
##' with media pre-conversion, better checking of preconditions such as the
##' input file existence, structured logging, and the use of a more modern
##' framework for user feedback. This version includes Rcpp optimizations
##' for improved performance on large batches of files.
##'
##' The native file type of this function is "wav" files (in "pcm_s16le"
##' format), SUNs "au", NIST, or CSL formats (kay or NSP extension). Input
##' signal conversion, when needed, is done by
##' [libavcodec](https://ffmpeg.org/libavcodec.html) and the excellent [av]
##' wrapper package.
##'
##' @note
##' This function is not considered computationally expensive enough to require caching of
##' results if applied to many signals. However, if the number of signals it will be applied to
##' is *very* long, then caching of results may be warranted.
##'
##' @inheritParams acfana
##' @param windowSize window size in milliseconds
##'
##' @return If `toFile` is `FALSE`, the function returns a list of [AsspDataObj]
##'   objects. If `toFile` is `TRUE`, the number (integer) of successfully
##'   processed and stored output files is returned.
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##'
##' @useDynLib superassp, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @examples
##'# get path to audio file
##'path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
##'
##'# calculate zcr values
##'res <- zcrana(path2wav, toFile=FALSE)
##'
##'# plot zcr values
##'plot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'       attr(res, 'startTime'),
##'     res[["ZCR[Hz]"]],
##'     type='l',
##'     xlab='time (s)',
##'     ylab='Zero Crossing Rates (Hz)')
##'
##' @export
zcrana <- function(listOfFiles = NULL,
                   beginTime = 0,
                   centerTime = FALSE,
                   endTime = 0,
                   windowShift = 5,
                   windowSize = 25,
                   toFile = TRUE,
                   explicitExt = "zcr",
                   outputDirectory = NULL,
                   assertLossless = NULL,
                   logToFile = FALSE,
                   convertOverwrites=FALSE,
                   keepConverted=FALSE,
                   verbose = TRUE){

  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"zcr",explicitExt)
  newTracknames <- c("ZCR[Hz]")
  nativeFiletypes <- c("wav","au","kay","nist","nsp")


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
    cli::cli_inform("Applying {.fun {funName}} to {cli::no(n_files)} recording{?s}")
  }

  # Use unified load-and-process helper (works for all file formats)
  result <- processMediaFiles_LoadAndProcess(
    listOfFiles = listOfFiles,
    beginTime = beginTime,
    endTime = endTime,
    nativeFiletypes = nativeFiletypes,
    fname = "zcrana",
    toFile = toFile,
    verbose = verbose,
    centerTime = centerTime,
    windowShift = windowShift,
    windowSize = windowSize,
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
attr(zcrana,"ext") <-  "zcr"
attr(zcrana,"tracks") <-  c("ZCR[Hz]")
attr(zcrana,"outputType") <-  "SSFF"
attr(zcrana,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(zcrana,"suggestCaching") <-  FALSE

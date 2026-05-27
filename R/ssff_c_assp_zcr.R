##' Track short-term zero-crossing rate
##'
##' Computes the average of the short-term positive and negative zero-crossing
##' rates of audio signals using the *libassp* C library
##' \insertCite{s5h}{superassp}. ZCR is a simple, fast measure correlated with
##' spectral centroid and useful for voicing detection and fricative
##' classification.
##'
##' @inheritParams trk_acf
##' @param windowSize Numeric. Analysis window size in milliseconds. Default 25 ms.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
##'   \describe{
##'     \item{\code{ZCR[Hz]}}{REAL32, Hz, n_frames x 1 column.
##'       Average zero-crossing rate (positive and negative crossings combined)
##'       per frame, expressed as a rate in Hz.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 200 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @details
##' The ZCR is reported in Hz (crossings per second), averaged over the positive
##' and negative zero-crossing rates within each analysis window.
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##'
##' @seealso [wrassp::zcrana]
##'
##' @useDynLib superassp, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @examples
##'# get path to audio file
##'path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
##'
##'# calculate zcr values
##'res <- trk_zcr(path2wav, toFile=FALSE)
##'
##'# plot zcr values
##'plot(seq(0, n_records(res) - 1) / sample_rate(res) +
##'       attr(res, 'startTime'),
##'     res[["ZCR[Hz]"]],
##'     type='l',
##'     xlab='time (s)',
##'     ylab='Zero Crossing Rates (Hz)')
##'
##' @export
trk_zcr <- function(listOfFiles,
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
    format_apply_msg(funName, n_files, beginTime, endTime)
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
attr(trk_zcr,"ext") <-  "zcr"
attr(trk_zcr,"tracks") <-  c("ZCR[Hz]")
attr(trk_zcr,"outputType") <-  "SSFF"
attr(trk_zcr,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(trk_zcr,"suggestCaching") <-  FALSE

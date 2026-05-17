##' Differentiate an audio waveform
##'
##' Applies a first-order finite-difference filter to audio signals using the
##' *libassp* C library \insertCite{s5h}{superassp}. Forward, backward, and
##' central difference modes are supported. Useful for pre-emphasising high
##' frequencies or computing the derivative of a waveform before further
##' analysis.
##'
##' @param listOfFiles Character vector of audio file paths. Any format supported by
##'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
##' @param computeBackwardDifference Logical. Use backward difference instead of
##'   forward. Default \code{FALSE}.
##' @param computeCentralDifference Logical. Use central difference instead of
##'   forward. Default \code{FALSE}.
##' @param channel Integer. Audio channel to process (1-based). Default \code{1L}.
##' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
##' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
##' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
##'   count written (invisibly). If \code{FALSE}, return an \code{AsspDataObj}.
##'   Default \code{TRUE}.
##' @param explicitExt Character. Output file extension. Default \code{"dif"}.
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
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track name
##'   preserved from libassp output (typically the same label as the input audio
##'   channel), containing INT16 or REAL32 differentiated sample values.
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @details
##' At most one of \code{computeBackwardDifference} or \code{computeCentralDifference}
##' should be \code{TRUE}. If both are \code{FALSE} (default), forward difference is used.
##'
##' @author Fredrik Nylén
##'
##' @export
##' @useDynLib superassp, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @references
##'   \insertAllCited{}
##'
##' @examples
##' path2wav <- list.files(system.file("samples", "sustained", package = "superassp"),
##'                        pattern = glob2rx("a1.wav"), full.names = TRUE)
##' res <- trk_afdiff(path2wav, toFile = FALSE)
trk_afdiff <- function(listOfFiles,
                       computeBackwardDifference = FALSE,
                       computeCentralDifference = FALSE,
                       channel = 1L,
                       beginTime = 0,
                       endTime = 0,
                       toFile = TRUE,
                       explicitExt = "dif",
                       outputDirectory = NULL,
                       assertLossless = NULL,
                       logToFile = FALSE,
                       keepConverted = FALSE,
                       convertOverwrites = FALSE,
                       verbose = TRUE) {

  explicitExt     <- if (is.null(explicitExt)) "dif" else explicitExt
  nativeFiletypes <- c("wav", "au", "kay", "nist", "nsp")

  currCall <- rlang::current_call()
  funName  <- rlang::call_name(currCall)

  knownLossless <- c(assertLossless, knownLossless())

  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime   <- if (is.null(endTime))   0.0 else endTime

  n_files <- length(listOfFiles)

  if (length(beginTime) > 1 && length(beginTime) != n_files) {
    cli::cli_abort("The {.par beginTime} must be length 1 or match {.par listOfFiles} length.")
  }
  if (length(endTime) > 1 && length(endTime) != n_files) {
    cli::cli_abort("The {.par endTime} must be length 1 or match {.par listOfFiles} length.")
  }

  beginTime <- fast_recycle_times(beginTime, n_files)
  endTime   <- fast_recycle_times(endTime,   n_files)

  makeOutputDirectory(outputDirectory, logToFile, funName)

  if (verbose) {
    format_apply_msg(funName, n_files, beginTime, endTime)
  }

  result <- processMediaFiles_LoadAndProcess(
    listOfFiles    = listOfFiles,
    beginTime      = beginTime,
    endTime        = endTime,
    nativeFiletypes = nativeFiletypes,
    fname          = "afdiff",
    toFile         = toFile,
    verbose        = verbose,
    computeBackwardDifference = computeBackwardDifference,
    computeCentralDifference  = computeCentralDifference,
    channel        = as.integer(channel),
    explicitExt    = explicitExt,
    outputDirectory = outputDirectory
  )

  externalRes <- result$externalRes
  toClear     <- character(0)

  if (n_files == 1) externalRes <- externalRes[[1]]

  cleanupConvertedInputMediaFiles(toClear, keepConverted, verbose)

  return(externalRes)
}

attr(trk_afdiff, "ext")             <- "dif"
attr(trk_afdiff, "tracks")          <- character(0)
attr(trk_afdiff, "outputType")      <- "SSFF"
attr(trk_afdiff, "nativeFiletypes") <- c("wav", "au", "kay", "nist", "nsp")
attr(trk_afdiff, "suggestCaching")  <- FALSE

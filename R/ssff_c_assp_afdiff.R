##' Differentiation of audio signals
##'
##' @description Computes finite differences of audio signals listed in
##' `listOfFiles` using algorithms implemented in *libassp*
##' \insertCite{s5h}{superassp}. Supports forward (default), backward, and
##' central difference methods. Input signals not in a natively supported
##' format are converted before processing; the conversion process will
##' display warnings about input files that are not in known losslessly
##' encoded formats.
##'
##' The results are written to an SSFF formatted file with the base name of
##' the input file and the extension specified by `explicitExt` (default:
##' `"dif"`).
##'
##' @details Forward, backward, and central difference options correspond to
##' first-order finite-difference approximations of the derivative. At most
##' one of `computeBackwardDifference` or `computeCentralDifference` should
##' be `TRUE`; if both are `FALSE` (the default), a forward difference is
##' computed.
##'
##' Native file types: WAV (`pcm_s16le`), Sun AU, NIST, CSL (kay/nsp).
##' Conversion via [libavcodec](https://ffmpeg.org/libavcodec.html) /
##' [av::av_audio_convert].
##'
##' @param listOfFiles vector of file paths to be processed
##' @param computeBackwardDifference compute backward difference instead of
##'   forward (default: `FALSE`)
##' @param computeCentralDifference compute central difference instead of
##'   forward (default: `FALSE`)
##' @param channel audio channel to process (1-based, default: `1L`)
##' @param beginTime start of processed interval in seconds (0 = file start)
##' @param endTime end of processed interval in seconds (0 = file end)
##' @param toFile write results to file (`TRUE`) or return `AsspDataObj`
##'   (`FALSE`)
##' @param explicitExt output file extension (default: `"dif"`)
##' @param outputDirectory directory for output files (NULL = same as input)
##' @param assertLossless additional file extensions to treat as lossless
##' @param logToFile write log to `outputDirectory` instead of console
##' @param keepConverted keep intermediate converted files
##' @param convertOverwrites allow conversion to overwrite existing files
##' @param verbose display progress messages
##'
##' @return Number of files written (`toFile=TRUE`) or an `AsspDataObj` /
##'   list thereof (`toFile=FALSE`).
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

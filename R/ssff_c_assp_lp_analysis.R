.lp_ana_impl <- function(listOfFiles,
                         beginTime = 0.0,
                         centerTime = FALSE,
                         endTime = 0.0,
                         windowShift = 5.0,
                         windowSize = 20.0,
                         effectiveLength = TRUE,
                         window = "BLACKMAN",
                         analysisOrder = NULL,
                         preemphasis = -0.95,
                         toFile = FALSE,
                         explicitExt = NULL,
                         outputDirectory = NULL,
                         assertLossless = NULL,
                         logToFile = FALSE,
                         keepConverted = FALSE,
                         convertOverwrites = FALSE,
                         verbose = TRUE,
                         lpType,
                         fileExt,
                         newTracknames) {

  explicitExt <- ifelse(is.null(explicitExt), fileExt, explicitExt)
  nativeFiletypes <- c("wav", "au", "kay", "nist", "nsp")

  if (!isAsspWindowType(toupper(window))) {
    cli::cli_abort(c(
      "WindowFunction of type {.val {window}} is not supported!",
      "i" = "Accepted window types for routines implemented in *libassp* are {.field {AsspWindowTypes()}}."
    ))
  }

  funName <- paste0(fileExt, "ana")

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
    fname          = "rfcana",
    toFile         = toFile,
    verbose        = verbose,
    centerTime     = centerTime,
    windowShift    = windowShift,
    windowSize     = windowSize,
    effectiveLength = effectiveLength,
    window         = window,
    order          = ifelse(is.null(analysisOrder), 0L, as.integer(analysisOrder)),
    preemphasis    = preemphasis,
    lpType         = lpType,
    explicitExt    = explicitExt,
    outputDirectory = outputDirectory
  )

  externalRes <- result$externalRes
  toClear     <- character(0)

  if (!toFile && !is.null(newTracknames)) {
    n_tracks <- length(names(externalRes[[1]]))
    if (n_tracks != length(newTracknames)) {
      cli::cli_abort(c(
        "Wrong number of track names supplied:",
        "i" = "Track{?s} named: {.field {names(externalRes[[1]])}}"
      ))
    }
    externalRes <- fast_rename_tracks(externalRes, newTracknames)
  }

  if (n_files == 1) externalRes <- externalRes[[1]]

  cleanupConvertedInputMediaFiles(toClear, keepConverted, verbose)

  return(externalRes)
}


##' Track LP reflection coefficients
##'
##' Linear Prediction analysis of audio signals using the autocorrelation method
##' and Durbin recursion, implemented in the *libassp* C library
##' \insertCite{s5h}{superassp}. Returns per-frame RMS amplitudes and reflection
##' coefficients. Use \code{trk_rfc} when a lattice-filter or PARCOR
##' representation of the vocal tract is needed.
##'
##' @usage trk_rfc(listOfFiles = NULL,
##'   beginTime = 0.0,
##'   centerTime = FALSE,
##'   endTime = 0.0,
##'   windowShift = 5.0,
##'   windowSize = 20.0,
##'   effectiveLength = TRUE,
##'   window = 'BLACKMAN',
##'   analysisOrder = NULL,
##'   preemphasis = -0.95,
##'   toFile = FALSE,
##'   explicitExt = NULL,
##'   outputDirectory = NULL,
##'   assertLossless = NULL,
##'   logToFile = FALSE,
##'   keepConverted = FALSE,
##'   convertOverwrites = FALSE,
##'   verbose = TRUE)
##'
##' @param listOfFiles Character vector of audio file paths. Any format supported by
##'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
##' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
##' @param centerTime Numeric or logical. Single-frame analysis time point in seconds;
##'   overrides \code{beginTime}, \code{endTime}, and \code{windowShift}. Default \code{FALSE}.
##' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
##' @param windowShift Numeric. Frame shift in milliseconds; sets output frame rate
##'   (1000 / windowShift Hz). Default 5 ms.
##' @param windowSize Numeric. Analysis window size in milliseconds. Default 20 ms.
##' @param effectiveLength Logical. Make window size effective rather than exact.
##'   Default \code{FALSE}.
##' @param window Character. Analysis window function type. Default \code{"BLACKMAN"}.
##'   See [superassp::AsspWindowTypes].
##' @param analysisOrder Integer. LP order; \code{NULL} or 0 defaults to
##'   sample rate in kHz + 3. Default \code{NULL}.
##' @param preemphasis Numeric. Pre-emphasis factor. Default -0.95.
##' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
##'   count written (invisibly). If \code{FALSE}, return an \code{AsspDataObj}.
##'   Default \code{FALSE}.
##' @param explicitExt Character. Output file extension override. Default \code{NULL}
##'   (uses \code{"rfc"}).
##' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
##'   writes alongside the input file.
##' @param assertLossless Character vector of additional file extensions to treat as
##'   losslessly encoded.
##' @param logToFile Logical. Write processing log to a file in \code{outputDirectory}
##'   rather than the console. Default \code{FALSE}.
##' @param keepConverted Logical. Retain intermediate transcoded files. Default \code{FALSE}.
##' @param convertOverwrites Logical. Allow transcoding to overwrite existing files.
##'   Default \code{FALSE}.
##' @param verbose Logical. Print per-file progress. Default \code{FALSE}.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
##'   \describe{
##'     \item{\code{RMS[dB]}}{REAL32, dB, n_frames x 1. RMS amplitude of the input frame.}
##'     \item{\code{gain[dB]}}{REAL32, dB, n_frames x 1. RMS amplitude of the LP residual.}
##'     \item{\code{RFC}}{REAL32, dimensionless, n_frames x \code{analysisOrder} columns.
##'       Reflection (PARCOR) coefficients, one per LP order.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 200 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @details
##' \code{analysisOrder = NULL} selects an order equal to sample rate in kHz + 3.
##' Pre-emphasis is applied before LP analysis and affects the residual gain.
##'
##' @seealso [wrassp::rfcana] [superassp::AsspWindowTypes] [av::av_audio_convert]
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
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
##' res <- trk_rfc(path2wav, toFile = FALSE)
##' matplot(seq(0, n_records(res) - 1) / sample_rate(res) +
##'           attr(res, "startTime"),
##'         res$RFC, type = "l",
##'         xlab = "time (s)", ylab = "Reflection coefficient values")
trk_rfc <- function(listOfFiles,
                       beginTime = 0.0,
                       centerTime = FALSE,
                       endTime = 0.0,
                       windowShift = 5.0,
                       windowSize = 20.0,
                       effectiveLength = TRUE,
                       window = "BLACKMAN",
                       analysisOrder = NULL,
                       preemphasis = -0.95,
                       toFile = FALSE,
                       explicitExt = NULL,
                       outputDirectory = NULL,
                       assertLossless = NULL,
                       logToFile = FALSE,
                       keepConverted = FALSE,
                       convertOverwrites = FALSE,
                       verbose = TRUE) {
  .lp_ana_impl(listOfFiles = listOfFiles,
               beginTime = beginTime, centerTime = centerTime,
               endTime = endTime, windowShift = windowShift,
               windowSize = windowSize, effectiveLength = effectiveLength,
               window = window, analysisOrder = analysisOrder,
               preemphasis = preemphasis, toFile = toFile,
               explicitExt = explicitExt, outputDirectory = outputDirectory,
               assertLossless = assertLossless, logToFile = logToFile,
               keepConverted = keepConverted,
               convertOverwrites = convertOverwrites, verbose = verbose,
               lpType = "RFC", fileExt = "rfc",
               newTracknames = c("RMS[dB]", "gain[dB]", "RFC"))
}
attr(trk_rfc, "ext")             <- "rfc"
attr(trk_rfc, "tracks")          <- c("RMS[dB]", "gain[dB]", "RFC")
attr(trk_rfc, "outputType")      <- "SSFF"
attr(trk_rfc, "nativeFiletypes") <- c("wav", "au", "kay", "nist", "nsp")
attr(trk_rfc, "suggestCaching")  <- FALSE


##' Track LP-derived vocal tract area function coefficients
##'
##' Linear Prediction analysis of audio signals using the autocorrelation method
##' and Durbin recursion, implemented in the *libassp* C library
##' \insertCite{s5h}{superassp}. Returns per-frame RMS amplitudes and vocal
##' tract area function (ARF) coefficients derived from the reflection
##' coefficients. Useful for vocal tract modelling applications.
##'
##' @usage trk_arf(listOfFiles = NULL,
##'   beginTime = 0.0,
##'   centerTime = FALSE,
##'   endTime = 0.0,
##'   windowShift = 5.0,
##'   windowSize = 20.0,
##'   effectiveLength = TRUE,
##'   window = 'BLACKMAN',
##'   analysisOrder = NULL,
##'   preemphasis = -0.95,
##'   toFile = FALSE,
##'   explicitExt = NULL,
##'   outputDirectory = NULL,
##'   assertLossless = NULL,
##'   logToFile = FALSE,
##'   keepConverted = FALSE,
##'   convertOverwrites = FALSE,
##'   verbose = TRUE)
##'
##' @inheritParams trk_rfc
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
##'   \describe{
##'     \item{\code{RMS[dB]}}{REAL32, dB, n_frames x 1. RMS amplitude of the input frame.}
##'     \item{\code{gain[dB]}}{REAL32, dB, n_frames x 1. RMS amplitude of the LP residual.}
##'     \item{\code{ARF}}{REAL32, dimensionless area ratios, n_frames x \code{analysisOrder}
##'       columns. Area function coefficients derived from reflection coefficients.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 200 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @details
##' ARF coefficients are derived from the LP reflection coefficients via the
##' standard area ratio transformation. See \code{trk_rfc} for parameter details.
##'
##' @seealso [wrassp::rfcana] [superassp::AsspWindowTypes] [av::av_audio_convert]
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
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
##' res <- trk_arf(path2wav, toFile = FALSE)
##' matplot(seq(0, n_records(res) - 1) / sample_rate(res) +
##'           attr(res, "startTime"),
##'         res$ARF, type = "l",
##'         xlab = "time (s)", ylab = "Area function")
trk_arf <- function(listOfFiles,
                       beginTime = 0.0,
                       centerTime = FALSE,
                       endTime = 0.0,
                       windowShift = 5.0,
                       windowSize = 20.0,
                       effectiveLength = TRUE,
                       window = "BLACKMAN",
                       analysisOrder = NULL,
                       preemphasis = -0.95,
                       toFile = FALSE,
                       explicitExt = NULL,
                       outputDirectory = NULL,
                       assertLossless = NULL,
                       logToFile = FALSE,
                       keepConverted = FALSE,
                       convertOverwrites = FALSE,
                       verbose = TRUE) {
  .lp_ana_impl(listOfFiles = listOfFiles,
               beginTime = beginTime, centerTime = centerTime,
               endTime = endTime, windowShift = windowShift,
               windowSize = windowSize, effectiveLength = effectiveLength,
               window = window, analysisOrder = analysisOrder,
               preemphasis = preemphasis, toFile = toFile,
               explicitExt = explicitExt, outputDirectory = outputDirectory,
               assertLossless = assertLossless, logToFile = logToFile,
               keepConverted = keepConverted,
               convertOverwrites = convertOverwrites, verbose = verbose,
               lpType = "ARF", fileExt = "arf",
               newTracknames = c("RMS[dB]", "gain[dB]", "ARF"))
}
attr(trk_arf, "ext")             <- "arf"
attr(trk_arf, "tracks")          <- c("RMS[dB]", "gain[dB]", "ARF")
attr(trk_arf, "outputType")      <- "SSFF"
attr(trk_arf, "nativeFiletypes") <- c("wav", "au", "kay", "nist", "nsp")
attr(trk_arf, "suggestCaching")  <- FALSE


##' Track LP-derived log area ratios
##'
##' Linear Prediction analysis of audio signals using the autocorrelation method
##' and Durbin recursion, implemented in the *libassp* C library
##' \insertCite{s5h}{superassp}. Returns per-frame RMS amplitudes and log area
##' ratio (LAR) coefficients. LARs are a log-domain reparameterisation of
##' reflection coefficients that can be more numerically stable near the unit
##' circle.
##'
##' @usage trk_lar(listOfFiles = NULL,
##'   beginTime = 0.0,
##'   centerTime = FALSE,
##'   endTime = 0.0,
##'   windowShift = 5.0,
##'   windowSize = 20.0,
##'   effectiveLength = TRUE,
##'   window = 'BLACKMAN',
##'   analysisOrder = NULL,
##'   preemphasis = -0.95,
##'   toFile = FALSE,
##'   explicitExt = NULL,
##'   outputDirectory = NULL,
##'   assertLossless = NULL,
##'   logToFile = FALSE,
##'   keepConverted = FALSE,
##'   convertOverwrites = FALSE,
##'   verbose = TRUE)
##'
##' @inheritParams trk_rfc
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
##'   \describe{
##'     \item{\code{RMS[dB]}}{REAL32, dB, n_frames x 1. RMS amplitude of the input frame.}
##'     \item{\code{gain[dB]}}{REAL32, dB, n_frames x 1. RMS amplitude of the LP residual.}
##'     \item{\code{LAR}}{REAL32, dimensionless, n_frames x \code{analysisOrder} columns.
##'       Log area ratios derived from reflection coefficients.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 200 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @details
##' LAR coefficients are defined as log((1 + k) / (1 - k)) where k is the
##' reflection coefficient. See \code{trk_rfc} for parameter details.
##'
##' @seealso [wrassp::rfcana] [superassp::AsspWindowTypes] [av::av_audio_convert]
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
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
##' res <- trk_lar(path2wav, toFile = FALSE)
##' matplot(seq(0, n_records(res) - 1) / sample_rate(res) +
##'           attr(res, "startTime"),
##'         res$LAR, type = "l",
##'         xlab = "time (s)", ylab = "Log area ratios")
trk_lar <- function(listOfFiles,
                       beginTime = 0.0,
                       centerTime = FALSE,
                       endTime = 0.0,
                       windowShift = 5.0,
                       windowSize = 20.0,
                       effectiveLength = TRUE,
                       window = "BLACKMAN",
                       analysisOrder = NULL,
                       preemphasis = -0.95,
                       toFile = FALSE,
                       explicitExt = NULL,
                       outputDirectory = NULL,
                       assertLossless = NULL,
                       logToFile = FALSE,
                       keepConverted = FALSE,
                       convertOverwrites = FALSE,
                       verbose = TRUE) {
  .lp_ana_impl(listOfFiles = listOfFiles,
               beginTime = beginTime, centerTime = centerTime,
               endTime = endTime, windowShift = windowShift,
               windowSize = windowSize, effectiveLength = effectiveLength,
               window = window, analysisOrder = analysisOrder,
               preemphasis = preemphasis, toFile = toFile,
               explicitExt = explicitExt, outputDirectory = outputDirectory,
               assertLossless = assertLossless, logToFile = logToFile,
               keepConverted = keepConverted,
               convertOverwrites = convertOverwrites, verbose = verbose,
               lpType = "LAR", fileExt = "lar",
               newTracknames = c("RMS[dB]", "gain[dB]", "LAR"))
}
attr(trk_lar, "ext")             <- "lar"
attr(trk_lar, "tracks")          <- c("RMS[dB]", "gain[dB]", "LAR")
attr(trk_lar, "outputType")      <- "SSFF"
attr(trk_lar, "nativeFiletypes") <- c("wav", "au", "kay", "nist", "nsp")
attr(trk_lar, "suggestCaching")  <- FALSE


##' Track LP filter coefficients
##'
##' Linear Prediction analysis of audio signals using the autocorrelation method
##' and Durbin recursion, implemented in the *libassp* C library
##' \insertCite{s5h}{superassp}. Returns per-frame RMS amplitudes and direct-form
##' LP filter (predictor) coefficients. Use \code{trk_lpc} when direct-form LP
##' coefficients are needed for synthesis or spectral estimation.
##'
##' @usage trk_lpc(listOfFiles = NULL,
##'   beginTime = 0.0,
##'   centerTime = FALSE,
##'   endTime = 0.0,
##'   windowShift = 5.0,
##'   windowSize = 20.0,
##'   effectiveLength = TRUE,
##'   window = 'BLACKMAN',
##'   analysisOrder = NULL,
##'   preemphasis = -0.95,
##'   toFile = FALSE,
##'   explicitExt = NULL,
##'   outputDirectory = NULL,
##'   assertLossless = NULL,
##'   logToFile = FALSE,
##'   keepConverted = FALSE,
##'   convertOverwrites = FALSE,
##'   verbose = TRUE)
##'
##' @inheritParams trk_rfc
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
##'   \describe{
##'     \item{\code{RMS[dB]}}{REAL32, dB, n_frames x 1. RMS amplitude of the input frame.}
##'     \item{\code{gain[dB]}}{REAL32, dB, n_frames x 1. RMS amplitude of the LP residual.}
##'     \item{\code{LPC}}{REAL32, dimensionless, n_frames x \code{analysisOrder} columns.
##'       Direct-form LP predictor coefficients a_1 … a_p.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 200 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @details
##' Coefficients are the direct-form LP predictor coefficients (a_1 … a_p) from
##' the Durbin recursion. See \code{trk_rfc} for parameter details. The LPC
##' spectrum can be evaluated by \code{trk_lps_spectrum}.
##'
##' @seealso [wrassp::rfcana] [superassp::AsspWindowTypes] [av::av_audio_convert]
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
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
##' res <- trk_lpc(path2wav, toFile = FALSE)
##' matplot(seq(0, n_records(res) - 1) / sample_rate(res) +
##'           attr(res, "startTime"),
##'         res$LPC, type = "l",
##'         xlab = "time (s)", ylab = "LP filter coefficients")
trk_lpc <- function(listOfFiles,
                       beginTime = 0.0,
                       centerTime = FALSE,
                       endTime = 0.0,
                       windowShift = 5.0,
                       windowSize = 20.0,
                       effectiveLength = TRUE,
                       window = "BLACKMAN",
                       analysisOrder = NULL,
                       preemphasis = -0.95,
                       toFile = FALSE,
                       explicitExt = NULL,
                       outputDirectory = NULL,
                       assertLossless = NULL,
                       logToFile = FALSE,
                       keepConverted = FALSE,
                       convertOverwrites = FALSE,
                       verbose = TRUE) {
  .lp_ana_impl(listOfFiles = listOfFiles,
               beginTime = beginTime, centerTime = centerTime,
               endTime = endTime, windowShift = windowShift,
               windowSize = windowSize, effectiveLength = effectiveLength,
               window = window, analysisOrder = analysisOrder,
               preemphasis = preemphasis, toFile = toFile,
               explicitExt = explicitExt, outputDirectory = outputDirectory,
               assertLossless = assertLossless, logToFile = logToFile,
               keepConverted = keepConverted,
               convertOverwrites = convertOverwrites, verbose = verbose,
               lpType = "LPC", fileExt = "lpc",
               newTracknames = c("RMS[dB]", "gain[dB]", "LPC"))
}
attr(trk_lpc, "ext")             <- "lpc"
attr(trk_lpc, "tracks")          <- c("RMS[dB]", "gain[dB]", "LPC")
attr(trk_lpc, "outputType")      <- "SSFF"
attr(trk_lpc, "nativeFiletypes") <- c("wav", "au", "kay", "nist", "nsp")
attr(trk_lpc, "suggestCaching")  <- FALSE

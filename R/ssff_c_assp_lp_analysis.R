.lp_ana_impl <- function(listOfFiles = NULL,
                         beginTime = 0.0,
                         centerTime = FALSE,
                         endTime = 0.0,
                         windowShift = 5.0,
                         windowSize = 20.0,
                         effectiveLength = TRUE,
                         window = "BLACKMAN",
                         analysisOrder = NULL,
                         preemphasis = -0.95,
                         toFile = TRUE,
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


##' Linear Prediction analysis with reflection coefficient output
##'
##' Linear Prediction analysis of `listOfFiles` using the
##' autocorrelation method and the Durbin recursion.
##' Calculates the RMS amplitudes of the input and residual signal in dB,
##' and reflection coefficients using algorithms implemented in
##' *libassp* \insertCite{s5h}{superassp}.
##' Input signals not in a natively supported format are converted before
##' analysis. The conversion process will display warnings about input files
##' that are not in known losslessly encoded formats.
##'
##' The results are written to an SSFF formatted file with the base
##' name of the input file and extension `.rfc` with tracks
##' `RMS[dB]`, `gain[dB]`, and `RFC`.
##'
##' @usage trk_rfcana(listOfFiles = NULL,
##'   beginTime = 0.0,
##'   centerTime = FALSE,
##'   endTime = 0.0,
##'   windowShift = 5.0,
##'   windowSize = 20.0,
##'   effectiveLength = TRUE,
##'   window = 'BLACKMAN',
##'   analysisOrder = NULL,
##'   preemphasis = -0.95,
##'   toFile = TRUE,
##'   explicitExt = NULL,
##'   outputDirectory = NULL,
##'   assertLossless = NULL,
##'   logToFile = FALSE,
##'   keepConverted = FALSE,
##'   convertOverwrites = FALSE,
##'   verbose = TRUE)
##'
##' @details Re-write of [wrassp::rfcana] with media pre-conversion, structured
##' logging, and Rcpp optimisations for large batches.
##'
##' Native file types: WAV (`pcm_s16le`), Sun AU, NIST, CSL (kay/nsp).
##' Conversion via [libavcodec](https://ffmpeg.org/libavcodec.html) /
##' [av::av_audio_convert].
##'
##' @param listOfFiles vector of file paths to be processed by function
##' @param beginTime start of analysed interval in seconds (0 = file start)
##' @param centerTime single-frame analysis time point (s); overrides
##'   `beginTime`, `endTime`, and `windowShift`
##' @param endTime end of analysed interval in seconds (0 = file end)
##' @param windowShift frame shift in ms
##' @param windowSize analysis window size in ms
##' @param effectiveLength make window size effective rather than exact
##' @param window analysis window function type (see [superassp::AsspWindowTypes])
##' @param analysisOrder LP order; NULL/0 defaults to sampleRate(kHz) + 3
##' @param preemphasis pre-emphasis factor (default: -0.95)
##' @param toFile write results to file (`TRUE`) or return `AsspDataObj` (`FALSE`)
##' @param explicitExt override the default output file extension
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
##' res <- trk_rfcana(path2wav, toFile = FALSE)
##' matplot(seq(0, numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'           attr(res, "startTime"),
##'         res$RFC, type = "l",
##'         xlab = "time (s)", ylab = "Reflection coefficient values")
trk_rfcana <- function(listOfFiles = NULL,
                       beginTime = 0.0,
                       centerTime = FALSE,
                       endTime = 0.0,
                       windowShift = 5.0,
                       windowSize = 20.0,
                       effectiveLength = TRUE,
                       window = "BLACKMAN",
                       analysisOrder = NULL,
                       preemphasis = -0.95,
                       toFile = TRUE,
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
attr(trk_rfcana, "ext")             <- "rfc"
attr(trk_rfcana, "tracks")          <- c("RMS[dB]", "gain[dB]", "RFC")
attr(trk_rfcana, "outputType")      <- "SSFF"
attr(trk_rfcana, "nativeFiletypes") <- c("wav", "au", "kay", "nist", "nsp")
attr(trk_rfcana, "suggestCaching")  <- FALSE


##' Linear Prediction analysis with area function output
##'
##' Linear Prediction analysis of `listOfFiles` using the
##' autocorrelation method and the Durbin recursion.
##' Calculates the RMS amplitudes of the input and residual signal in dB,
##' and area function coefficients using algorithms implemented in
##' *libassp* \insertCite{s5h}{superassp}.
##' Input signals not in a natively supported format are converted before
##' analysis. The conversion process will display warnings about input files
##' that are not in known losslessly encoded formats.
##'
##' The results are written to an SSFF formatted file with the base
##' name of the input file and extension `.arf` with tracks
##' `RMS[dB]`, `gain[dB]`, and `ARF`.
##'
##' @usage trk_arfana(listOfFiles = NULL,
##'   beginTime = 0.0,
##'   centerTime = FALSE,
##'   endTime = 0.0,
##'   windowShift = 5.0,
##'   windowSize = 20.0,
##'   effectiveLength = TRUE,
##'   window = 'BLACKMAN',
##'   analysisOrder = NULL,
##'   preemphasis = -0.95,
##'   toFile = TRUE,
##'   explicitExt = NULL,
##'   outputDirectory = NULL,
##'   assertLossless = NULL,
##'   logToFile = FALSE,
##'   keepConverted = FALSE,
##'   convertOverwrites = FALSE,
##'   verbose = TRUE)
##'
##' @details Re-write of [wrassp::rfcana] with area function output,
##' media pre-conversion, structured logging, and Rcpp optimisations.
##'
##' @inheritParams trk_rfcana
##'
##' @return Number of files written (`toFile=TRUE`) or an `AsspDataObj` /
##'   list thereof (`toFile=FALSE`).
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
##' res <- trk_arfana(path2wav, toFile = FALSE)
##' matplot(seq(0, numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'           attr(res, "startTime"),
##'         res$ARF, type = "l",
##'         xlab = "time (s)", ylab = "Area function")
trk_arfana <- function(listOfFiles = NULL,
                       beginTime = 0.0,
                       centerTime = FALSE,
                       endTime = 0.0,
                       windowShift = 5.0,
                       windowSize = 20.0,
                       effectiveLength = TRUE,
                       window = "BLACKMAN",
                       analysisOrder = NULL,
                       preemphasis = -0.95,
                       toFile = TRUE,
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
attr(trk_arfana, "ext")             <- "arf"
attr(trk_arfana, "tracks")          <- c("RMS[dB]", "gain[dB]", "ARF")
attr(trk_arfana, "outputType")      <- "SSFF"
attr(trk_arfana, "nativeFiletypes") <- c("wav", "au", "kay", "nist", "nsp")
attr(trk_arfana, "suggestCaching")  <- FALSE


##' Linear Prediction analysis with log area ratio output
##'
##' Linear Prediction analysis of `listOfFiles` using the
##' autocorrelation method and the Durbin recursion.
##' Calculates the RMS amplitudes of the input and residual signal in dB,
##' and log area ratios using algorithms implemented in
##' *libassp* \insertCite{s5h}{superassp}.
##' Input signals not in a natively supported format are converted before
##' analysis. The conversion process will display warnings about input files
##' that are not in known losslessly encoded formats.
##'
##' The results are written to an SSFF formatted file with the base
##' name of the input file and extension `.lar` with tracks
##' `RMS[dB]`, `gain[dB]`, and `LAR`.
##'
##' @usage trk_larana(listOfFiles = NULL,
##'   beginTime = 0.0,
##'   centerTime = FALSE,
##'   endTime = 0.0,
##'   windowShift = 5.0,
##'   windowSize = 20.0,
##'   effectiveLength = TRUE,
##'   window = 'BLACKMAN',
##'   analysisOrder = NULL,
##'   preemphasis = -0.95,
##'   toFile = TRUE,
##'   explicitExt = NULL,
##'   outputDirectory = NULL,
##'   assertLossless = NULL,
##'   logToFile = FALSE,
##'   keepConverted = FALSE,
##'   convertOverwrites = FALSE,
##'   verbose = TRUE)
##'
##' @details Re-write of [wrassp::rfcana] with log area ratio output,
##' media pre-conversion, structured logging, and Rcpp optimisations.
##'
##' @inheritParams trk_rfcana
##'
##' @return Number of files written (`toFile=TRUE`) or an `AsspDataObj` /
##'   list thereof (`toFile=FALSE`).
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
##' res <- trk_larana(path2wav, toFile = FALSE)
##' matplot(seq(0, numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'           attr(res, "startTime"),
##'         res$LAR, type = "l",
##'         xlab = "time (s)", ylab = "Log area ratios")
trk_larana <- function(listOfFiles = NULL,
                       beginTime = 0.0,
                       centerTime = FALSE,
                       endTime = 0.0,
                       windowShift = 5.0,
                       windowSize = 20.0,
                       effectiveLength = TRUE,
                       window = "BLACKMAN",
                       analysisOrder = NULL,
                       preemphasis = -0.95,
                       toFile = TRUE,
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
attr(trk_larana, "ext")             <- "lar"
attr(trk_larana, "tracks")          <- c("RMS[dB]", "gain[dB]", "LAR")
attr(trk_larana, "outputType")      <- "SSFF"
attr(trk_larana, "nativeFiletypes") <- c("wav", "au", "kay", "nist", "nsp")
attr(trk_larana, "suggestCaching")  <- FALSE


##' Linear Prediction analysis with LP filter coefficient output
##'
##' Linear Prediction analysis of `listOfFiles` using the
##' autocorrelation method and the Durbin recursion.
##' Calculates the RMS amplitudes of the input and residual signal in dB,
##' and LP filter coefficients using algorithms implemented in
##' *libassp* \insertCite{s5h}{superassp}.
##' Input signals not in a natively supported format are converted before
##' analysis. The conversion process will display warnings about input files
##' that are not in known losslessly encoded formats.
##'
##' The results are written to an SSFF formatted file with the base
##' name of the input file and extension `.lpc` with tracks
##' `RMS[dB]`, `gain[dB]`, and `LPC`.
##'
##' @usage trk_lpcana(listOfFiles = NULL,
##'   beginTime = 0.0,
##'   centerTime = FALSE,
##'   endTime = 0.0,
##'   windowShift = 5.0,
##'   windowSize = 20.0,
##'   effectiveLength = TRUE,
##'   window = 'BLACKMAN',
##'   analysisOrder = NULL,
##'   preemphasis = -0.95,
##'   toFile = TRUE,
##'   explicitExt = NULL,
##'   outputDirectory = NULL,
##'   assertLossless = NULL,
##'   logToFile = FALSE,
##'   keepConverted = FALSE,
##'   convertOverwrites = FALSE,
##'   verbose = TRUE)
##'
##' @details Re-write of [wrassp::rfcana] with LP filter coefficient output,
##' media pre-conversion, structured logging, and Rcpp optimisations.
##'
##' @inheritParams trk_rfcana
##'
##' @return Number of files written (`toFile=TRUE`) or an `AsspDataObj` /
##'   list thereof (`toFile=FALSE`).
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
##' res <- trk_lpcana(path2wav, toFile = FALSE)
##' matplot(seq(0, numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'           attr(res, "startTime"),
##'         res$LPC, type = "l",
##'         xlab = "time (s)", ylab = "LP filter coefficients")
trk_lpcana <- function(listOfFiles = NULL,
                       beginTime = 0.0,
                       centerTime = FALSE,
                       endTime = 0.0,
                       windowShift = 5.0,
                       windowSize = 20.0,
                       effectiveLength = TRUE,
                       window = "BLACKMAN",
                       analysisOrder = NULL,
                       preemphasis = -0.95,
                       toFile = TRUE,
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
attr(trk_lpcana, "ext")             <- "lpc"
attr(trk_lpcana, "tracks")          <- c("RMS[dB]", "gain[dB]", "LPC")
attr(trk_lpcana, "outputType")      <- "SSFF"
attr(trk_lpcana, "nativeFiletypes") <- c("wav", "au", "kay", "nist", "nsp")
attr(trk_lpcana, "suggestCaching")  <- FALSE

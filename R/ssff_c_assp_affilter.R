##' Apply a digital filter to audio signals
##'
##' Filters audio waveforms using FIR or Butterworth IIR designs implemented in
##' the *libassp* C library \insertCite{s5h}{superassp}. Supports high-pass,
##' low-pass, and band-pass configurations. At least one of \code{highPass} or
##' \code{lowPass} must be specified.
##'
##' @param listOfFiles Character vector of audio file paths. Any format supported by
##'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
##' @param highPass Numeric. High-pass cutoff frequency in Hz. \code{NULL} disables
##'   high-pass filtering. Default \code{NULL}.
##' @param lowPass Numeric. Low-pass cutoff frequency in Hz. \code{NULL} disables
##'   low-pass filtering. Default \code{NULL}.
##' @param stopBand Numeric. FIR stop-band attenuation in dB (Kaiser window design).
##'   Default 96.0.
##' @param transition Numeric. FIR transition band width in Hz. Default 250.0.
##' @param useIIR Logical. Use Butterworth IIR filter instead of FIR. Default \code{FALSE}.
##' @param numIIRsections Integer. Number of 2nd-order IIR sections (filter order).
##'   Default \code{4L}.
##' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
##' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
##' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
##'   count written (invisibly). If \code{FALSE}, return an \code{AsspDataObj}.
##'   Default \code{TRUE}.
##' @param explicitExt Character. Output file extension. Default \code{"flt"}.
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
##'   preserved from libassp output (the filtered waveform, same label as the
##'   input channel), containing INT16 or REAL32 sample values.
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @details
##' Filter mode is determined by \code{highPass}/\code{lowPass}: supply only
##' \code{highPass} for a high-pass filter, only \code{lowPass} for a low-pass
##' filter, or both for a band-pass filter. \code{stopBand} and \code{transition}
##' govern the FIR Kaiser-window design; set \code{useIIR = TRUE} to switch to
##' a Butterworth IIR design instead.
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
##' # High-pass filter above 100 Hz
##' res <- trk_affilter(path2wav, highPass = 100, toFile = FALSE)
trk_affilter <- function(listOfFiles,
                         highPass = NULL,
                         lowPass = NULL,
                         stopBand = 96.0,
                         transition = 250.0,
                         useIIR = FALSE,
                         numIIRsections = 4L,
                         beginTime = 0,
                         endTime = 0,
                         toFile = TRUE,
                         explicitExt = "flt",
                         outputDirectory = NULL,
                         assertLossless = NULL,
                         logToFile = FALSE,
                         keepConverted = FALSE,
                         convertOverwrites = FALSE,
                         verbose = TRUE) {

  if (is.null(highPass) && is.null(lowPass)) {
    cli::cli_abort("At least one of {.par highPass} or {.par lowPass} must be specified.")
  }

  explicitExt     <- if (is.null(explicitExt)) "flt" else explicitExt
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

  # Build param list, skipping NULL cutoff values so C option parser
  # does not receive R NULL where it expects REAL
  filter_params <- list(
    stopBand       = as.double(stopBand),
    transition     = as.double(transition),
    useIIR         = useIIR,
    numIIRsections = as.integer(numIIRsections),
    explicitExt    = explicitExt,
    outputDirectory = outputDirectory
  )
  if (!is.null(highPass)) filter_params$highPass <- as.double(highPass)
  if (!is.null(lowPass))  filter_params$lowPass  <- as.double(lowPass)

  result <- do.call(
    processMediaFiles_LoadAndProcess,
    c(list(listOfFiles     = listOfFiles,
           beginTime       = beginTime,
           endTime         = endTime,
           nativeFiletypes = nativeFiletypes,
           fname           = "affilter",
           toFile          = toFile,
           verbose         = verbose),
      filter_params)
  )

  externalRes <- result$externalRes
  toClear     <- character(0)

  if (n_files == 1) externalRes <- externalRes[[1]]

  cleanupConvertedInputMediaFiles(toClear, keepConverted, verbose)

  return(externalRes)
}

attr(trk_affilter, "ext")             <- "flt"
attr(trk_affilter, "tracks")          <- character(0)
attr(trk_affilter, "outputType")      <- "SSFF"
attr(trk_affilter, "nativeFiletypes") <- c("wav", "au", "kay", "nist", "nsp")
attr(trk_affilter, "suggestCaching")  <- FALSE

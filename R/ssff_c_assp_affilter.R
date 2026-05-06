##' Digital filtering of audio signals
##'
##' @description Applies a FIR or IIR digital filter to audio signals listed
##' in `listOfFiles` using algorithms implemented in *libassp*
##' \insertCite{s5h}{superassp}. Supports high-pass, low-pass, and band-pass
##' configurations. At least one of `highPass` or `lowPass` must be set.
##' Input signals not in a natively supported format are converted before
##' processing; the conversion process will display warnings about input files
##' that are not in known losslessly encoded formats.
##'
##' The results are written to an SSFF formatted file with the base name of
##' the input file and the extension specified by `explicitExt`
##' (default: `"flt"`).
##'
##' @details Filter mode is determined by the combination of `highPass` and
##' `lowPass`:
##' - `highPass` only → high-pass filter
##' - `lowPass` only → low-pass filter
##' - both → band-pass filter
##'
##' `stopBand` and `transition` control the FIR filter design (Kaiser window).
##' Set `useIIR = TRUE` to use a Butterworth IIR filter instead, in which
##' case `numIIRsections` controls the filter order (number of 2nd-order
##' sections).
##'
##' Native file types: WAV (`pcm_s16le`), Sun AU, NIST, CSL (kay/nsp).
##' Conversion via [libavcodec](https://ffmpeg.org/libavcodec.html) /
##' [av::av_audio_convert].
##'
##' @param listOfFiles vector of file paths to be processed
##' @param highPass high-pass cutoff frequency in Hz; `NULL` disables
##'   high-pass (default: `NULL`)
##' @param lowPass low-pass cutoff frequency in Hz; `NULL` disables
##'   low-pass (default: `NULL`)
##' @param stopBand FIR stop-band attenuation in dB (default: 96.0)
##' @param transition FIR transition band width in Hz (default: 250.0)
##' @param useIIR use Butterworth IIR filter instead of FIR (default: `FALSE`)
##' @param numIIRsections number of 2nd-order IIR sections (default: 4L)
##' @param beginTime start of processed interval in seconds (0 = file start)
##' @param endTime end of processed interval in seconds (0 = file end)
##' @param toFile write results to file (`TRUE`) or return `AsspDataObj`
##'   (`FALSE`)
##' @param explicitExt output file extension (default: `"flt"`)
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

##' Estimate formant frequencies and bandwidths (Rcpp-optimized)
##'
##' Formant estimation of the signal(s) in `listOfFiles`. Raw
##' resonance frequency and bandwidth values are obtained by root-solving of the
##' Linear Prediction polynomial from the autocorrelation method and the
##' Split-Levinson-Algorithm (SLA). Resonances are then classified as formants
##' using the so-called Pisarenko frequencies (by-product of the SLA) and a
##' formant frequency range table derived from the nominal \ifelse{html}{\out{F<sub>1</sub>}}{\eqn{F_1}}. The
##' latter may have to be increased by about 12% for female voices (see
##' `nominalF1` and `gender` parameters). This function uses the *libassp* C library
##'  \insertCite{s5h}{superassp} for the DSP work.
##'
##' Input signals not in a natively supported file format will be converted
##' before the autocorrelation functions are computed. The conversion process
##' will display warnings about input files that are not in known losslessly
##' encoded formats.
##'
##' Default output is in SSFF binary format, with tracks containing the
##' estimated mid formant frequency of each formant (track *F[Hz]*, one column per
##' formant) and the associated formant bandwidth  (track *B[Hz]*, one column per
##' formant). If `toFile` is `TRUE`, the results will be written to a file with the
##' same name as the input file, but with an extension *.fms*.
##'
##' @details The function is a re-write of the [wrassp::forest] function, but
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
##'
##' @inheritParams fo
##' @param effectiveLength make window size effective rather than exact
##' @param nominalF1 = The nominal (assumed)  \ifelse{html}{\out{F<sub>1</sub>}}{\eqn{F_1}} frequency (default: 500.0 Hz)
##' @param gender = Use gender specific parameters? Permitted codes are  "f"[emale], "m"[ale] or "u"[nknown]. When "f", the effective window length is set to 12.5 ms and the nominal \ifelse{html}{\out{F<sub>1</sub>}}{\eqn{F_1}} to 560 Hz.
##' @param estimate insert rough frequency estimates of missing formants? By default, the frequency is set to zero.
##' @param order decrease default LPC filter order by 2 (one resonance less)
##' @param incrOrder increase default LPC filter order by 2 (one resonance more)
##' @param numFormants = The number of formants to identify. Defaults to 4, and the maximum value is 8 or half the LPC filter order)
##' @param window = <type>: set analysis window function to <type> (default: BLACKMAN)
##' @param preemphasis = <val>: set pre-emphasis factor to <val> (-1 <= val <= 0)
##' (default: dependent on sample rate and nominal \ifelse{html}{\out{F<sub>1</sub>}}{\eqn{F_1}})
##' @param windowSize window size in milliseconds
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##'
##' @return If `toFile` is `FALSE`, the function returns a list of [AsspDataObj]
##'   objects. If `toFile` is `TRUE`, the number (integer) of successfully
##'   processed and stored output files is returned.
##'
##' @seealso [wrassp::acfana]
##' @seealso [superassp::AsspWindowTypes]
##' @seealso [av::av_audio_convert]
##'
##' @useDynLib superassp, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
##'
##'
##' # calculate formant values
##' res <- trk_forest(path2wav, toFile=FALSE)
##'
##' # plot formant values
##' matplot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'           attr(res, 'startTime'),
##'         res[["F[Hz]"]],
##'         type='l',
##'         xlab='time (s)',
##'         ylab='Formant frequency (Hz)')
##'
##' @export
trk_forest <- function(listOfFiles,
                   beginTime = 0.0,
                   endTime = 0.0,
                   windowShift = 5.0,
                   windowSize = 20.0,
                   effectiveLength = TRUE,
                   nominalF1 = 500,
                   gender = 'm',
                   estimate = FALSE,
                   order = 0,
                   incrOrder = 0,
                   numFormants = 4,
                   window = 'BLACKMAN',
                   preemphasis = -0.8,
                   toFile = TRUE,
                   explicitExt = "fms",
                   outputDirectory = NULL,
                   assertLossless = NULL,
                   logToFile = FALSE,
                   convertOverwrites=FALSE,
                   keepConverted=FALSE,
                   verbose = TRUE) {

  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"fms",explicitExt)
  newTracknames <- c("F[Hz]","B[Hz]")
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
    fname = "trk_forest",
    toFile = toFile,
    verbose = verbose,
    windowShift = windowShift,
    windowSize = windowSize,
    effectiveLength = effectiveLength,
    nominalF1 = nominalF1,
    gender = gender,
    estimate = estimate,
    order = as.integer(order),
    incrOrder = as.integer(incrOrder),
    numFormants = as.integer(numFormants),
    window = window,
    preemphasis = preemphasis,
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

attr(trk_forest,"ext") <-  "fms"
attr(trk_forest,"tracks") <-  c("Fi[Hz]","Bi[Hz]")
attr(trk_forest,"outputType") <-  "SSFF"
attr(trk_forest,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(trk_forest,"suggestCaching") <-  FALSE

##' Finds the f0 using the K.Schaefer-Vincent periodicity detection algorithm (Rcpp-optimized)
##'
##' Applies Schäefer-Vincent periodicity analysis \insertCite{Schäfer-Vincent.1983.10.1159/000261691}{superassp}
##' to find \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}} (in Hz)
##' along signals listed in `listOfFiles`. Input signals not in a file format natively
##' supported will be converted before the autocorrelation functions are
##' computed. The conversion process will display warnings about input files
##' that are not in known losslessly encoded formats.
##'
##' The results will be will be written to an SSFF formated file with the base
##' name of the input file and extension *.fo* in a track *fo[Hz]*.
##'
##' @details The function is a re-write of the [wrassp::ksvF0] function, but
##' with media pre-conversion, better checking of preconditions such as the
##' input file existence, structured logging, and the use of a more modern
##' framework for user feedback. This version includes Rcpp optimizations
##' for improved performance on large batches of files.
##'
##' Optionally, location and type of the signal extrema on
##' which the \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}} data are based, may be stored in a label
##' file. The name of this file will consist of the basename of the `.fo` file and the extension '.prd'.
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
##' is *very* long, then caching of results may be warranted.
##'
##' @inheritParams acfana
##' @param gender = <code>  set gender-specific \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}} ranges; <code> may be:
##' "f[emale]" (80.0 - 640.0 Hz)
##' "m[ale]" (50.0 - 400.0 Hz)
##' "u[nknown]" (default; 50.0 - 600.0 Hz)
##' @param maxF = <freq>: set maximum \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}} value to <freq> Hz (default: 500.0)
##' @param minF = <freq>: set minimum \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}} value to <freq> Hz (default: 50.0)
##' @param minAmp = <amp>: set amplitude threshold for voiced samples to <amp> (default: 100)
##' @param maxZCR maximum zero crossing rate in Hz (for voicing detection)
##'
##' @return The number of successfully written files (if `toFile=TRUE`), or a vector of `AsspDataObj` objects (if `toFile=FALSE`).
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##'
##' @aliases foana fo_ksv
##'
##' @references
##'   \insertAllCited{}
##'
##' @aliases foana fo_ksv fo ksvfo
##'
##' @seealso \code{\link{pitch}} for a tracker of pitch
##' @useDynLib superassp, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @examples
##' # get path to audio file
##'path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
##'
##'# calculate fundamental frequency contour
##'res <- ksvfo(path2wav, toFile=FALSE)
##'
##'# plot the fundamental frequency contour
##'plot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'       attr(res, 'startTime'),
##'     res[["fo[Hz]"]],
##'     type='l',
##'     xlab='time (s)',
##'     ylab=expression(paste(f[o]," frequency (Hz)")))
##'
##' @export
fo <- ksvfo <- foana <- fo_ksv <- function(listOfFiles = NULL,
                                           beginTime = 0.0,
                                           endTime = 0.0,
                                           windowShift = 5.0,
                                           gender = 'u',
                                           maxF = 600,
                                           minF = 50,
                                           minAmp = 50,
                                           maxZCR = 3000.0,
                                           toFile = TRUE,
                                           explicitExt = "fo",
                                           outputDirectory = NULL,
                                           assertLossless = NULL,
                                           logToFile = FALSE,
                                           convertOverwrites=FALSE,
                                           keepConverted=FALSE,
                                           verbose = TRUE) {

  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"fo",explicitExt)
  newTracknames <- c("fo[Hz]")
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
    fname = "f0ana",
    toFile = toFile,
    verbose = verbose,
    windowShift = windowShift,
    gender = gender,
    maxF = maxF,
    minF = minF,
    minAmp = minAmp,
    maxZCR = maxZCR,
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

attr(ksvfo,"ext") <-  "fo"
attr(ksvfo,"tracks") <-  c("fo[Hz]")
attr(ksvfo,"outputType") <-  "SSFF"
attr(ksvfo,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(ksvfo,"suggestCaching") <-  FALSE

attr(foana,"ext") <-  "fo"
attr(foana,"tracks") <-  c("fo[Hz]")
attr(foana,"outputType") <-  "SSFF"
attr(foana,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(foana,"suggestCaching") <-  FALSE

attr(fo_ksv,"ext") <-  "fo"
attr(fo_ksv,"tracks") <-  c("fo[Hz]")
attr(fo_ksv,"outputType") <-  "SSFF"
attr(fo_ksv,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(fo_ksv,"suggestCaching") <-  FALSE

attr(fo,"ext") <-  "fo"
attr(fo,"tracks") <-  c("fo[Hz]")
attr(fo,"outputType") <-  "SSFF"
attr(fo,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(fo,"suggestCaching") <-  FALSE

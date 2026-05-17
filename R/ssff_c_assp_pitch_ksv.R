##' Track fundamental frequency using the KSV periodicity detector
##'
##' Estimates the fundamental frequency \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}}
##' using the Schäfer-Vincent periodicity detection algorithm
##' \insertCite{Schäfer-Vincent.1983.10.1159/000261691}{superassp} implemented in
##' the *libassp* C library \insertCite{s5h}{superassp}. This extremum-based
##' method is fast and works directly on the waveform without spectral analysis.
##'
##' @param listOfFiles Character vector of audio file paths. Any format supported by
##'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
##' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
##' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
##' @param windowShift Numeric. Frame shift in milliseconds; sets output frame rate
##'   (1000 / windowShift Hz). Default 5 ms.
##' @param gender Character. Gender-specific \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}}
##'   search range: \code{"f"} (female, 80–640 Hz), \code{"m"} (male, 50–400 Hz),
##'   \code{"u"} (unknown, 50–600 Hz, default).
##' @param maxF Numeric. Maximum \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}} in Hz.
##'   Default 600.
##' @param minF Numeric. Minimum \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}} in Hz.
##'   Default 50.
##' @param minAmp Numeric. Minimum waveform amplitude threshold for voiced frames.
##'   Default 50.
##' @param maxZCR Numeric. Maximum zero-crossing rate in Hz for voicing detection.
##'   Default 3000.0.
##' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
##'   count written (invisibly). If \code{FALSE}, return an \code{AsspDataObj}.
##'   Default \code{TRUE}.
##' @param explicitExt Character. Output file extension. Default \code{"fo"}.
##' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
##'   writes alongside the input file.
##' @param assertLossless Character vector of additional file extensions to treat as
##'   losslessly encoded.
##' @param logToFile Logical. Write processing log to a file in \code{outputDirectory}
##'   rather than the console. Default \code{FALSE}.
##' @param convertOverwrites Logical. Allow transcoding to overwrite existing files.
##'   Default \code{FALSE}.
##' @param keepConverted Logical. Retain intermediate transcoded files. Default \code{FALSE}.
##' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
##'   \describe{
##'     \item{\code{fo[Hz]}}{REAL32, Hz, n_frames x 1 column.
##'       Estimated fundamental frequency; 0 indicates unvoiced frames.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 200 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @details
##' \code{gender} sets the default \code{minF}/\code{maxF} search range but is
##' overridden by explicit \code{minF}/\code{maxF} values. \code{minAmp} and
##' \code{maxZCR} control voicing detection independently of the pitch range.
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
##' @seealso [wrassp::ksvF0]
##' @useDynLib superassp, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @export
##'
##' @examples
##' # get path to audio file
##'path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
##'
##'# calculate fundamental frequency contour
##'res <- trk_pitch_ksv(path2wav, toFile=FALSE)
##'
##'# plot the fundamental frequency contour
##'plot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'       attr(res, 'startTime'),
##'     res[["fo[Hz]"]],
##'     type='l',
##'     xlab='time (s)',
##'     ylab=expression(paste(f[o]," frequency (Hz)")))
##'
trk_pitch_ksv <- function(listOfFiles = NULL,
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
    format_apply_msg(funName, n_files, beginTime, endTime)
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

attr(trk_pitch_ksv,"ext") <-  "fo"
attr(trk_pitch_ksv,"tracks") <-  c("fo[Hz]")
attr(trk_pitch_ksv,"outputType") <-  "SSFF"
attr(trk_pitch_ksv,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(trk_pitch_ksv,"suggestCaching") <-  FALSE

##' @rdname trk_pitch_ksv
##' @usage NULL
##' @export
trk_ksvfo <- function(...) {
  lifecycle::deprecate_warn("2.5.3", "trk_ksvfo()", "trk_pitch_ksv()")
  trk_pitch_ksv(...)
}

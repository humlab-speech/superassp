##' Track short-term autocorrelation function
##'
##' Computes the short-term autocorrelation function (ACF) of audio signals
##' using the *libassp* C library \insertCite{s5h}{superassp}. Useful as a
##' front-end feature for voicing detection and LP-based analysis. Prefer
##' \code{trk_acf} over manual lag computation when frame-synchronous output
##' in SSFF format is needed.
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
##' @param effectiveLength Logical. Make window size effective rather than exact. Default \code{TRUE}.
##' @param window Character. Analysis window function type. Default \code{"BLACKMAN"}.
##'   See [superassp::AsspWindowTypes] for supported types.
##' @param analysisOrder Integer. Number of lag coefficients per frame. \code{0} sets
##'   order to sample rate in kHz + 3 (e.g. 19 for 16 kHz audio). Default 0.
##' @param energyNormalization Logical. Compute energy-normalised ACF. Default \code{FALSE}.
##' @param lengthNormalization Logical. Compute length-normalised ACF. Default \code{FALSE}.
##' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
##'   count written (invisibly). If \code{FALSE}, return an \code{AsspDataObj}.
##'   Default \code{TRUE}.
##' @param explicitExt Character. Output file extension. Default \code{"acf"}.
##' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
##'   writes alongside the input file.
##' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
##' @param assertLossless Character vector of additional file extensions to treat as
##'   losslessly encoded.
##' @param logToFile Logical. Write processing log to a file in \code{outputDirectory}
##'   rather than the console. Default \code{FALSE}.
##' @param keepConverted Logical. Retain intermediate transcoded files. Default \code{FALSE}.
##' @param convertOverwrites Logical. Allow transcoding to overwrite existing files.
##'   Default \code{FALSE}.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
##'   \describe{
##'     \item{\code{ACF}}{REAL32, dimensionless, n_frames x \code{analysisOrder} columns.
##'       Autocorrelation coefficients at lags 0 … analysisOrder-1.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 200 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @details
##' \code{analysisOrder = 0} selects an order equal to the sample rate in kHz + 3.
##' Energy normalisation divides each frame's ACF by its lag-0 value.
##' Length normalisation divides by frame length. Both can be combined.
##'
##' @seealso [wrassp::acfana]
##' @seealso [superassp::AsspWindowTypes]
##' @seealso [av::av_audio_convert]
##'
##' @useDynLib superassp, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @examples
##' # get path to audio file
#' path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
#'
#' # calculate short-term autocorrelation
#' res <- trk_acf(path2wav, toFile=FALSE)
#'
#' # plot short-term autocorrelation values
#' matplot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
#'         attr(res, 'startTime'),
#'         res$acf,
#'         type='l',
#'         xlab='time (s)',
#'         ylab='Short-term autocorrelation values')
##'
##' @export
##' @references
##'   \insertAllCited{}
##'

trk_acf <- function(listOfFiles,
                   beginTime = 0,
                   centerTime = FALSE,
                   endTime = 0,
                   windowShift = 5,
                   windowSize = 20,
                   effectiveLength = TRUE,
                   window = "BLACKMAN",
                   analysisOrder = 0,
                   energyNormalization = FALSE,
                   lengthNormalization = FALSE,
                   toFile = TRUE,
                   explicitExt = "acf",
                   outputDirectory = NULL,
                   assertLossless = NULL,
                   logToFile = FALSE,
                   keepConverted=FALSE,
                   convertOverwrites=FALSE,
                   verbose = TRUE) {



  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"acf",explicitExt)
  newTracknames <- "ACF"  ## Only used if SSFF tracks needs to be renamed from the called function (in C) before returning the SSFF track obj
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

  knownLossless <- c(assertLossless,knownLossless()) #Use the user asserted information about lossless encoding, in addition to what is already known by superassp

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
    fname = "acfana",
    toFile = toFile,
    verbose = verbose,
    centerTime = centerTime,
    windowShift = windowShift,
    windowSize = windowSize,
    effectiveLength = effectiveLength,
    window = window,
    analysisOrder = as.integer(analysisOrder),
    energyNormalization = energyNormalization,
    lengthNormalization = lengthNormalization,
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

# Function attributes
attr(trk_acf,"ext") <-  "acf"
attr(trk_acf,"tracks") <-  c("ACF")
attr(trk_acf,"outputType") <-  "SSFF"
attr(trk_acf,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(trk_acf,"suggestCaching") <-  FALSE

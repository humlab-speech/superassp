##' Track formant frequencies and bandwidths (FOREST)
##'
##' Estimates vocal tract resonance (formant) frequencies and their bandwidths
##' using the FOREST algorithm from the *libassp* C library
##' \insertCite{s5h}{superassp}. Root-solving of the LP polynomial is guided by
##' Pisarenko frequencies from the Split-Levinson Algorithm (SLA) to classify
##' resonances as formants. Prefer \code{trk_formant_forest} when broad
##' compatibility and speed on large corpora are priorities.
##'
##' @param listOfFiles Character vector of audio file paths. Any format supported by
##'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
##' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
##' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
##' @param windowShift Numeric. Frame shift in milliseconds; sets output frame rate
##'   (1000 / windowShift Hz). Default 5 ms.
##' @param windowSize Numeric. Analysis window size in milliseconds. Default 20 ms.
##' @param effectiveLength Logical. Make window size effective rather than exact.
##'   Default \code{TRUE}.
##' @param nominalF1 Numeric. Assumed F1 frequency in Hz used to build the formant
##'   range table. Increase by ~12% for female voices. Default 500.0 Hz.
##' @param gender Character. Gender-specific parameter preset: \code{"f"} (female,
##'   sets window to 12.5 ms and nominalF1 to 560 Hz), \code{"m"} (male), or
##'   \code{"u"} (unknown, default).
##' @param estimate Logical. Insert rough frequency estimates for missing formants
##'   rather than returning zero. Default \code{FALSE}.
##' @param order Integer. Decrease the default LPC filter order by 2 (one fewer
##'   resonance). Default 0 (no change).
##' @param incrOrder Integer. Increase the default LPC filter order by 2 (one more
##'   resonance). Default 0 (no change).
##' @param numFormants Integer. Number of formants to track (maximum 8 or half the
##'   LPC order). Default 4.
##' @param window Character. Analysis window function type. Default \code{"BLACKMAN"}.
##'   See [superassp::AsspWindowTypes].
##' @param preemphasis Numeric. Pre-emphasis factor (-1 <= val <= 0); default is
##'   sample-rate- and nominalF1-dependent.
##' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
##'   count written (invisibly). If \code{FALSE}, return an \code{AsspDataObj}.
##'   Default \code{TRUE}.
##' @param explicitExt Character. Output file extension. Default \code{"fms"}.
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
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
##'   \describe{
##'     \item{\code{F[Hz]}}{REAL32, Hz, n_frames x \code{numFormants} columns.
##'       Estimated centre frequency of each formant (0 = missing).}
##'     \item{\code{B[Hz]}}{REAL32, Hz, n_frames x \code{numFormants} columns.
##'       Bandwidth of each formant.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 200 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @details
##' The \code{gender} preset overrides \code{windowSize} and \code{nominalF1}.
##' Set \code{estimate = TRUE} to fill missing formant slots with rough estimates
##' rather than zeros, which can help downstream processing.
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##'
##' @seealso [wrassp::forest]
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
##' res <- trk_formant_forest(path2wav, toFile=FALSE)
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
trk_formant_forest <- function(listOfFiles,
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
    format_apply_msg(funName, n_files, beginTime, endTime)
  }

  # Use unified load-and-process helper (works for all file formats)
  result <- processMediaFiles_LoadAndProcess(
    listOfFiles = listOfFiles,
    beginTime = beginTime,
    endTime = endTime,
    nativeFiletypes = nativeFiletypes,
    fname = "forest",
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

  # Set attributes on returned object(s) for proper template expansion
  if(!toFile) {
    if(n_files == 1) {
      attr(externalRes, "func") <- "trk_formant_forest"
      attr(externalRes, "tracks") <- c("Fi[Hz]", "Bi[Hz]")
    } else {
      for(i in seq_along(externalRes)) {
        attr(externalRes[[i]], "func") <- "trk_formant_forest"
        attr(externalRes[[i]], "tracks") <- c("Fi[Hz]", "Bi[Hz]")
      }
    }
  }

  #### Cleanup ####
  cleanupConvertedInputMediaFiles(toClear, keepConverted, verbose)

  return(externalRes)
}

attr(trk_formant_forest,"ext") <-  "fms"
attr(trk_formant_forest,"tracks") <-  c("F[Hz]","B[Hz]")
attr(trk_formant_forest,"outputType") <-  "SSFF"
attr(trk_formant_forest,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(trk_formant_forest,"suggestCaching") <-  FALSE

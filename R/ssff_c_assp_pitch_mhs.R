##' Track pitch using the Modified Harmonic Sieve algorithm
##'
##' Estimates pitch in Hz using Michel Scheffers' Modified Harmonic Sieve (MHS)
##' algorithm implemented in the *libassp* C library \insertCite{s5h}{superassp}.
##' MHS operates in the frequency domain and is robust to noise; it is a
##' complementary alternative to the waveform-based \code{trk_ksvfo}.
##'
##' @param listOfFiles Character vector of audio file paths. Any format supported by
##'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
##' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
##' @param centerTime Numeric or logical. Single-frame analysis time point in seconds;
##'   overrides \code{beginTime}, \code{endTime}, and \code{windowShift}. Default \code{FALSE}.
##' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
##' @param windowShift Numeric. Frame shift in milliseconds; sets output frame rate
##'   (1000 / windowShift Hz). Default 5 ms.
##' @param gender Character. Gender-specific pitch search range: \code{"f"} (female),
##'   \code{"m"} (male), \code{"u"} (unknown, default).
##' @param maxF Numeric. Maximum pitch in Hz. Default 600.0.
##' @param minF Numeric. Minimum pitch in Hz. Default 50.0.
##' @param minAmp Numeric. Minimum signal amplitude threshold. Default 50.0.
##' @param minAC1 Numeric. Minimum first autocorrelation coefficient. Default 0.25.
##' @param minRMS Numeric. Minimum RMS amplitude in dB for voiced detection. Default 18.0.
##' @param maxZCR Numeric. Maximum zero-crossing rate in Hz for voiced detection.
##'   Default 3000.0.
##' @param minProb Numeric. Minimum harmonic sieve fit quality (0–1) for a frame to be
##'   considered voiced. Default 0.52.
##' @param plainSpectrum Logical. Use plain (non-pre-emphasised) spectrum. Default \code{FALSE}.
##' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
##'   count written (invisibly). If \code{FALSE}, return an \code{AsspDataObj}.
##'   Default \code{TRUE}.
##' @param explicitExt Character. Output file extension. Default \code{"pit"}.
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
##'     \item{\code{pitch[Hz]}}{REAL32, Hz, n_frames x 1 column.
##'       Estimated pitch frequency; 0 indicates unvoiced frames.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 200 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @details
##' Voicing is determined by joint thresholds on \code{minAmp}, \code{minAC1},
##' \code{minRMS}, \code{maxZCR}, and \code{minProb}. Increase \code{minProb}
##' to reduce false voiced decisions in noisy conditions.
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##'
##' @aliases mhspitch pitch_mhs pitch
##'
##' @seealso [wrassp::mhsF0]
##' @seealso [superassp::trk_ksvfo]
##'
##' @export
##' @useDynLib superassp, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @references
##'   \insertAllCited{}
##'
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
##'
##' # calculate short-term autocorrelation
##' res <- pitch(path2wav, toFile=FALSE)
##'
##' # plot fundamental frequency contour
##' plot(seq(0,numRecs.AsspDataObj(res) - 1) / rate.AsspDataObj(res) +
##'        attr(res, 'startTime'),
##'      res[["pitch[Hz]"]],
##'      type='l',
##'      xlab='time (s)',
##'      ylab="Pitch (Hz)")
##'
trk_pitch_mhs <- function(listOfFiles = NULL,
                                             beginTime = 0.0,
                                             centerTime = FALSE,
                                             endTime = 0.0,
                                             windowShift = 5.0,
                                             gender = 'u',
                                             maxF = 600.0,
                                             minF = 50.0,
                                             minAmp = 50.0,
                                             minAC1 = 0.25,
                                             minRMS = 18.0,
                                             maxZCR = 3000.0,
                                             minProb = 0.52,
                                             plainSpectrum = FALSE,
                                             toFile = TRUE,
                                             explicitExt = "pit",
                                             outputDirectory = NULL,
                                             assertLossless = NULL,
                                             logToFile = FALSE,
                                             convertOverwrites=FALSE,
                                             keepConverted=FALSE,
                                             verbose = TRUE){

  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"pit",explicitExt)
  newTracknames <- c("pitch[Hz]")
  nativeFiletypes <- c("wav","au","kay","nist","nsp")

  if(! gender %in% c('u','f','m')) cli::cli_abort(c("Incorrect specification of gender: {.val {gender}}",
                                                  "i"="Valid  specifications are {.val {c('u','f','m')}} for 'unspecified', and cis 'male', and 'female' genders, respectively."))

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
    fname = "mhspitch",
    toFile = toFile,
    verbose = verbose,
    centerTime = centerTime,
    windowShift = windowShift,
    gender = gender,
    maxF = maxF,
    minF = minF,
    minAmp = minAmp,
    minAC1 = minAC1,
    minRMS = minRMS,
    maxZCR = maxZCR,
    minProb = minProb,
    plainSpectrum = plainSpectrum,
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

attr(trk_pitch_mhs,"ext") <-  "pit"
attr(trk_pitch_mhs,"tracks") <-  c("pitch[Hz]")
attr(trk_pitch_mhs,"outputType") <-  "SSFF"
attr(trk_pitch_mhs,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(trk_pitch_mhs,"suggestCaching") <-  FALSE

pitch_mhs <- pitch <- trk_pitch_mhs

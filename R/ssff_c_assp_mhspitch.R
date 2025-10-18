##' Pitch analysis using a Modified Harmonic Sieve (Rcpp-optimized)
##'
##' This function finds the pitch (in Hz) along files in <listOfFile> using
##' Michel Scheffers' Modified Harmonic Sieve algorithm implmented in *libassp* \insertCite{s5h}{superassp}.
##' Input signals not in a file format natively
##' supported will be converted before the autocorrelation functions are
##' computed. The conversion process will display warnings about input files
##' that are not in known losslessly encoded formats.
##'
##' The results will be will be written to an SSFF formated file with the base
##' name of the input file and extension *.pit* in a track *pitch[Hz]*.
##'
##' @details The function is a re-write of the [wrassp::mhsF0] function, but
##' with media pre-conversion, better checking of preconditions such as the
##' input file existence, structured logging, and the use of a more modern
##' framework for user feedback. This version includes Rcpp optimizations
##' for improved performance on large batches of files.
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
##' is *very* large, then caching of results may be warranted.
##'
##' @inheritParams fo
##' @param centerTime sets single frame analysis time point (in seconds)
##' @param minAmp = <amp>:  minimum signal amplitude (default: 50)
##' @param minAC1 = <freq>: minimum 1st correlation coefficient (default: 0.250)
##' @param minRMS = <num>:  minimum RMS amplitude in dB (default: 18.0)
##' @param maxZCR = <freq>: maximum zero crossing rate in Hz (default: 3000)
##' @param minProb = <num>: minimum quality value of \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}} fit (default: 0.520)
##' @param plainSpectrum plain spectrum
##'
##' @return The number of successfully written files (if `toFile=TRUE`), or a vector of `AsspDataObj` objects (if `toFile=FALSE`).
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##'
##' @aliases mhspitch pitch_mhs pitch
##'
##' @seealso \code{\link{ksv_fo}} for an algorithm for tracking the fundamental frequency \ifelse{html}{\out{f<sub>o</sub>}}{\eqn{f_o}}.
##'
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
##' @export
pitch <- mhspitch <- pitch_mhs <-function(listOfFiles = NULL,
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
    cli::cli_inform("Applying {.fun {funName}} to {cli::no(n_files)} recording{?s}")
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

attr(mhspitch,"ext") <-  "pit"
attr(mhspitch,"tracks") <-  c("pitch[Hz]")
attr(mhspitch,"outputType") <-  "SSFF"
attr(mhspitch,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(mhspitch,"suggestCaching") <-  FALSE

attr(pitch_mhs,"ext") <-  "pit"
attr(pitch_mhs,"tracks") <-  c("pitch[Hz]")
attr(pitch_mhs,"outputType") <-  "SSFF"
attr(pitch_mhs,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(pitch_mhs,"suggestCaching") <-  FALSE

attr(pitch,"ext") <-  "pit"
attr(pitch,"tracks") <-  c("pitch[Hz]")
attr(pitch,"outputType") <-  "SSFF"
attr(pitch,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(pitch,"suggestCaching") <-  FALSE

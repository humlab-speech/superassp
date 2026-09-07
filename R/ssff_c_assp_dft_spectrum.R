##' Track short-term DFT power spectrum
##'
##' Computes a short-term power spectrum via the Fast Fourier Transform using
##' the *libassp* C library \insertCite{s5h}{superassp}. Produces an
##' unsmoothed narrow-band spectrum from 0 Hz to the Nyquist rate. Prefer this
##' function when raw spectral detail is needed; use \code{trk_css_spectrum} or
##' \code{trk_lps_spectrum} for smoothed spectral envelopes.
##'
##' @inheritParams trk_lps_spectrum
##' @param bandwidth Numeric. Effective analysis bandwidth in Hz. Default 0
##'   yields the minimum bandwidth determined by the FFT length.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
##'   \describe{
##'     \item{\code{DFT[dB]}}{REAL32, dB power, n_frames x (FFT_length/2 + 1) columns.
##'       Power spectral amplitude from 0 Hz to the Nyquist rate.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 200 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @details
##' The FFT length is determined by \code{resolution} unless overridden by
##' \code{fftLength}. \code{bandwidth} widens the effective analysis window,
##' trading spectral resolution for reduced side-lobe leakage.
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##'
##' @seealso [wrassp::dftSpectrum]
##' @seealso [superassp::AsspWindowTypes]
##' @seealso [av::av_audio_convert]
##'
##' @useDynLib superassp, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("samples", "sustained", package = "superassp"),
##'                        pattern = glob2rx("a1.wav"),
##'                        full.names = TRUE)
##'
##' # calculate dft spectrum
##' res <- trk_dft_spectrum(path2wav, toFile=FALSE)
##'
##' # plot spectral values at midpoint of signal
##' plot(res[["DFT[dB]"]][dim(res[["DFT[dB]"]])[1]/2,],
##'      type='l',
##'      xlab='spectral value index',
##'      ylab='spectral value')
##'
##' @usage trk_dft_spectrum(listOfFiles, beginTime = 0, centerTime = FALSE, endTime = 0, resolution = 40, fftLength = 0, windowShift = 5, window = "BLACKMAN", bandwidth = 0, toFile = FALSE, explicitExt = "dft", outputDirectory = NULL, assertLossless = NULL, logToFile = FALSE, keepConverted = FALSE, convertOverwrites = FALSE, verbose = TRUE)
##' @param beginTime Start time for the extracted portion in seconds. Default: NULL (beginning of signal). Note: uses `beginTime`/`endTime` (seconds) matching DSP function conventions, unlike [read_audio()] which uses `begin`/`end`.
##' @param centerTime Numeric or logical. Single-frame analysis time point in seconds; overrides \code{beginTime}, \code{endTime}, and \code{windowShift}. Default \code{FALSE}.
##' @param endTime The end time of the section of the sound files that should be analysed (in seconds). Use 0 for end of file.
##' @param resolution Numeric. Target FFT frequency resolution in Hz; the FFT length is set to the smallest power-of-2 meeting this target. Default 40.0.
##' @param fftLength Integer. Explicit FFT length in points; overrides \code{resolution}. Default 0 (use \code{resolution}).
##' @param windowShift Numeric. Frame shift in milliseconds; sets output frame rate (\code{1000 / windowShift} Hz). Default 5.0 ms (200 Hz). Must be strictly less than 32 ms (the 512-sample analysis window at 16 kHz). Values other than the training default (5 ms) may slightly reduce accuracy.
##' @param window Character. Analysis window function type. Default \code{"BLACKMAN"}. See [superassp::AsspWindowTypes] for supported types.
##' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the count written. If \code{FALSE}, return an \code{AsspDataObj} (single file only). Default \code{TRUE}.
##' @param explicitExt By default, a character "d" will be prepended to the file name suffix when writing the output to file. The user can also specify an explicit extension which will be used instead.
##' @param outputDirectory The directory where the slice file should be stored. If not defiled (NULL), the sparse slice file will placed in the same folder as the media file.
##' @param assertLossless Character vector of additional file extensions to treat as losslessly encoded.
##' @param logToFile Logical. Write processing log to a file in \code{outputDirectory} rather than the console. Default \code{FALSE}.
##' @param keepConverted Logical. Retain intermediate transcoded files. Default \code{FALSE}.
##' @param convertOverwrites Logical. Allow transcoding to overwrite existing files. Default \code{FALSE}.
##' @param verbose Logical. Show a progress bar (sequential path) or a progress-aware parallel apply (`pbapply`/`pbmcapply`, if installed).
##' @export
##'
##' @references
##'   \insertAllCited{}
##'
'trk_dft_spectrum' <- function(listOfFiles,
                          beginTime = 0.0,
                          centerTime = FALSE,
                          endTime = 0.0,
                          resolution = 40.0,
                          fftLength = 0,
                          windowShift = 5.0,
                          window = 'BLACKMAN',
                          bandwidth = 0.0,
                          toFile = FALSE,
                          explicitExt = "dft",
                          outputDirectory = NULL,
                          assertLossless = NULL,
                          logToFile = FALSE,
                          keepConverted=FALSE,
                          convertOverwrites=FALSE,
                          verbose = TRUE){

  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"dft",explicitExt)
  newTracknames <- "DFT[dB]"
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
    fname = "spectrum",
    toFile = toFile,
    verbose = verbose,
    centerTime = centerTime,
    resolution = resolution,
    fftLength = as.integer(fftLength),
    windowShift = windowShift,
    window = window,
    bandwidth = bandwidth,
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
attr(trk_dft_spectrum,"ext") <-  "dft"
attr(trk_dft_spectrum,"tracks") <-  c("DFT[dB]")
attr(trk_dft_spectrum,"outputType") <-  "SSFF"
attr(trk_dft_spectrum,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(trk_dft_spectrum,"suggestCaching") <-  FALSE

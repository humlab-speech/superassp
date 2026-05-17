##' Track short-term cepstral coefficients
##'
##' Computes the short-term real cepstrum of audio signals via FFT using the
##' *libassp* C library \insertCite{s5h,Oppenheim.2004.10.1109/msp.2004.1328092,Childers.1977.10.1109/proc.1977.10747}{superassp}.
##' Output coefficients index quefrency from 0 to FFT_length/2 in steps of
##' 0.5 ms (at the default 40 Hz resolution). Useful for pitch-period
##' detection and spectral tilt estimation.
##'
##' @inheritParams trk_acf
##' @param resolution Numeric. Target FFT frequency resolution in Hz; the FFT
##'   length is set to the smallest power-of-2 meeting this target. Default 40.0.
##' @param fftLength Integer. Explicit FFT length in points; overrides
##'   \code{resolution}. Default 0 (use \code{resolution}).
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
##'   \describe{
##'     \item{\code{C[dB]}}{REAL32, dB amplitude, n_frames x (FFT_length/2 + 1) columns.
##'       Each column corresponds to a quefrency of col_index × 0.5 ms.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 200 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @details
##' The number of coefficients per frame equals FFT_length / 2 + 1 (one-sided,
##' not mirrored). Use \code{fftLength} to fix the transform size; otherwise
##' \code{resolution} governs it. Pre-emphasis is not applied by this function;
##' apply \code{trk_afdiff} first if needed.
##'
##' @author Raphael Winkelmann
##' @author Lasse Bombien
##' @author Fredrik Nylén
##'
##' @seealso [wrassp::cepstrum]
##' @seealso [superassp::AsspWindowTypes]
##' @seealso [av::av_audio_convert]
##'
##' @useDynLib superassp, .registration = TRUE
##' @importFrom Rcpp sourceCpp
##' @examples
##' # get path to audio file
##' path2wav <- list.files(system.file("samples","sustained", package = "superassp"), pattern = glob2rx("a1.wav"), full.names = TRUE)
##'
##' # calulate cepstrum
##' res <- trk_cepstrum(path2wav, toFile=FALSE)
##'
##' # plot cepstral values at midpoint of signal
##' plot(y=res[["C[dB]"]][400,],
##'     x=seq(1,ncol(res[["C[dB]"]])),
##'     type='l',
##'     xlab='Quefrency (ms)',
##'     ylab='Amplitude (dB)')
##'
##' @export
##'
##' @references
##'   \insertAllCited{}
##'
trk_cepstrum <- function(listOfFiles,
                       beginTime = 0.0,
                       centerTime = FALSE,
                       endTime = 0.0,
                       resolution = 40.0,
                       fftLength = 0,
                       windowShift = 5.0,
                       window = 'BLACKMAN',
                    toFile = TRUE,
                    explicitExt = "cep",
                    outputDirectory = NULL,
                    assertLossless = NULL,
                    logToFile = FALSE,
                    keepConverted=FALSE,
                    convertOverwrites=FALSE,
                    verbose = TRUE){

  ## Initial constants -- specific to this function
  explicitExt <- ifelse(is.null(explicitExt),"cep",explicitExt)
  newTracknames <- "C[dB]"
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
    spectrumType = 'CEP',
    resolution = resolution,
    fftLength = as.integer(fftLength),
    windowShift = windowShift,
    window = window,
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
attr(trk_cepstrum,"ext") <-  "cep"
attr(trk_cepstrum,"tracks") <-  c("C[dB]")
attr(trk_cepstrum,"outputType") <-  "SSFF"
attr(trk_cepstrum,"nativeFiletypes") <-  c("wav","au","kay","nist","nsp")
attr(trk_cepstrum,"suggestCaching") <-  FALSE

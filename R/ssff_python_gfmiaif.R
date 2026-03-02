#' Compute vocal tract and glottis source-filter separation using GFM-IAIF
#'
#' GFM-IAIF \insertCite{perrotin2019gfmiaif}{superassp} (Glottal Flow Model-based
#' Iterative Adaptive Inverse Filtering) is an algorithm for source-filter separation
#' of speech signals. It extends the classical IAIF algorithm with improved pre-emphasis,
#' allowing extraction of wide-band glottis response that incorporates both glottal
#' formant and spectral tilt characteristics.
#'
#' The algorithm decomposes speech into three components:
#' \itemize{
#'   \item \strong{Vocal tract filter}: Models formant structure and resonances (default: 48th order)
#'   \item \strong{Glottis source filter}: Models glottal excitation characteristics (3rd order, highly recommended)
#'   \item \strong{Lip radiation filter}: Models radiation at the lips (2nd order, d=0.99)
#' }
#'
#' Implemented in native C++ (no Python dependency).
#'
#' @param listOfFiles A character vector of file paths to audio files, or an AVAudio S7 object.
#' @param beginTime The start time of the section of the sound file that should be processed (in seconds).
#'   Default is 0.0 (start of file).
#' @param centerTime If TRUE, use window center times for timestamps. If FALSE, use window start times.
#'   Default is FALSE.
#' @param endTime The end time of the section of the sound file that should be processed (in seconds).
#'   Default is 0.0 (end of file).
#' @param windowShift The frame shift in milliseconds. Default is 10.0 ms.
#' @param windowSize The frame size in milliseconds. Default is 32.0 ms (~512 samples at 16kHz).
#'   Must be large enough to accommodate nv+1 samples.
#' @param nv Vocal tract LPC order. Default is 48. Controls spectral resolution of formant structure.
#' @param ng Glottis LPC order. Default is 3, and highly recommended.
#' @param d Leaky integration coefficient for lip radiation modeling. Default is 0.99.
#' @param window Window function type. Default is "HANN". Options: "HANN", "HAMMING", "BLACKMAN".
#' @param explicitExt The file extension for the output SSFF file. Default is "gfm".
#' @param outputDirectory The directory where output files should be written.
#'   If NULL (default), files are written to the same directory as the input files.
#' @param toFile If TRUE (default), write results to SSFF files and return file count.
#'   If FALSE, return AsspDataObj for single file processing.
#' @param verbose If TRUE (default), print progress messages.
#'
#' @return
#'   If \code{toFile = TRUE}: Returns the number of successfully processed files.
#'   If \code{toFile = FALSE}: Returns an AsspDataObj with tracks:
#'   \itemize{
#'     \item \strong{av_0 to av_<nv>}: Vocal tract LP coefficients (REAL64 format), nv+1 columns
#'     \item \strong{ag_0 to ag_<ng>}: Glottis LP coefficients (REAL64 format), ng+1 columns
#'     \item \strong{al_0, al_1}: Lip radiation LP coefficients (REAL64 format), 2 columns
#'   }
#'
#' @details
#' \strong{GFM-IAIF Algorithm:}
#'
#' The algorithm performs six steps:
#' \enumerate{
#'   \item Pre-frame addition: Adds mean-normalized ramp to reduce edge effects
#'   \item Lip radiation cancellation: Applies leaky integration
#'   \item Gross glottis estimation: Iterative 1st-order LPC (ng iterations)
#'   \item Gross vocal tract estimation: nv-order LPC after glottis removal
#'   \item Fine glottis estimation: ng-order LPC after vocal tract removal
#'   \item Fine vocal tract estimation: Final nv-order LPC refinement
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' result <- trk_gfmiaif("speech.wav", toFile = FALSE)
#'
#' # Lower vocal tract order for faster processing
#' result <- trk_gfmiaif("speech.wav", nv = 24, windowShift = 5.0, toFile = FALSE)
#'
#' # Batch processing to SSFF files
#' files <- c("file1.wav", "file2.wav", "file3.mp3")
#' trk_gfmiaif(files, toFile = TRUE)
#' }
#'
#' @references
#'   \insertAllCited{}
#'
#' @export
trk_gfmiaif <- function(listOfFiles,
                       beginTime = 0.0,
                       centerTime = FALSE,
                       endTime = 0.0,
                       windowShift = 10.0,
                       windowSize = 32.0,
                       nv = 48L,
                       ng = 3L,
                       d = 0.99,
                       window = "HANN",
                       explicitExt = "gfm",
                       outputDirectory = NULL,
                       toFile = TRUE,
                       verbose = TRUE) {

  # Validate inputs
  if (length(listOfFiles) > 1 && !toFile) {
    stop("length(listOfFiles) is > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.",
         call. = FALSE)
  }

  nv <- as.integer(nv)
  ng <- as.integer(ng)

  if (nv < 1 || nv > 100) stop("nv must be between 1 and 100", call. = FALSE)
  if (ng < 1 || ng > 10) stop("ng must be between 1 and 10", call. = FALSE)

  if (ng != 3 && verbose) {
    warning("GFM-IAIF is designed for ng=3. Using ng=", ng, " may not be optimal.",
            call. = FALSE)
  }

  if (d < 0.9 || d > 0.999) stop("d must be between 0.9 and 0.999", call. = FALSE)

  # Map window type
  window_type <- tolower(window)
  if (!window_type %in% c("hann", "hamming", "blackman")) {
    stop("window must be one of: HANN, HAMMING, BLACKMAN", call. = FALSE)
  }

  window_shift_sec <- windowShift / 1000.0
  window_size_sec <- windowSize / 1000.0

  # Create file/time dataframe
  tryCatch({
    fileBeginEnd <- data.frame(
      listOfFiles = listOfFiles,
      beginTime = beginTime,
      endTime = endTime,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    stop("The beginTime and endTime must either be a single value or the same length as listOfFiles",
         call. = FALSE)
  })

  # Check files exist
  filesEx <- file.exists(listOfFiles)
  if (!all(filesEx)) {
    stop("Unable to find the sound file(s): ", paste(listOfFiles[!filesEx], collapse = ", "),
         call. = FALSE)
  }

  outListOfFiles <- c()

  for (i in seq_len(nrow(fileBeginEnd))) {
    origSoundFile <- normalizePath(fileBeginEnd[i, "listOfFiles"], mustWork = TRUE)
    bt <- fileBeginEnd[i, "beginTime"]
    et <- fileBeginEnd[i, "endTime"]

    if (verbose) {
      message(sprintf("Processing file %d/%d: %s", i, nrow(fileBeginEnd), basename(origSoundFile)))
    }

    # Load audio via av
    audio_data <- av::read_audio_bin(
      audio = origSoundFile,
      start_time = if (bt > 0) bt else NULL,
      end_time = if (et > 0) et else NULL,
      channels = 1
    )
    sample_rate <- attr(audio_data, "sample_rate")

    # Convert INT32 to float64
    audio_float <- as.numeric(audio_data) / 2147483647.0

    # Call C++ GFM-IAIF
    result <- gfmiaif_cpp(audio_float, as.integer(sample_rate),
                           window_shift_sec, window_size_sec,
                           nv, ng, d, window_type)

    av_matrix <- result$av
    ag_matrix <- result$ag
    al_matrix <- result$al
    timestamps <- as.numeric(result$timestamps)
    frame_rate_hz <- as.numeric(result$sample_rate_hz)
    n_frames <- as.integer(result$n_frames)

    if (verbose) {
      message(sprintf("  Processed %d frames (%.2f seconds)", n_frames, max(timestamps)))
    }

    # Adjust timestamps if centerTime is FALSE
    if (!centerTime) {
      timestamps <- timestamps - (window_size_sec / 2.0)
    }

    # Offset timestamps by beginTime (since audio was loaded from bt)
    timestamps <- timestamps + bt
    # Clamp to non-negative (floating point can make first frame slightly negative)
    timestamps <- pmax(timestamps, 0.0)

    # Create AsspDataObj
    outDataObj <- list()
    attr(outDataObj, "trackFormats") <- c(
      rep("REAL64", nv + 1),
      rep("REAL64", ng + 1),
      rep("REAL64", 2)
    )
    attr(outDataObj, "sampleRate") <- as.numeric(frame_rate_hz)
    attr(outDataObj, "origFreq") <- as.numeric(sample_rate)
    attr(outDataObj, "startTime") <- as.numeric(timestamps[1])
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(n_frames)
    class(outDataObj) <- "AsspDataObj"

    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2)

    # Add vocal tract tracks
    for (j in 0:nv) {
      track_name <- sprintf("av_%d", j)
      outDataObj <- addTrack(outDataObj, track_name,
                            as.matrix(av_matrix[, j + 1]),
                            "REAL64")
    }

    # Add glottis tracks
    for (j in 0:ng) {
      track_name <- sprintf("ag_%d", j)
      outDataObj <- addTrack(outDataObj, track_name,
                            as.matrix(ag_matrix[, j + 1]),
                            "REAL64")
    }

    # Add lip radiation tracks
    for (j in 0:1) {
      track_name <- sprintf("al_%d", j)
      outDataObj <- addTrack(outDataObj, track_name,
                            as.matrix(al_matrix[, j + 1]),
                            "REAL64")
    }

    if (!is.AsspDataObj(outDataObj)) {
      stop("The AsspDataObj created by trk_gfmiaif is invalid.", call. = FALSE)
    }

    # Output file path
    ssff_file <- sub("\\.[^.]+$", paste0(".", explicitExt), origSoundFile)
    if (!is.null(outputDirectory)) {
      ssff_file <- file.path(outputDirectory, basename(ssff_file))
    }
    attr(outDataObj, "filePath") <- as.character(ssff_file)

    if (toFile) {
      write.AsspDataObj(dobj = outDataObj, file = ssff_file)
      outListOfFiles <- c(outListOfFiles, TRUE)
      if (verbose) message(sprintf("  -> Written to: %s", ssff_file))
    }
  }

  if (toFile) return(length(outListOfFiles)) else return(outDataObj)
}

# Set function attributes
attr(trk_gfmiaif, "ext") <- "gfm"
attr(trk_gfmiaif, "tracks") <- function(nv = 48, ng = 3) {
  c(paste0("av_", 0:nv), paste0("ag_", 0:ng), paste0("al_", 0:1))
}
attr(trk_gfmiaif, "outputType") <- "SSFF"

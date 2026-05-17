#' Decompose speech into vocal tract, glottis, and lip radiation filters using GFM-IAIF
#'
#' Performs source-filter separation using the Glottal Flow Model-based Iterative
#' Adaptive Inverse Filtering (GFM-IAIF) algorithm \insertCite{perrotin2019gfmiaif}{superassp}.
#' GFM-IAIF extends classical IAIF with a wide-band glottis model that captures both
#' glottal formant and spectral tilt, yielding LP coefficient tracks for the vocal
#' tract, glottis, and lip radiation filters. Prefer this over \code{trk_covarep_iaif()}
#' when you need frame-level LP coefficients rather than a sample-domain glottal waveform.
#'
#' @param listOfFiles Character vector of audio file paths. Any format supported by
#'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
#' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
#' @param centerTime Logical. If \code{TRUE}, timestamps refer to window centres;
#'   if \code{FALSE}, to window starts. Default \code{FALSE}.
#' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
#' @param windowShift Numeric. Frame shift in milliseconds; sets output frame rate
#'   (1000 / windowShift Hz). Default 10.0 ms.
#' @param windowSize Numeric. Analysis window length in milliseconds. Must be large enough
#'   for \code{nv + 1} samples at the audio sample rate. Default 32.0 ms.
#' @param nv Integer. Vocal tract LPC order (1–100). Higher values capture more spectral
#'   detail but increase compute. Default 48.
#' @param ng Integer. Glottis LPC order. The value 3 is recommended by the original
#'   authors; departing from it may degrade estimates. Default 3.
#' @param d Numeric. Leaky integration coefficient for lip radiation (0.9–0.999).
#'   Default 0.99.
#' @param window Character. Analysis window type: \code{"HANN"}, \code{"HAMMING"}, or
#'   \code{"BLACKMAN"}. Default \code{"HANN"}.
#' @param explicitExt Character. Output file extension. Default \code{"gfm"}.
#' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
#'   writes alongside the input file.
#' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
#'   count written. If \code{FALSE}, return an \code{AsspDataObj} (single file only).
#'   Default \code{TRUE}.
#' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
#'
#' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
#'   \describe{
#'     \item{\code{av_0} … \code{av_\{nv\}}}{REAL64, vocal tract LP coefficients,
#'       n_frames × (nv + 1). Dimensionless.}
#'     \item{\code{ag_0} … \code{ag_\{ng\}}}{REAL64, glottis LP coefficients,
#'       n_frames × (ng + 1). Dimensionless.}
#'     \item{\code{al_0}, \code{al_1}}{REAL64, lip radiation LP coefficients,
#'       n_frames × 2. Dimensionless.}
#'   }
#'   Frame rate: \code{1000 / windowShift} Hz (default 100 Hz).
#'   If \code{toFile = TRUE}: integer count of files written.
#'
#' @details
#' \code{ng = 3} is strongly recommended. When \code{toFile = FALSE}, only single-file
#' input is permitted. The \code{ng} warning is issued automatically for other values.
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
#' @references \insertAllCited{}
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

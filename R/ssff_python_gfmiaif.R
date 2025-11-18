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
#' GFM-IAIF performs iterative refinement to separate these components through inverse filtering,
#' providing frame-by-frame linear prediction coefficients for each filter.
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
#'   Typical range: 24-48 for speech at 16kHz.
#' @param ng Glottis LPC order. Default is 3., and highly recommended as the
#'   algorithm is specifically designed assuming a 3rd-order filter captures glottis-related timbre
#'   variations (tenseness, effort, breathiness).
#' @param d Leaky integration coefficient for lip radiation modeling. Default is 0.99.
#'   Range: 0.95-0.99. Higher values emphasize high frequencies more strongly.
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
#'   If \code{toFile = FALSE}: Returns an AsspDataObj with tracks (columns depend on nv and ng):
#'   \itemize{
#'     \item \strong{av_0 to av_<nv>}: Vocal tract LP coefficients (REAL64 format), nv+1 columns
#'     \item \strong{ag_0 to ag_<ng>}: Glottis LP coefficients (REAL64 format), ng+1 columns
#'     \item \strong{al_0, al_1}: Lip radiation LP coefficients (REAL64 format), 2 columns
#'   }
#'
#'   Total tracks: (nv+1) + (ng+1) + 2. With defaults (nv=48, ng=3): 49 + 4 + 2 = 55 tracks.
#'
#' @details
#' \strong{GFM-IAIF Algorithm:}
#'
#' The algorithm performs six steps:
#' \enumerate{
#'   \item Pre-frame addition: Adds mean-normalized ramp to reduce edge effects
#'   \item Lip radiation cancellation: Applies leaky integration (1/\[1 - d·z^(-1)\])
#'   \item Gross glottis estimation: Iterative 1st-order LPC (ng iterations)
#'   \item Gross vocal tract estimation: nv-order LPC after glottis removal
#'   \item Fine glottis estimation: ng-order LPC after vocal tract removal
#'   \item Fine vocal tract estimation: Final nv-order LPC refinement
#' }
#'
#' \strong{Installation Requirements:}
#'
#' GFM-IAIF requires Python with numpy and scipy installed. Install with:
#'
#' \code{install_gfmiaif()}
#'
#' Or manually:
#'
#' \code{reticulate::py_install(c("numpy", "scipy"))}
#'
#' For maximum performance (5-10x speedup), also install numba:
#'
#' \code{reticulate::py_install("numba")}
#'
#' \strong{Performance:}
#'
#' \itemize{
#'   \item Single frame (512 samples): ~0.25 ms
#'   \item Real-time factor (16kHz): ~130x
#'   \item With Numba JIT: additional 5-10x speedup for LPC
#' }
#'
#' \strong{Parameter Guidelines:}
#'
#' \itemize{
#'   \item \strong{windowSize}: Should be 20-40ms for speech. Longer windows provide
#'     better frequency resolution but worse time resolution.
#'   \item \strong{nv}: Higher order = more detailed formant structure, but may overfit.
#'     Use nv=24 for low sample rates, nv=48 for 16kHz.
#'   \item \strong{ng}: Keep at ng=3 unless you have specific reasons. This is a core
#'     design assumption of GFM-IAIF.
#'   \item \strong{d}: Values 0.97-0.99 are standard. Lower values reduce high-frequency emphasis.
#' }
#'
#' @examples
#' \dontrun{
#' # Install GFM-IAIF first
#' install_gfmiaif()
#'
#' # Basic usage with default parameters
#' result <- trk_gfmiaif("speech.wav", toFile = FALSE)
#'
#' # Extract vocal tract coefficients (first 49 columns)
#' vocal_tract <- result$av_0  # Plus av_1, av_2, ..., av_48
#'
#' # Lower vocal tract order for faster processing
#' result <- trk_gfmiaif("speech.wav",
#'                      nv = 24,
#'                      windowShift = 5.0,
#'                      toFile = FALSE)
#'
#' # Batch processing to SSFF files
#' files <- c("file1.wav", "file2.wav", "file3.mp3")
#' trk_gfmiaif(files, toFile = TRUE)
#'
#' # Time windowing
#' result <- trk_gfmiaif("long_recording.wav",
#'                      beginTime = 1.0,
#'                      endTime = 3.0,
#'                      toFile = FALSE)
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

  # Validate parameters
  nv <- as.integer(nv)
  ng <- as.integer(ng)

  if (nv < 1 || nv > 100) {
    stop("nv must be between 1 and 100", call. = FALSE)
  }

  if (ng < 1 || ng > 10) {
    stop("ng must be between 1 and 10", call. = FALSE)
  }

  if (ng != 3 && verbose) {
    warning("GFM-IAIF is designed for ng=3. Using ng=", ng, " may not be optimal.",
            call. = FALSE)
  }

  if (d < 0.9 || d > 0.999) {
    stop("d must be between 0.9 and 0.999", call. = FALSE)
  }

  # Calculate minimum frame length needed
  min_frame_length <- nv + 2  # Safety margin
  min_window_size_ms <- (min_frame_length / 16.0)  # Assuming 16kHz

  if (windowSize < min_window_size_ms && verbose) {
    warning(sprintf("windowSize (%.1f ms) may be too small for nv=%d. Consider >= %.1f ms",
                   windowSize, nv, min_window_size_ms),
            call. = FALSE)
  }

  # Check GFM-IAIF Python module availability
  if (!reticulate::py_module_available("numpy")) {
    stop(
      "trk_gfmiaif() requires numpy.\n\n",
      "Install with: install_gfmiaif()\n",
      "Or manually: reticulate::py_install('numpy')",
      call. = FALSE
    )
  }

  if (!reticulate::py_module_available("scipy")) {
    stop(
      "trk_gfmiaif() requires scipy.\n\n",
      "Install with: install_gfmiaif()\n",
      "Or manually: reticulate::py_install('scipy')",
      call. = FALSE
    )
  }

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

  # Check that all files exist
  filesEx <- file.exists(listOfFiles)
  if (!all(filesEx)) {
    filesNotExists <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s): ", paste(filesNotExists, collapse = ", "),
         call. = FALSE)
  }

  # Get paths to Python module
  python_dir <- system.file("python", package = "superassp")
  gfmiaif_module_path <- file.path(python_dir, "gfmiaif")

  if (!dir.exists(gfmiaif_module_path)) {
    stop("GFM-IAIF Python module not found in package installation.",
         call. = FALSE)
  }

  # Add module to Python path
  reticulate::py_run_string(sprintf("import sys; sys.path.insert(0, '%s')",
                                   gsub("\\\\", "/", python_dir)))

  # Import modules
  gfmiaif <- reticulate::import("gfmiaif", delay_load = FALSE)
  np <- reticulate::import("numpy", delay_load = FALSE)

  # Map window type to lowercase
  window_type <- tolower(window)
  if (!window_type %in% c("hann", "hamming", "blackman")) {
    stop("window must be one of: HANN, HAMMING, BLACKMAN", call. = FALSE)
  }

  # Convert window shift/size to seconds
  window_shift_sec <- windowShift / 1000.0
  window_size_sec <- windowSize / 1000.0

  # Vector of successfully processed files
  outListOfFiles <- c()

  # Process each file
  for (i in seq_len(nrow(fileBeginEnd))) {
    origSoundFile <- normalizePath(fileBeginEnd[i, "listOfFiles"], mustWork = TRUE)
    bt <- fileBeginEnd[i, "beginTime"]
    et <- fileBeginEnd[i, "endTime"]

    if (verbose) {
      message(sprintf("Processing file %d/%d: %s", i, nrow(fileBeginEnd), basename(origSoundFile)))
    }

    # Load audio using av package
    audio_data <- av::read_audio_bin(
      audio = origSoundFile,
      start_time = if (bt > 0) bt else NULL,
      end_time = if (et > 0) et else NULL,
      channels = 1  # GFM-IAIF requires mono
    )

    # Extract audio properties
    sample_rate <- attr(audio_data, "sample_rate")

    # Convert to float64 numpy array
    # av::read_audio_bin returns INT32, normalize to [-1, 1]
    audio_float <- as.numeric(audio_data) / 2147483647.0  # INT32_MAX
    audio_array <- np$array(audio_float, dtype = "float64")

    # Run GFM-IAIF frame-based processing
    result <- gfmiaif$gfmiaif_frame_based(
      audio = audio_array,
      sample_rate = as.integer(sample_rate),
      window_shift = as.numeric(window_shift_sec),
      window_size = as.numeric(window_size_sec),
      nv = as.integer(nv),
      ng = as.integer(ng),
      d = as.numeric(d),
      window_type = window_type
    )

    # Extract results
    av_matrix <- result$av  # shape: (n_frames, nv+1)
    ag_matrix <- result$ag  # shape: (n_frames, ng+1)
    al_matrix <- result$al  # shape: (n_frames, 2)
    timestamps <- as.numeric(result$timestamps)
    frame_rate_hz <- as.numeric(result$sample_rate_hz)
    n_frames <- as.integer(result$n_frames)

    if (verbose) {
      message(sprintf("  Processed %d frames (%.2f seconds)", n_frames, max(timestamps)))
    }

    # Adjust timestamps if centerTime is FALSE (use start times)
    if (!centerTime) {
      frame_shift_sec <- as.numeric(result$frame_shift_sec)
      frame_size_sec <- as.numeric(result$frame_size_sec)
      timestamps <- timestamps - (frame_size_sec / 2.0)
    }

    # Create AsspDataObj
    outDataObj <- list()
    attr(outDataObj, "trackFormats") <- c(
      rep("REAL64", nv + 1),  # av coefficients
      rep("REAL64", ng + 1),  # ag coefficients
      rep("REAL64", 2)        # al coefficients
    )
    attr(outDataObj, "sampleRate") <- as.numeric(frame_rate_hz)
    attr(outDataObj, "origFreq") <- as.numeric(sample_rate)
    attr(outDataObj, "startTime") <- as.numeric(timestamps[1])
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(n_frames)
    class(outDataObj) <- "AsspDataObj"

    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2)  # binary

    # Add vocal tract tracks (av_0 to av_<nv>)
    for (j in 0:nv) {
      track_name <- sprintf("av_%d", j)
      outDataObj <- addTrack(outDataObj, track_name,
                            as.matrix(av_matrix[, j + 1]),  # R is 1-indexed
                            "REAL64")
    }

    # Add glottis tracks (ag_0 to ag_<ng>)
    for (j in 0:ng) {
      track_name <- sprintf("ag_%d", j)
      outDataObj <- addTrack(outDataObj, track_name,
                            as.matrix(ag_matrix[, j + 1]),
                            "REAL64")
    }

    # Add lip radiation tracks (al_0, al_1)
    for (j in 0:1) {
      track_name <- sprintf("al_%d", j)
      outDataObj <- addTrack(outDataObj, track_name,
                            as.matrix(al_matrix[, j + 1]),
                            "REAL64")
    }

    # Validate AsspDataObj
    if (!is.AsspDataObj(outDataObj)) {
      stop("The AsspDataObj created by trk_gfmiaif is invalid.", call. = FALSE)
    }

    # Prepare output file path
    ssff_file <- sub("\\.[^.]+$", paste0(".", explicitExt), origSoundFile)
    if (!is.null(outputDirectory)) {
      ssff_file <- file.path(outputDirectory, basename(ssff_file))
    }

    attr(outDataObj, "filePath") <- as.character(ssff_file)

    # Write to file or return
    if (toFile) {
      write.AsspDataObj(dobj = outDataObj, file = ssff_file)
      outListOfFiles <- c(outListOfFiles, TRUE)

      if (verbose) {
        message(sprintf("  → Written to: %s", ssff_file))
      }
    }
  }

  # Return based on toFile flag
  if (toFile) {
    return(length(outListOfFiles))
  } else {
    return(outDataObj)
  }
}

# Set function attributes
attr(trk_gfmiaif, "ext") <- "gfm"
attr(trk_gfmiaif, "tracks") <- function(nv = 48, ng = 3) {
  c(
    paste0("av_", 0:nv),
    paste0("ag_", 0:ng),
    paste0("al_", 0:1)
  )
}
attr(trk_gfmiaif, "outputType") <- "SSFF"

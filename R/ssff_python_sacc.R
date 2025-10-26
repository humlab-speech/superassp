#' Compute pitch using SAcC (Subband Autocorrelation Classification)
#'
#' SAcC is a robust pitch tracking algorithm developed by Dan Ellis at Columbia University
#' \insertCite{ellis2014sacc}{superassp}. The algorithm works by filtering audio into 24 subbands,
#' computing autocorrelation features, projecting onto principal components, and classifying
#' with a neural network followed by Viterbi decoding for temporal continuity.
#'
#' SAcC is particularly robust to noise and handles speech well. It was originally designed
#' for telephone speech analysis (8kHz) and is effective for both clean and noisy conditions.
#'
#' The algorithm processes audio at **8kHz** (automatic resampling from other rates) with
#' **10ms frame shifts** (100 Hz frame rate).
#'
#' @param listOfFiles A character vector of file paths to audio files, or an AVAudio S7 object.
#' @param beginTime The start time of the section of the sound file that should be processed (in seconds).
#'   Default is 0.0 (start of file).
#' @param endTime The end time of the section of the sound file that should be processed (in seconds).
#'   Default is 0.0 (end of file).
#' @param hmm_vp HMM voicing penalty parameter (0-1). Higher values penalize unvoiced states.
#'   Default is 0.9. Values closer to 1 result in more voiced frames.
#' @param dither_level Level of dithering noise to add to avoid digital zeros.
#'   Default is 1e-3 (0.001 times signal standard deviation). Set to 0 to disable.
#' @param explicitExt The file extension for the output SSFF file. Default is "sacc".
#' @param outputDirectory The directory where output files should be written.
#'   If NULL (default), files are written to the same directory as the input files.
#' @param toFile If TRUE (default), write results to SSFF files and return file count.
#'   If FALSE, return AsspDataObj for single file processing.
#' @param verbose If TRUE (default), print progress messages.
#'
#' @return
#'   If \code{toFile = TRUE}: Returns the number of successfully processed files.
#'   If \code{toFile = FALSE}: Returns an AsspDataObj with two tracks:
#'   \itemize{
#'     \item \strong{F0}: Fundamental frequency in Hz (REAL32 format, 0 = unvoiced)
#'     \item \strong{prob_voiced}: Probability of voicing 0-1 (REAL32 format)
#'   }
#'
#' @details
#' \strong{SAcC Algorithm Stages:}
#' \enumerate{
#'   \item \strong{Subband Filtering}: Audio is filtered into 24 frequency bands using
#'     gammatone filters (ERB-spaced, 100-800 Hz range)
#'   \item \strong{Autocorrelation}: Each subband is analyzed with normalized autocorrelation
#'     over 25ms windows with 10ms shifts
#'   \item \strong{PCA}: Autocorrelation features are projected onto 10 principal components
#'     per subband (240 features total)
#'   \item \strong{Neural Network}: Multi-layer perceptron classifies frames into 67 pitch
#'     bins + 1 unvoiced state
#'   \item \strong{Viterbi Decoding}: HMM smooths the pitch track over time for temporal
#'     continuity
#' }
#'
#' \strong{SAcC Specifications:}
#' \itemize{
#'   \item \strong{Sample Rate}: Processes at 8kHz (automatic resampling)
#'   \item \strong{Frame Rate}: 10ms frames (100 Hz)
#'   \item \strong{Frequency Range}: ~80-500 Hz (67 pitch bins)
#'   \item \strong{Subbands}: 24 gammatone filters (6 bands/octave)
#'   \item \strong{Features}: 240 PCA features (24 subbands × 10 components)
#'   \item \strong{Classifier}: MLP with 100 hidden units
#' }
#'
#' \strong{Installation Requirements:}
#'
#' SAcC requires Python with numpy, scipy, and soundfile packages. Install with:
#'
#' \code{install_sacc()}
#'
#' Or manually:
#'
#' \code{reticulate::py_install(c("numpy", "scipy", "soundfile"))}
#'
#' \strong{Performance Notes:}
#'
#' SAcC is computationally intensive compared to simpler algorithms like RAPT or SWIPE:
#' \itemize{
#'   \item For a 3-second audio file: ~500-800ms processing time
#'   \item Batch processing benefits from reduced overhead
#'   \item Most time spent in subband filtering and autocorrelation
#' }
#'
#' @examples
#' \dontrun{
#' # Install SAcC first
#' install_sacc()
#'
#' # Basic usage
#' result <- trk_sacc("speech.wav", toFile = FALSE)
#'
#' # Extract pitch and voicing probability
#' f0 <- result$F0  # Hz, 0 = unvoiced
#' prob_voiced <- result$prob_voiced  # 0-1
#'
#' # Time windowing
#' result <- trk_sacc("long_recording.wav",
#'                   beginTime = 1.0,
#'                   endTime = 3.0,
#'                   toFile = FALSE)
#'
#' # Adjust voicing penalty (more permissive)
#' result <- trk_sacc("noisy_speech.wav",
#'                   hmm_vp = 0.7,
#'                   toFile = FALSE)
#'
#' # Batch processing to SSFF files
#' files <- c("file1.wav", "file2.wav", "file3.mp3")
#' trk_sacc(files, toFile = TRUE)
#'
#' # S7 AVAudio dispatch (automatic)
#' audio <- read_avaudio("speech.wav", sample_rate = 8000)
#' result <- trk_sacc(audio, toFile = FALSE)
#' }
#'
#' @references
#' Ellis, D. P. W. (2014). Subband Autocorrelation Classification (SAcC) pitch tracker.
#' \url{http://labrosa.ee.columbia.edu/projects/SAcC/}
#'
#' @export
trk_sacc <- function(listOfFiles,
                     beginTime = 0.0,
                     endTime = 0.0,
                     hmm_vp = 0.9,
                     dither_level = 1e-3,
                     explicitExt = "sacc",
                     outputDirectory = NULL,
                     toFile = TRUE,
                     verbose = TRUE) {

  # Validate inputs
  if (length(listOfFiles) > 1 && !toFile) {
    stop("length(listOfFiles) is > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.",
         call. = FALSE)
  }

  # Check SAcC availability
  if (!sacc_available()) {
    stop(
      "trk_sacc() requires Python with numpy, scipy, and soundfile.\n\n",
      "Install with: install_sacc()\n",
      "Or manually: reticulate::py_install(c('numpy', 'scipy', 'soundfile'))\n",
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

  # Add SAcC Python module to path
  sacc_path <- system.file("python", "calc_sbpca", "python", package = "superassp")
  if (sacc_path == "") {
    stop("Unable to find SAcC Python module in package installation", call. = FALSE)
  }

  # Add to Python path
  sys <- reticulate::import("sys", delay_load = FALSE)
  if (!sacc_path %in% sys$path) {
    sys$path <- c(sacc_path, sys$path)
  }

  # Import required modules
  np <- reticulate::import("numpy", delay_load = FALSE)
  scipy <- reticulate::import("scipy", delay_load = FALSE)
  SAcC_module <- reticulate::import("SAcC", delay_load = FALSE)

  # Create SAcC configuration
  config <- SAcC_module$default_config()
  config$hmm_vp <- as.numeric(hmm_vp)
  config$dither_level <- as.numeric(dither_level)
  config$write_time <- 1L  # Include timestamps
  config$write_pitch <- 1L  # Include pitch
  config$write_pvx <- 1L  # Include prob(voiced)
  config$verbose <- 0L  # Suppress Python verbosity

  # Create SAcC extractor object
  sacc_extractor <- SAcC_module$SAcC(config)

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

    tryCatch({
      # Load audio using av package (handles all formats, time windowing)
      audio_data <- av::read_audio_bin(
        audio = origSoundFile,
        start_time = if (bt > 0) bt else NULL,
        end_time = if (et > 0) et else NULL,
        channels = 1,  # SAcC requires mono
        sample_rate = 8000  # SAcC processes at 8kHz
      )

      # Extract audio properties
      sample_rate <- attr(audio_data, "sample_rate")

      # Verify sample rate is 8kHz (av should have resampled)
      if (sample_rate != 8000) {
        warning(sprintf("Expected 8kHz audio, got %d Hz. Results may be inaccurate.", sample_rate))
      }

      # Convert to float64 numpy array (required by SAcC)
      # av::read_audio_bin returns INT32, normalize to [-1, 1]
      audio_float <- as.numeric(audio_data) / 2147483647.0  # INT32_MAX
      audio_array <- np$array(audio_float, dtype = "float64")

      # Run SAcC pitch tracking
      # Returns numpy array with columns: [time, pitch_hz, prob_voiced]
      features <- sacc_extractor$sacc(audio_array, as.integer(sample_rate))

      # Extract columns
      # Note: Python is 0-indexed, R is 1-indexed
      timestamps <- as.numeric(features[, 1])  # time in seconds
      pitch_hz <- as.numeric(features[, 2])    # F0 in Hz
      prob_voiced <- as.numeric(features[, 3]) # P(voiced)

      # Create data frame
      inTable <- data.frame(
        F0 = pitch_hz,          # Keep as REAL32
        prob_voiced = prob_voiced  # REAL32 (0-1)
      )

      # Calculate SSFF sample rate (frame rate in Hz)
      # SAcC uses 10ms frames = 100 Hz
      ssff_sample_rate <- 100.0

      # Adjust if we have actual timestamps
      if (length(timestamps) > 1) {
        frame_shift_sec <- timestamps[2] - timestamps[1]
        ssff_sample_rate <- 1.0 / frame_shift_sec
      }

      # Create AsspDataObj
      outDataObj <- list()
      attr(outDataObj, "trackFormats") <- c("REAL32", "REAL32")
      attr(outDataObj, "sampleRate") <- as.numeric(ssff_sample_rate)
      attr(outDataObj, "origFreq") <- as.numeric(sample_rate)
      attr(outDataObj, "startTime") <- as.numeric(if (bt > 0) bt else 0.0)
      attr(outDataObj, "startRecord") <- as.integer(1)
      attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
      class(outDataObj) <- "AsspDataObj"

      AsspFileFormat(outDataObj) <- "SSFF"
      AsspDataFormat(outDataObj) <- as.integer(2)  # binary

      # Add F0 track
      outDataObj <- addTrack(outDataObj, "F0",
                            as.matrix(inTable$F0),
                            "REAL32")

      # Add prob_voiced track
      outDataObj <- addTrack(outDataObj, "prob_voiced",
                            as.matrix(inTable$prob_voiced),
                            "REAL32")

      # Validate AsspDataObj
      if (!is.AsspDataObj(outDataObj)) {
        stop("The AsspDataObj created by trk_sacc is invalid.", call. = FALSE)
      }

      # Prepare output file path
      ssff_file <- sub("\\.[^.]+$", paste0(".", explicitExt), origSoundFile)
      if (!is.null(outputDirectory)) {
        ssff_file <- file.path(outputDirectory, basename(ssff_file))
      }

      attr(outDataObj, "filePath") <- as.character(ssff_file)

      # Write to file or accumulate
      if (toFile) {
        write.AsspDataObj(dobj = outDataObj, file = ssff_file)
        outListOfFiles <- c(outListOfFiles, TRUE)

        if (verbose) {
          message(sprintf("  → Written to: %s", ssff_file))
        }
      }

    }, error = function(e) {
      warning(sprintf("Error processing %s: %s", basename(origSoundFile), e$message))
      if (toFile) {
        outListOfFiles <- c(outListOfFiles, FALSE)
      }
    })
  }

  # Return based on toFile flag
  if (toFile) {
    return(sum(outListOfFiles))  # Count of successful files
  } else {
    return(outDataObj)
  }
}

# Set function attributes (for compatibility with existing infrastructure)
attr(trk_sacc, "ext") <- "sacc"
attr(trk_sacc, "tracks") <- c("F0", "prob_voiced")
attr(trk_sacc, "outputType") <- "SSFF"

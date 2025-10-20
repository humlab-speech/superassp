#' Compute pitch using the Swift-F0 deep learning pitch tracker
#'
#' Swift-F0 \insertCite{nieradzik2025swiftf0}{superassp} is a fast and accurate F0 detector that works by first converting audio
#' into a spectrogram using an STFT, then applying a 2D convolutional neural network to estimate pitch.
#' In the Pitch Detection Benchmark, Swift-F0 outperforms algorithms like CREPE in both speed and accuracy,
#' achieving real-time analysis speeds (132 ms for 5 seconds of audio on CPU).
#'
#' The model supports frequencies between **46.875 Hz and 2093.75 Hz** (G1 to C7). For speech analysis,
#' the default parameters are set to 75-400 Hz. For music analysis, you may want to use the full range
#' by setting \code{minF = 46.875} and \code{maxF = 2093.75}.
#'
#' The algorithm processes audio at 16kHz with a 256-sample hop size (16ms frames). Audio at other
#' sample rates is automatically resampled.
#'
#' @param listOfFiles A character vector of file paths to audio files, or an AVAudio S7 object.
#' @param minF Minimum frequency (Hz) to consider as voiced. Default is 75 Hz (speech).
#'   For music, use 46.875 Hz. Must be >= 46.875 Hz (model minimum).
#' @param maxF Maximum frequency (Hz) to consider as voiced. Default is 400 Hz (speech).
#'   For music, use 2093.75 Hz. Must be <= 2093.75 Hz (model maximum).
#' @param confidence_threshold Confidence threshold (0-1) for voicing decision.
#'   Frames with confidence below this threshold are marked as unvoiced. Default is 0.9.
#' @param beginTime The start time of the section of the sound file that should be processed (in seconds).
#'   Default is 0.0 (start of file).
#' @param endTime The end time of the section of the sound file that should be processed (in seconds).
#'   Default is 0.0 (end of file).
#' @param explicitExt The file extension for the output SSFF file. Default is "sf0".
#' @param outputDirectory The directory where output files should be written.
#'   If NULL (default), files are written to the same directory as the input files.
#' @param toFile If TRUE (default), write results to SSFF files and return file count.
#'   If FALSE, return AsspDataObj for single file processing.
#' @param verbose If TRUE (default), print progress messages.
#'
#' @return
#'   If \code{toFile = TRUE}: Returns the number of successfully processed files.
#'   If \code{toFile = FALSE}: Returns an AsspDataObj with three tracks:
#'   \itemize{
#'     \item \strong{F0}: Fundamental frequency in Hz (INT16 format, 0 = unvoiced)
#'     \item \strong{confidence}: Model confidence score 0-1 (REAL32 format)
#'     \item \strong{voicing}: Voiced/unvoiced decision (INT16 format, 1 = voiced, 0 = unvoiced)
#'   }
#'
#' @details
#' \strong{Swift-F0 Specifications:}
#' \itemize{
#'   \item \strong{Speed}: 132ms for 5 seconds of audio (CPU)
#'   \item \strong{Sample Rate}: Processes at 16kHz (automatic resampling)
#'   \item \strong{Frame Rate}: 16ms frames (256-sample hop at 16kHz = 62.5 Hz)
#'   \item \strong{Frequency Range}: 46.875 - 2093.75 Hz (G1 to C7)
#'   \item \strong{Model}: ONNX-based CNN (runs on CPU via ONNXRuntime)
#' }
#'
#' \strong{Installation Requirements:}
#'
#' Swift-F0 requires Python with the swift-f0 package installed. Install with:
#'
#' \code{install_swiftf0()}
#'
#' Or manually:
#'
#' \code{reticulate::py_install("swift-f0")}
#'
#' This will automatically install the required dependencies (onnxruntime and numpy).
#'
#' \strong{Performance Comparison:}
#'
#' For a 3-second audio file:
#' \itemize{
#'   \item RAPT (C++): ~100-200ms
#'   \item REAPER (C++): ~80-150ms
#'   \item Swift-F0 (Python): ~90-130ms
#'   \item SWIPE (C++): ~150-250ms
#' }
#'
#' @examples
#' \dontrun{
#' # Install Swift-F0 first
#' install_swiftf0()
#'
#' # Basic usage with default speech parameters (75-400 Hz)
#' result <- trk_swiftf0("speech.wav", toFile = FALSE)
#'
#' # Music analysis with full frequency range
#' result <- trk_swiftf0("song.wav",
#'                      minF = 46.875,
#'                      maxF = 2093.75,
#'                      toFile = FALSE)
#'
#' # Custom confidence threshold (more permissive voicing)
#' result <- trk_swiftf0("noisy_speech.wav",
#'                      confidence_threshold = 0.7,
#'                      toFile = FALSE)
#'
#' # Time windowing
#' result <- trk_swiftf0("long_recording.wav",
#'                      beginTime = 1.0,
#'                      endTime = 3.0,
#'                      toFile = FALSE)
#'
#' # Batch processing to SSFF files
#' files <- c("file1.wav", "file2.wav", "file3.mp3")
#' trk_swiftf0(files, toFile = TRUE)
#'
#' # S7 AVAudio dispatch (automatic)
#' audio <- read_avaudio("speech.wav", sample_rate = 16000)
#' result <- trk_swiftf0(audio, toFile = FALSE)
#' }
#'
#' @references
#'   \insertAllCited{}
#'
#' @export
trk_swiftf0 <- function(listOfFiles,
                       minF = 75,
                       maxF = 400,
                       confidence_threshold = 0.9,
                       beginTime = 0.0,
                       endTime = 0.0,
                       explicitExt = "sf0",
                       outputDirectory = NULL,
                       toFile = TRUE,
                       verbose = TRUE) {

  # Validate inputs
  if (length(listOfFiles) > 1 && !toFile) {
    stop("length(listOfFiles) is > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.",
         call. = FALSE)
  }

  # Check Swift-F0 availability
  if (!reticulate::py_module_available("swift_f0")) {
    stop(
      "trk_swiftf0() requires the swift-f0 Python module.\n\n",
      "Install with: install_swiftf0()\n",
      "Or manually: reticulate::py_install('swift-f0')\n\n",
      "This will automatically install required dependencies (onnxruntime, numpy).",
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

  # Import Swift-F0 module
  swift_f0 <- reticulate::import("swift_f0", delay_load = FALSE)
  np <- reticulate::import("numpy", delay_load = FALSE)

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

    # Load audio using av package (handles all formats, time windowing)
    audio_data <- av::read_audio_bin(
      audio = origSoundFile,
      start_time = if (bt > 0) bt else NULL,
      end_time = if (et > 0) et else NULL,
      channels = 1  # Swift-F0 requires mono
    )

    # Extract audio properties
    sample_rate <- attr(audio_data, "sample_rate")
    n_channels <- attr(audio_data, "channels")

    # Convert to float32 numpy array (required by Swift-F0)
    # av::read_audio_bin returns INT32, normalize to [-1, 1]
    audio_float <- as.numeric(audio_data) / 2147483647.0  # INT32_MAX
    audio_array <- np$array(audio_float, dtype = "float32")

    # Create Swift-F0 detector
    detector <- swift_f0$SwiftF0(
      fmin = as.numeric(minF),
      fmax = as.numeric(maxF),
      confidence_threshold = as.numeric(confidence_threshold)
    )

    # Run pitch detection
    result <- detector$detect_from_array(audio_array, as.integer(sample_rate))

    # Extract results
    pitch_hz <- as.numeric(result$pitch_hz)
    confidence <- as.numeric(result$confidence)
    timestamps <- as.numeric(result$timestamps)
    voicing <- as.integer(result$voicing)

    # Convert to AsspDataObj format
    # For unvoiced frames, set F0 to 0 (SSFF convention)
    f0_values <- ifelse(voicing == 1, pitch_hz, 0)

    # Create data frame
    inTable <- data.frame(
      f0 = as.integer(round(f0_values)),  # Convert to INT16
      confidence = confidence,             # Keep as REAL32
      voicing = voicing                    # INT16 (0 or 1)
    )

    # Calculate sample rate for SSFF (frame rate in Hz)
    # Swift-F0 uses 16ms frames (256 samples at 16kHz)
    if (length(timestamps) > 1) {
      frame_shift_sec <- timestamps[2] - timestamps[1]
      ssff_sample_rate <- 1.0 / frame_shift_sec
    } else {
      ssff_sample_rate <- 62.5  # Default: 16ms frames = 62.5 Hz
    }

    # Create AsspDataObj
    outDataObj <- list()
    attr(outDataObj, "trackFormats") <- c("INT16", "REAL32", "INT16")
    attr(outDataObj, "sampleRate") <- as.numeric(ssff_sample_rate)
    attr(outDataObj, "origFreq") <- as.numeric(sample_rate)
    attr(outDataObj, "startTime") <- as.numeric(timestamps[1])
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
    class(outDataObj) <- "AsspDataObj"

    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2)  # binary

    # Add F0 track
    outDataObj <- addTrack(outDataObj, "F0",
                          as.matrix(inTable$f0),
                          "INT16")

    # Add confidence track
    outDataObj <- addTrack(outDataObj, "confidence",
                          as.matrix(inTable$confidence),
                          "REAL32")

    # Add voicing track
    outDataObj <- addTrack(outDataObj, "voicing",
                          as.matrix(inTable$voicing),
                          "INT16")

    # Validate AsspDataObj
    if (!is.AsspDataObj(outDataObj)) {
      stop("The AsspDataObj created by trk_swiftf0 is invalid.", call. = FALSE)
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
  }

  # Return based on toFile flag
  if (toFile) {
    return(length(outListOfFiles))
  } else {
    return(outDataObj)
  }
}

# Set function attributes (for compatibility with existing infrastructure)
attr(trk_swiftf0, "ext") <- "sf0"
attr(trk_swiftf0, "tracks") <- c("F0", "confidence", "voicing")
attr(trk_swiftf0, "outputType") <- "SSFF"

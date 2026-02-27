##' Helper functions for converting av audio data to Python-compatible formats
##'
##' This file contains functions that enable memory-based processing for Python DSP
##' procedures by converting audio data loaded with the av package directly into
##' numpy arrays compatible with libraries like librosa, without writing intermediate
##' files to disk.
##'
##' @name av_python_helpers
##' @keywords internal
NULL

#' Convert av audio data to Python numpy array
#'
#' This function takes audio data read by av::read_audio_bin and converts it to
#' a Python numpy array compatible with librosa and other Python DSP libraries.
#' This enables memory-based processing without the "convert → store → read → DSP" loop.
#'
#' @param audio_data Integer vector from av::read_audio_bin (s32le format)
#' @param sample_rate Sample rate in Hz
#' @param channels Number of audio channels
#'
#' @return A Python numpy array (float64) normalized to \[-1.0, 1.0\] range
#'
#' @examples
#' \dontrun{
#' # Read audio using av
#' info <- av::av_media_info("audio.wav")
#' audio_raw <- av::read_audio_bin("audio.wav",
#'                                  channels = info$audio$channels,
#'                                  sample_rate = info$audio$sample_rate)
#'
#' # Convert to Python numpy array
#' audio_np <- av_to_python_audio(audio_raw,
#'                                 sample_rate = info$audio$sample_rate,
#'                                 channels = info$audio$channels)
#'
#' # Now use with Python DSP (e.g., librosa, pysptk)
#' # The audio_np is already in memory - no file I/O needed!
#' }
av_to_python_audio <- function(audio_data, sample_rate, channels = 1) {

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required but not installed. Please install it with: install.packages('reticulate')")
  }

  # Check if audio_data is valid
  if (length(audio_data) == 0) {
    stop("Empty audio_data provided to av_to_python_audio", call. = FALSE)
  }

  # av::read_audio_bin returns 32-bit signed integers (s32le format)
  # We need to convert to float64 in range [-1.0, 1.0] for librosa

  # Convert from 32-bit int range to float
  # av returns s32le: range is -2147483648 to 2147483647 (2^31)
  audio_float <- as.numeric(audio_data) / 2147483648.0

  # Check conversion succeeded
  if (length(audio_float) == 0) {
    stop("Audio conversion to float failed", call. = FALSE)
  }

  # De-interleave channels if multi-channel (though most Python DSP expects mono)
  if (channels > 1) {
    n_frames <- length(audio_float) / channels
    # For now, take first channel (librosa default behavior)
    # More sophisticated handling could average or let user choose
    audio_float <- audio_float[seq(1, length(audio_float), by = channels)]
  }

  # Convert to Python numpy array
  # Use convert=TRUE to ensure proper conversion when passed through reticulate
  np <- reticulate::import("numpy", convert = TRUE)
  audio_np <- np$array(audio_float, dtype = np$float64)

  # Verify numpy array was created
  if (is.null(audio_np) || length(audio_np) == 0) {
    stop("Failed to create numpy array from audio data", call. = FALSE)
  }

  return(audio_np)
}


#' Load audio file with av and pass to Python as numpy array
#'
#' This function combines av file loading with Python numpy conversion,
#' providing a complete memory-based pipeline for Python DSP processing.
#'
#' @param file_path Path to the audio/video file
#' @param start_time Start time in seconds (default 0)
#' @param end_time End time in seconds (default NULL for end of file)
#' @param target_sample_rate Target sample rate (default NULL to keep original)
#'
#' @return A list with:
#'   - audio_np: Python numpy array (float64) of audio samples
#'   - sample_rate: Sample rate in Hz
#'   - original_file: Original file path
#'
#' @examples
#' \dontrun{
#' # Load audio and get Python-ready numpy array
#' result <- av_load_for_python("myfile.mp4", start_time = 10, end_time = 20)
#'
#' # Pass directly to Python DSP without intermediate files
#' py <- reticulate::import_main()
#' py$audio <- result$audio_np
#' py$fs <- result$sample_rate
#'
#' # Run Python DSP code
#' reticulate::py_run_string("
#'   import pysptk
#'   f0 = pysptk.trk_swipe(audio, fs=fs, hopsize=int(0.005*fs))
#' ")
#' }
av_load_for_python <- function(file_path, start_time = 0, end_time = NULL,
                                target_sample_rate = NULL) {

  if (Sys.getenv("SUPERASSP_DEBUG") != "") {
    message("DEBUG: av_load_for_python called")
    message("DEBUG: file_path = ", file_path)
    message("DEBUG: start_time = ", start_time)
    message("DEBUG: end_time = ", ifelse(is.null(end_time), "NULL", end_time))
  }

  if (!requireNamespace("av", quietly = TRUE)) {
    stop("Package 'av' is required but not installed. Please install it with: install.packages('av')")
  }

  # Read audio info first
  info <- av::av_media_info(file_path)

  if (length(info$audio) == 0) {
    stop("No audio stream found in file: ", file_path)
  }

  audio_info <- info$audio
  original_sample_rate <- audio_info$sample_rate
  channels <- audio_info$channels

  # Use original sample rate if not specified
  if (is.null(target_sample_rate)) {
    target_sample_rate <- original_sample_rate
  }

  # Calculate time window
  duration <- info$duration
  if (is.null(end_time)) {
    end_time <- duration
  }

  # Validate time window
  if (start_time < 0) start_time <- 0
  if (end_time > duration) end_time <- duration
  if (start_time >= end_time) {
    stop("Invalid time window: start_time (", start_time,
         ") >= end_time (", end_time, ")")
  }

  # Read audio data using av
  audio_data <- av::read_audio_bin(file_path,
                                    channels = channels,
                                    start_time = start_time,
                                    end_time = end_time,
                                    sample_rate = target_sample_rate)

  # Check if audio data is empty
  if (length(audio_data) == 0) {
    stop("Failed to read audio data from file: ", file_path,
         "\nTime window: ", start_time, " to ", end_time, " seconds",
         call. = FALSE)
  }

  # Debug: check audio_data before conversion
  if (Sys.getenv("SUPERASSP_DEBUG") != "") {
    message("DEBUG: audio_data length = ", length(audio_data))
    message("DEBUG: channels = ", channels)
  }

  # Convert to Python numpy array
  audio_np <- av_to_python_audio(audio_data, target_sample_rate, channels)

  # Debug: check audio_np after conversion
  if (Sys.getenv("SUPERASSP_DEBUG") != "") {
    message("DEBUG: audio_np length = ", length(audio_np))
    message("DEBUG: audio_np class = ", paste(class(audio_np), collapse=", "))
  }

  return(list(
    audio_np = audio_np,
    sample_rate = target_sample_rate,
    original_file = file_path,
    start_time = start_time,
    end_time = end_time
  ))
}


#' Process media files with Python DSP using memory-based audio loading
#'
#' This is the Python equivalent of processMediaFiles_LoadAndProcess. It loads
#' audio files with av and passes them as numpy arrays to Python DSP functions,
#' eliminating the "convert → store → read → DSP" pattern.
#'
#' @param listOfFiles Character vector of input file paths
#' @param beginTime Numeric vector of begin times (seconds)
#' @param endTime Numeric vector of end times (seconds, 0 = end of file)
#' @param python_function Character name of Python function to call
#' @param ... Additional parameters to pass to the Python function
#'
#' @return List of results from Python function
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # Define a Python DSP function
#' reticulate::py_run_string("
#' def compute_swipe(audio_np, fs, hopsize, fmin, fmax):
#'     import pysptk
#'     f0 = pysptk.trk_swipe(audio_np, fs=fs, hopsize=hopsize, min=fmin, max=fmax)
#'     return f0
#' ")
#'
#' # Process files using memory-based approach
#' results <- processMediaFiles_Python(
#'   listOfFiles = c("file1.mp4", "file2.wav"),
#'   beginTime = c(0, 5),
#'   endTime = c(0, 15),
#'   python_function = "compute_swipe",
#'   hopsize = 220,
#'   fmin = 70,
#'   fmax = 200
#' )
#' }
processMediaFiles_Python <- function(listOfFiles, beginTime, endTime,
                                     python_function, ...) {

  n_files <- length(listOfFiles)

  # Ensure time vectors match file count
  if(length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if(length(endTime) == 1) endTime <- rep(endTime, n_files)

  # Normalize paths
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles))

  # Prepare results storage
  results <- vector("list", n_files)

  # Get Python function
  py <- reticulate::import_main()
  py_func <- py[[python_function]]

  # Process each file using memory-based approach
  for (i in seq_along(listOfFiles)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- if(endTime[i] == 0.0) NULL else endTime[i]

    # Load audio with av → convert to numpy (in memory!)
    audio_data <- av_load_for_python(
      file_path,
      start_time = bt,
      end_time = et,
      target_sample_rate = NULL
    )

    # Call Python function with audio numpy array
    # Pass additional parameters from ...
    results[[i]] <- py_func(
      audio_np = audio_data$audio_np,
      fs = audio_data$sample_rate,
      ...
    )
  }

  return(results)
}

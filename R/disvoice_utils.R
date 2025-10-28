# DisVoice Utility Functions
#
# Helper functions for converting between R/av audio formats and
# Python/numpy arrays for DisVoice processing.

#' Convert av Audio Data to NumPy Array
#'
#' Converts audio data loaded by av::read_audio_bin() from INT32 format
#' to float64 numpy array suitable for DisVoice processing.
#'
#' @param audio_data Audio data from av::read_audio_bin()
#' @param sample_rate Sample rate (extracted from attributes if NULL)
#'
#' @return List with `array` (numpy array) and `sample_rate` (integer)
#' @keywords internal
as_numpy_audio <- function(audio_data, sample_rate = NULL) {
  # Get sample rate from attributes if not provided
  if (is.null(sample_rate)) {
    sample_rate <- attr(audio_data, "sample_rate")
    if (is.null(sample_rate)) {
      stop("Could not determine sample rate from audio data")
    }
  }

  # Convert INT32 to float64 (-1.0 to 1.0 range)
  # av uses signed 32-bit integers: range [-2^31, 2^31-1]
  INT32_MAX <- 2147483647
  audio_float <- as.numeric(audio_data) / INT32_MAX

  # Convert to numpy array
  np <- reticulate::import("numpy", convert = FALSE)
  audio_np <- np$array(audio_float, dtype = "float64")

  list(
    array = audio_np,
    sample_rate = as.integer(sample_rate)
  )
}

#' Load Audio File as NumPy Array
#'
#' Convenience function that loads an audio file using av and
#' converts it to numpy format in one step.
#'
#' @param audio_path Path to audio file
#' @param channels Number of channels (1 = mono, 2 = stereo)
#'
#' @return List with `array` (numpy array) and `sample_rate` (integer)
#' @keywords internal
load_audio_as_numpy <- function(audio_path, channels = 1) {
  if (!file.exists(audio_path)) {
    stop(sprintf("Audio file not found: %s", audio_path))
  }

  # Load with av
  audio_data <- av::read_audio_bin(audio_path, channels = channels)

  # Convert to numpy
  as_numpy_audio(audio_data)
}

#' Load Audio File as Parselmouth Sound Object
#'
#' Convenience function that loads an audio file using av and
#' converts it to a Parselmouth Sound object for DisVoice processing.
#'
#' @param audio_path Path to audio file
#' @param channels Number of channels (1 = mono, 2 = stereo)
#'
#' @return List with `sound` (Parselmouth Sound object) and `sample_rate` (integer)
#' @keywords internal
load_audio_as_sound <- function(audio_path, channels = 1) {
  # Load as numpy first
  audio_np <- load_audio_as_numpy(audio_path, channels = channels)

  # Create Parselmouth Sound object
  py <- reticulate::import("parselmouth", convert = FALSE)
  sound <- py$Sound(audio_np$array, sampling_frequency = audio_np$sample_rate)

  list(
    sound = sound,
    sample_rate = audio_np$sample_rate
  )
}

#' Get Audio Duration
#'
#' Get duration of audio file in seconds using av.
#'
#' @param audio_path Path to audio file
#'
#' @return Duration in seconds
#' @keywords internal
get_audio_duration <- function(audio_path) {
  if (!file.exists(audio_path)) {
    stop(sprintf("Audio file not found: %s", audio_path))
  }

  info <- av::av_media_info(audio_path)
  info$duration
}

#' Convert NumPy Array to R Vector
#'
#' Converts a numpy array returned by DisVoice functions to an R vector.
#' Handles automatic conversion via reticulate.
#'
#' @param np_array NumPy array (Python object)
#'
#' @return R numeric vector
#' @keywords internal
numpy_to_r <- function(np_array) {
  # Check if it's a Python object
  if (inherits(np_array, "python.builtin.object")) {
    # Convert using reticulate
    reticulate::py_to_r(np_array)
  } else {
    # Already an R object
    as.numeric(np_array)
  }
}

#' Extract Voicing Track from TextGrid
#'
#' Converts a Parselmouth TextGrid object with voicing decisions (V/U intervals)
#' into a binary voicing track at regular time intervals.
#'
#' @param textgrid Parselmouth TextGrid object
#' @param times Numeric vector of time points (in seconds)
#'
#' @return Numeric vector: 1 = voiced, 0 = unvoiced
#' @keywords internal
extract_voicing_track <- function(textgrid, times) {
  # Get Parselmouth call function
  call <- reticulate::import("parselmouth.praat", convert = FALSE)$call

  # Get number of intervals in first tier
  n_intervals_py <- call(textgrid, "Get number of intervals", 1L)
  n_intervals <- reticulate::py_to_r(n_intervals_py)

  # Initialize voicing vector
  voicing <- numeric(length(times))

  # Iterate through intervals and mark voiced regions
  for (i in seq_len(n_intervals)) {
    label_py <- call(textgrid, "Get label of interval", 1L, as.integer(i))
    start_time_py <- call(textgrid, "Get start time of interval", 1L, as.integer(i))
    end_time_py <- call(textgrid, "Get end time of interval", 1L, as.integer(i))

    # Convert to R
    label <- reticulate::py_to_r(label_py)
    start_time <- reticulate::py_to_r(start_time_py)
    end_time <- reticulate::py_to_r(end_time_py)

    # Mark as voiced if label is "V"
    if (label == "V" || label == '"V"') {
      # Find times within this interval
      in_interval <- times >= start_time & times <= end_time
      voicing[in_interval] <- 1
    }
  }

  voicing
}

#' Create AsspDataObj from Tracks
#'
#' Creates an AsspDataObj structure from track data, compatible with
#' superassp conventions.
#'
#' @param tracks Named list of numeric vectors (track data)
#' @param times Numeric vector of time points (seconds)
#' @param sample_rate Track sample rate (Hz)
#' @param orig_freq Original audio sample rate (Hz)
#' @param start_time Start time (seconds)
#' @param audio_file Path to source audio file
#' @param track_format Character: track type descriptor
#'
#' @return AsspDataObj
#' @keywords internal
create_assp_data_obj_from_tracks <- function(tracks,
                                             times,
                                             sample_rate,
                                             orig_freq,
                                             start_time,
                                             audio_file,
                                             track_format) {
  # Determine number of records
  n_records <- length(times)

  # Create track data matrix
  track_names <- names(tracks)
  track_matrix <- do.call(cbind, tracks)
  colnames(track_matrix) <- track_names

  # Create AsspDataObj structure
  obj <- structure(
    list(
      trackFormats = rep(track_format, length(tracks)),
      sampleRate = sample_rate,
      origFreq = orig_freq,
      startTime = start_time,
      startRecord = 1,
      endRecord = n_records,
      duration = times[length(times)] - times[1],
      filePath = audio_file,
      tracks = track_matrix
    ),
    class = c("AsspDataObj", "list")
  )

  obj
}

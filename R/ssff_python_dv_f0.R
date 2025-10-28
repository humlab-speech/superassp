# Pitch (F0) Tracking with DisVoice/Parselmouth
#
# Optimized pitch extraction using in-memory Parselmouth processing.
# 1.5-2.5x faster than file-based Praat methods.
#
# File naming convention: ssff_python_dv_*.R
# - ssff: Returns SSFF-compatible data structure
# - python: Python backend via reticulate
# - dv: DisVoice library
# - f0: Fundamental frequency (pitch) method

#' Track Fundamental Frequency (F0) using DisVoice
#'
#' Extracts pitch (F0) contour from audio using DisVoice's optimized
#' Parselmouth implementation. This function uses in-memory processing
#' and is 1.5-2.5x faster than file-based Praat methods.
#'
#' @param audio_path Character string: path to audio file (any format supported by av)
#' @param frame_shift Numeric: frame shift in milliseconds (default: 10ms = 100 Hz)
#' @param min_f0 Numeric: minimum F0 in Hz (default: 75)
#' @param max_f0 Numeric: maximum F0 in Hz (default: 600)
#' @param include_voicing Logical: include voicing decision track (default: TRUE)
#' @param output_format Character: "AsspDataObj", "dataframe", or "list"
#' @param ... Additional arguments (reserved for future use)
#'
#' @return AsspDataObj with F0 track (and optionally voicing track),
#'         or data.frame/list depending on output_format
#'
#' @details
#' This function uses DisVoice's `extract_pitch_and_voicing()` which:
#' - Loads audio with av (in-memory, no temp files)
#' - Processes with Parselmouth Python API
#' - Returns F0 values, time stamps, and TextGrid with voicing decisions
#' - Converts to equal-interval track format
#'
#' The voicing track is binary: 1 = voiced, 0 = unvoiced
#'
#' **Performance**: Approximately 2x faster than file-based Praat methods
#' due to elimination of file I/O overhead.
#'
#' **Backend**: Python (DisVoice + Parselmouth)
#'
#' @references
#' \insertCite{Jadoul2018}{superassp}
#'
#' \insertCite{OrozcoArroyave2018}{superassp}
#'
#' @seealso
#' \code{\link{has_disvoice_support}} to check availability
#' \code{\link{install_disvoice_python}} to install dependencies
#' \code{\link{trk_dv_formants}} for formant tracking
#'
#' @export
#' @examples
#' \dontrun{
#' # Check if DisVoice is available
#' if (has_disvoice_support()) {
#'   # Basic usage
#'   f0_track <- trk_dv_f0("audio.wav")
#'
#'   # Custom parameters
#'   f0_track <- trk_dv_f0(
#'     "audio.wav",
#'     frame_shift = 5,  # 5ms = 200 Hz
#'     min_f0 = 60,
#'     max_f0 = 500,
#'     include_voicing = TRUE
#'   )
#'
#'   # Get as data.frame
#'   f0_df <- trk_dv_f0("audio.wav", output_format = "dataframe")
#' }
#' }
trk_dv_f0 <- function(audio_path,
                      frame_shift = 10,
                      min_f0 = 75,
                      max_f0 = 600,
                      include_voicing = TRUE,
                      output_format = c("AsspDataObj", "dataframe", "list"),
                      ...) {
  # Check DisVoice availability
  if (!has_disvoice_support()) {
    stop(
      "DisVoice Python support not available. ",
      "Install with: install_disvoice_python()"
    )
  }

  # Validate inputs
  if (!file.exists(audio_path)) {
    stop(sprintf("Audio file not found: %s", audio_path))
  }
  output_format <- match.arg(output_format)

  # Load audio with av and convert to Parselmouth Sound
  audio_sound <- load_audio_as_sound(audio_path, channels = 1)

  # Get DisVoice environment
  disvoice <- get_disvoice_env()

  # Calculate signal length for full contour extraction
  len_signal <- length(av::read_audio_bin(audio_path, channels = 1)) /
    audio_sound$sample_rate

  # Extract pitch and voicing using DisVoice
  result <- disvoice$praat_functions$extract_pitch_and_voicing(
    audio_filename = audio_sound$sound,
    time_stepF0 = frame_shift / 1000,  # Convert ms to seconds
    minf0 = as.integer(min_f0),
    maxf0 = as.integer(max_f0),
    maxVUVPeriod = 0.02,
    averageVUVPeriod = 0.01,
    len_signal = len_signal
  )

  # Extract components: (F0_array, time_array, textgrid)
  # Note: Python tuple is 0-indexed
  f0_values <- numpy_to_r(result[[0]])
  times <- numpy_to_r(result[[1]])
  textgrid <- result[[2]]

  # Create tracks list
  tracks <- list(f0 = f0_values)

  # Extract voicing if requested
  if (include_voicing) {
    voicing <- extract_voicing_track(textgrid, times)
    tracks$voicing <- voicing
  }

  # Return in requested format
  if (output_format == "list") {
    return(list(
      tracks = tracks,
      times = times,
      sample_rate = 1 / (frame_shift / 1000),
      audio_file = audio_path
    ))
  }

  if (output_format == "dataframe") {
    df <- data.frame(
      time = times,
      f0 = f0_values
    )
    if (include_voicing) {
      df$voicing <- voicing
    }
    return(df)
  }

  # Return as AsspDataObj (default)
  create_assp_data_obj_from_tracks(
    tracks = tracks,
    times = times,
    sample_rate = 1 / (frame_shift / 1000),
    orig_freq = audio_sound$sample_rate,
    start_time = 0,
    audio_file = audio_path,
    track_format = "F0"
  )
}

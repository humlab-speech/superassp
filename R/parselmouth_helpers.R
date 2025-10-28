#' Convert av audio data to parselmouth Sound object
#'
#' @description
#' Converts audio data loaded by \code{\link[av]{read_audio_bin}} to a
#' parselmouth Sound object for in-memory processing. This eliminates the
#' need for temporary WAV files.
#'
#' @param audio_data Integer vector of audio samples from \code{av::read_audio_bin}.
#'   Expected to be INT32 format in range [-2147483648, 2147483647].
#' @param sample_rate Integer sample rate in Hz. If NULL, extracted from
#'   audio_data attributes.
#' @param normalize Logical. Normalize audio to [-1, 1] range (default: TRUE).
#'   Required for parselmouth compatibility.
#'
#' @return A parselmouth.Sound object that can be used with all parselmouth
#'   functions (pitch analysis, formant extraction, etc.)
#'
#' @details
#' This function performs the conversion:
#' \enumerate{
#'   \item Extracts sample rate from av audio data attributes
#'   \item Converts INT32 samples to float64 normalized to [-1, 1]
#'   \item Creates numpy array via reticulate
#'   \item Constructs parselmouth.Sound object in memory
#' }
#'
#' The resulting Sound object can be used with any parselmouth/Praat function:
#' - \code{sound$to_pitch()} - Pitch analysis
#' - \code{sound$to_formant()} - Formant analysis
#' - \code{sound$to_intensity()} - Intensity analysis
#' - \code{sound$to_spectrogram()} - Spectrogram
#' - Any Praat command via \code{parselmouth.praat.call()}
#'
#' @examples
#' \dontrun{
#' # Load audio using av
#' audio_data <- av::read_audio_bin("speech.wav", channels = 1)
#'
#' # Convert to parselmouth Sound (in memory, no file I/O)
#' sound <- av_to_parselmouth_sound(audio_data)
#'
#' # Use with parselmouth functions
#' pitch <- sound$to_pitch(time_step = 0.01, pitch_floor = 75, pitch_ceiling = 600)
#' formants <- sound$to_formant_burg()
#' intensity <- sound$to_intensity()
#'
#' # Extract values
#' pm <- reticulate::import("parselmouth")
#' mean_pitch <- pm$praat$call(pitch, "Get mean", 0, 0, "Hertz")
#' }
#'
#' @seealso
#' \code{\link{av_load_for_parselmouth}} for complete workflow with time windowing
#'
#' @export
av_to_parselmouth_sound <- function(audio_data,
                                    sample_rate = NULL,
                                    normalize = TRUE) {

  # Check if parselmouth is available
  if (!reticulate::py_module_available("parselmouth")) {
    cli::cli_abort(c(
      "x" = "parselmouth module not available",
      "i" = "Install with: reticulate::py_install('praat-parselmouth')"
    ))
  }

  # Extract sample rate from attributes if not provided
  if (is.null(sample_rate)) {
    sample_rate <- attr(audio_data, "sample_rate")
    if (is.null(sample_rate)) {
      cli::cli_abort("sample_rate must be provided or present in audio_data attributes")
    }
  }

  # Import modules
  pm <- reticulate::import("parselmouth", delay_load = FALSE)
  np <- reticulate::import("numpy", convert = FALSE, delay_load = FALSE)

  # Convert INT32 to float
  # av::read_audio_bin returns INT32 in range [-2147483648, 2147483647]
  INT32_MAX <- 2147483647

  if (normalize) {
    # Normalize to [-1, 1] range (required for parselmouth)
    audio_float <- as.numeric(audio_data) / INT32_MAX
  } else {
    audio_float <- as.numeric(audio_data)
  }

  # Convert R vector to numpy array
  # parselmouth expects float64
  audio_np <- np$array(audio_float, dtype = "float64")

  # Create parselmouth Sound object
  sound <- pm$Sound(audio_np, sampling_frequency = as.integer(sample_rate))

  return(sound)
}


#' Load audio for parselmouth processing
#'
#' @description
#' Complete workflow to load audio using av package and convert to parselmouth
#' Sound object. Handles time windowing and format conversion automatically.
#'
#' @param file_path Character path to audio file (any format supported by av)
#' @param start_time Numeric start time in seconds (default: NULL = file start)
#' @param end_time Numeric end time in seconds (default: NULL = file end)
#' @param channels Integer number of channels (default: 1 = mono)
#' @param target_sample_rate Integer target sample rate in Hz (default: NULL = original)
#'
#' @return A parselmouth.Sound object ready for analysis
#'
#' @details
#' This function combines \code{av::read_audio_bin} and
#' \code{av_to_parselmouth_sound} into a single convenient call.
#'
#' Advantages over file-based approach:
#' \itemize{
#'   \item No temporary files created (pure in-memory processing)
#'   \item Time windowing handled by av (efficient)
#'   \item Supports all media formats (WAV, MP3, MP4, video, etc.)
#'   \item Automatic resampling if needed
#'   \item Cleaner code, no temp file cleanup required
#' }
#'
#' @examples
#' \dontrun{
#' # Load entire file
#' sound <- av_load_for_parselmouth("speech.wav")
#'
#' # Load with time windowing
#' sound <- av_load_for_parselmouth("speech.wav",
#'                                   start_time = 1.0,
#'                                   end_time = 3.0)
#'
#' # Load and resample
#' sound <- av_load_for_parselmouth("speech.mp3",
#'                                   target_sample_rate = 16000)
#'
#' # Use with parselmouth
#' pitch <- sound$to_pitch()
#' formants <- sound$to_formant_burg()
#' }
#'
#' @seealso
#' \code{\link{av_to_parselmouth_sound}} for the conversion step only
#'
#' @export
av_load_for_parselmouth <- function(file_path,
                                    start_time = NULL,
                                    end_time = NULL,
                                    channels = 1,
                                    target_sample_rate = NULL) {

  # Load audio using av package
  audio_data <- av::read_audio_bin(
    audio = file_path,
    start_time = start_time,
    end_time = end_time,
    channels = channels,
    sample_rate = target_sample_rate
  )

  # Convert to parselmouth Sound
  sound <- av_to_parselmouth_sound(audio_data)

  return(sound)
}


#' Check if parselmouth is available
#'
#' @description
#' Helper function to check if the parselmouth Python module is installed
#' and can be imported.
#'
#' @return Logical TRUE if parselmouth is available, FALSE otherwise
#'
#' @examples
#' \dontrun{
#' if (parselmouth_available()) {
#'   sound <- av_load_for_parselmouth("audio.wav")
#' } else {
#'   message("Install parselmouth: reticulate::py_install('praat-parselmouth')")
#' }
#' }
#'
#' @export
parselmouth_available <- function() {
  reticulate::py_module_available("parselmouth")
}

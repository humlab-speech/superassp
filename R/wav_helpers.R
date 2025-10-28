#' Write WAV file from raw audio data
#'
#' Internal helper function to write a minimal WAV file from raw audio samples.
#' This eliminates the need for external dependencies like tuneR.
#'
#' @param audio_data Integer vector of audio samples (from av::read_audio_bin)
#' @param sample_rate Integer sample rate in Hz
#' @param filename Character path to output WAV file
#' @param bit_depth Integer bits per sample (default: 16)
#'
#' @details
#' Creates a PCM WAV file with RIFF header. Supports:
#' - Mono audio only
#' - 16-bit or 32-bit PCM
#' - Standard WAV format compatible with all audio software
#'
#' The audio_data should be raw integer samples as returned by av::read_audio_bin
#' (INT32 format, -2147483648 to 2147483647 range).
#'
#' @return Invisibly returns TRUE on success
#'
#' @keywords internal
#' @noRd
write_wav_file <- function(audio_data, sample_rate, filename, bit_depth = 16) {

  # Validate inputs
  if (!is.numeric(audio_data)) {
    stop("audio_data must be numeric", call. = FALSE)
  }
  if (length(audio_data) == 0) {
    stop("audio_data is empty", call. = FALSE)
  }
  if (sample_rate <= 0) {
    stop("sample_rate must be positive", call. = FALSE)
  }
  if (!bit_depth %in% c(16, 32)) {
    stop("bit_depth must be 16 or 32", call. = FALSE)
  }

  # Convert audio data based on bit depth
  if (bit_depth == 16) {
    # Convert INT32 to INT16
    # av::read_audio_bin returns INT32 in range [-2147483648, 2147483647]
    # Scale to INT16 range [-32768, 32767]
    audio_scaled <- as.numeric(audio_data) / 65536  # 2^16
    audio_int <- as.integer(pmax(-32768, pmin(32767, audio_scaled)))
    bytes_per_sample <- 2
  } else {
    # Keep as INT32
    audio_int <- as.integer(audio_data)
    bytes_per_sample <- 4
  }

  # Calculate sizes
  num_samples <- length(audio_int)
  num_channels <- 1  # Mono only
  byte_rate <- sample_rate * num_channels * bytes_per_sample
  block_align <- num_channels * bytes_per_sample
  data_size <- num_samples * bytes_per_sample

  # Open binary connection
  con <- file(filename, "wb")
  on.exit(close(con), add = TRUE)

  # Write RIFF header
  writeBin(charToRaw("RIFF"), con)
  writeBin(as.integer(36 + data_size), con, size = 4, endian = "little")
  writeBin(charToRaw("WAVE"), con)

  # Write fmt chunk
  writeBin(charToRaw("fmt "), con)
  writeBin(as.integer(16), con, size = 4, endian = "little")  # Chunk size
  writeBin(as.integer(1), con, size = 2, endian = "little")   # Audio format (1 = PCM)
  writeBin(as.integer(num_channels), con, size = 2, endian = "little")
  writeBin(as.integer(sample_rate), con, size = 4, endian = "little")
  writeBin(as.integer(byte_rate), con, size = 4, endian = "little")
  writeBin(as.integer(block_align), con, size = 2, endian = "little")
  writeBin(as.integer(bit_depth), con, size = 2, endian = "little")

  # Write data chunk
  writeBin(charToRaw("data"), con)
  writeBin(as.integer(data_size), con, size = 4, endian = "little")

  # Write audio samples
  writeBin(audio_int, con, size = bytes_per_sample, endian = "little")

  invisible(TRUE)
}


#' Load audio for Python processing
#'
#' Internal helper to load audio using av package and prepare for Python functions.
#' Handles time windowing and returns both the audio data and a temporary WAV file path.
#'
#' @param file_path Character path to audio file
#' @param start_time Numeric start time in seconds (default: NULL = file start)
#' @param end_time Numeric end time in seconds (default: NULL = file end)
#' @param channels Integer number of channels (default: 1 = mono)
#' @param target_sample_rate Integer target sample rate (default: NULL = original)
#'
#' @return List with:
#'   \item{audio_data}{Raw audio samples}
#'   \item{sample_rate}{Sample rate in Hz}
#'   \item{temp_file}{Path to temporary WAV file (must be unlinked after use)}
#'
#' @keywords internal
#' @noRd
av_load_for_python <- function(file_path,
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

  # Get sample rate
  sample_rate <- attr(audio_data, "sample_rate")

  # Create temporary WAV file
  temp_wav <- tempfile(fileext = ".wav")

  # Write WAV file using our internal helper
  write_wav_file(
    audio_data = audio_data,
    sample_rate = sample_rate,
    filename = temp_wav,
    bit_depth = 16
  )

  return(list(
    audio_data = audio_data,
    sample_rate = sample_rate,
    temp_file = temp_wav
  ))
}

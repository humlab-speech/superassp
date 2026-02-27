#' AVAudio S7 Class
#'
#' An S7 class for representing audio data with metadata, wrapping the output
#' of `prep_recode()` or `av::read_audio_bin()`.
#'
#' @description
#' AVAudio is an S7 class that encapsulates audio sample data along with
#' essential metadata (sample rate, channels). It provides a modern object-oriented
#' interface for audio processing functions in superassp.
#'
#' @details
#' The AVAudio class stores:
#' - `samples`: Integer vector of audio samples (s32le format - 32-bit signed integers)
#' - `sample_rate`: Sample rate in Hz (integer)
#' - `channels`: Number of audio channels (integer)
#' - `file_path`: Optional source file path (character, NULL if generated)
#'
#' AVAudio objects can be created from:
#' 1. File paths via `prep_recode()` or `read_avaudio()`
#' 2. Raw audio data via `as_avaudio()`
#' 3. Existing av::read_audio_bin() output via `as_avaudio()`
#'
#' @examples
#' \dontrun{
#' # Create from file
#' audio <- read_avaudio("speech.wav")
#'
#' # Create from prep_recode output
#' audio_data <- prep_recode("speech.wav", format = "wav")
#' audio <- as_avaudio(audio_data)
#'
#' # Access properties
#' audio$sample_rate
#' audio$channels
#' length(audio$samples)
#'
#' # Use in DSP functions
#' f0 <- trk_rapt(audio)
#' features <- lst_voice_sauce(audio)
#' }
#'
#' @name AVAudio-class
#' @export
AVAudio <- S7::new_class(
  "AVAudio",
  properties = list(
    samples = S7::class_integer,
    sample_rate = S7::class_integer,
    channels = S7::class_integer,
    file_path = S7::new_property(
      class = S7::class_character,
      default = NA_character_
    )
  ),
  validator = function(self) {
    # Validate sample_rate
    if (length(self@sample_rate) != 1) {
      return("sample_rate must be a single integer")
    }
    if (self@sample_rate <= 0) {
      return("sample_rate must be positive")
    }

    # Validate channels
    if (length(self@channels) != 1) {
      return("channels must be a single integer")
    }
    if (self@channels <= 0) {
      return("channels must be positive (typically 1 or 2)")
    }

    # Validate samples length is consistent with channels
    if (length(self@samples) %% self@channels != 0) {
      return("samples length must be a multiple of channels")
    }

    # Validate file_path if present
    if (length(self@file_path) != 1) {
      return("file_path must be a single character or NA")
    }
  }
)

#' Create AVAudio Object from File
#'
#' Read an audio file and create an AVAudio object.
#'
#' @param file_path Character; path to audio file
#' @param format Character; output format (default: "wav")
#' @param sample_rate Integer; target sample rate in Hz (default: NULL, keep original)
#' @param channels Integer; number of channels (default: NULL, keep original)
#' @param start_time Numeric; start time in seconds (default: NULL)
#' @param end_time Numeric; end time in seconds (default: NULL)
#' @param ... Additional arguments passed to `prep_recode()`
#'
#' @return AVAudio object
#'
#' @examples
#' \dontrun{
#' # Read entire file
#' audio <- read_avaudio("speech.wav")
#'
#' # Read with time windowing
#' audio <- read_avaudio("speech.wav", start_time = 1.0, end_time = 3.0)
#'
#' # Read with resampling
#' audio <- read_avaudio("speech.wav", sample_rate = 16000)
#' }
#'
read_avaudio <- function(file_path,
                         format = "wav",
                         sample_rate = NULL,
                         channels = NULL,
                         start_time = NULL,
                         end_time = NULL,
                         ...) {

  if (!file.exists(file_path)) {
    stop("File not found: ", file_path, call. = FALSE)
  }

  # Use prep_recode to load audio
  audio_data <- prep_recode(
    listOfFiles = file_path,
    format = format,
    sample_rate = sample_rate,
    channels = channels,
    start_time = start_time,
    end_time = end_time,
    verbose = FALSE,
    ...
  )

  if (is.null(audio_data)) {
    stop("Failed to read audio file: ", file_path, call. = FALSE)
  }

  # Convert to AVAudio
  as_avaudio(audio_data, file_path = file_path)
}

#' Convert to AVAudio Object
#'
#' Convert audio data (from av::read_audio_bin or prep_recode) to AVAudio object.
#'
#' @param x Integer vector with audio samples (must have channels and sample_rate attributes)
#' @param file_path Character; optional source file path
#'
#' @return AVAudio object
#'
#' @examples
#' \dontrun{
#' # From av::read_audio_bin
#' audio_data <- av::read_audio_bin("speech.wav")
#' audio <- as_avaudio(audio_data)
#'
#' # From prep_recode
#' audio_data <- prep_recode("speech.wav", format = "wav")
#' audio <- as_avaudio(audio_data, file_path = "speech.wav")
#' }
#'
as_avaudio <- function(x, file_path = NA_character_) {

  if (!is.integer(x)) {
    stop("x must be an integer vector (audio samples)", call. = FALSE)
  }

  # Extract metadata
  sample_rate <- attr(x, "sample_rate", exact = TRUE)
  channels <- attr(x, "channels", exact = TRUE)

  if (is.null(sample_rate)) {
    stop("x must have 'sample_rate' attribute", call. = FALSE)
  }
  if (is.null(channels)) {
    stop("x must have 'channels' attribute", call. = FALSE)
  }

  # Ensure integer types
  sample_rate <- as.integer(sample_rate)
  channels <- as.integer(channels)

  # Create AVAudio object
  AVAudio(
    samples = x,
    sample_rate = sample_rate,
    channels = channels,
    file_path = as.character(file_path)
  )
}

#' Check if Object is AVAudio
#'
#' @param x Object to check
#' @return Logical; TRUE if x is an AVAudio object
is_avaudio <- function(x) {
  S7::S7_inherits(x, AVAudio)
}

#' Convert AVAudio to av::read_audio_bin Format
#'
#' Convert AVAudio object back to the format returned by av::read_audio_bin(),
#' which is an integer vector with channels and sample_rate attributes.
#'
#' @param audio AVAudio object
#' @return Integer vector with attributes (channels, sample_rate)
#'
#' @examples
#' \dontrun{
#' audio <- read_avaudio("speech.wav")
#' audio_vec <- avaudio_to_av(audio)
#'
#' # audio_vec is now compatible with av::read_audio_bin output
#' }
#'
avaudio_to_av <- function(audio) {
  if (!is_avaudio(audio)) {
    stop("audio must be an AVAudio object", call. = FALSE)
  }

  # Extract samples
  samples <- audio@samples

  # Add attributes
  attr(samples, "sample_rate") <- audio@sample_rate
  attr(samples, "channels") <- audio@channels

  samples
}

#' Convert AVAudio to Temporary WAV File
#'
#' Write AVAudio object to a temporary WAV file and return the path.
#' Useful for passing to external DSP functions that require file paths.
#'
#' @param audio AVAudio object
#' @param verbose Logical; show messages (default: FALSE)
#'
#' @return Character; path to temporary WAV file
#'
#' @details
#' The temporary file is created with `tempfile()` and will be automatically
#' deleted when the R session ends, unless explicitly deleted earlier with `unlink()`.
#'
#' @examples
#' \dontrun{
#' audio <- read_avaudio("speech.wav")
#' temp_path <- avaudio_to_tempfile(audio)
#'
#' # Use temp file with external tools
#' # ... processing ...
#'
#' # Clean up (optional - happens automatically at session end)
#' unlink(temp_path)
#' }
#'
avaudio_to_tempfile <- function(audio, verbose = FALSE) {
  if (!is_avaudio(audio)) {
    stop("audio must be an AVAudio object", call. = FALSE)
  }

  # If audio has source file, just return it if no modifications
  if (!is.na(audio@file_path) && file.exists(audio@file_path)) {
    # For now, always create temp file to ensure consistency
    # Future: optimize by returning source file when no modifications
  }

  # Create temporary file
  temp_file <- tempfile(fileext = ".wav")

  # Write WAV file using base R WAV writing
  # This is a simple WAV writer for PCM data
  .write_wav_file(temp_file, audio@samples, audio@sample_rate, audio@channels)

  temp_file
}

#' Write WAV File
#'
#' Internal function to write audio samples to WAV file format.
#'
#' @param filename Character; output file path
#' @param samples Integer vector; audio samples (s32le format)
#' @param sample_rate Integer; sample rate in Hz
#' @param channels Integer; number of channels
#' @keywords internal
.write_wav_file <- function(filename, samples, sample_rate, channels) {
  # Convert s32le samples to INT16 for WAV
  # Normalize from 32-bit to 16-bit range
  samples_16 <- as.integer(samples / 65536)

  # Clip to INT16 range
  samples_16[samples_16 > 32767] <- 32767L
  samples_16[samples_16 < -32768] <- -32768L

  # Calculate sizes
  n_samples <- length(samples_16)
  byte_rate <- sample_rate * channels * 2  # 2 bytes per sample (INT16)
  block_align <- channels * 2
  data_size <- n_samples * 2
  file_size <- 36 + data_size

  # Open file for binary writing
  con <- file(filename, "wb")
  on.exit(close(con))

  # RIFF header
  writeBin(charToRaw("RIFF"), con)
  writeBin(as.integer(file_size), con, size = 4, endian = "little")
  writeBin(charToRaw("WAVE"), con)

  # fmt chunk
  writeBin(charToRaw("fmt "), con)
  writeBin(as.integer(16), con, size = 4, endian = "little")  # Chunk size
  writeBin(as.integer(1), con, size = 2, endian = "little")   # PCM format
  writeBin(as.integer(channels), con, size = 2, endian = "little")
  writeBin(as.integer(sample_rate), con, size = 4, endian = "little")
  writeBin(as.integer(byte_rate), con, size = 4, endian = "little")
  writeBin(as.integer(block_align), con, size = 2, endian = "little")
  writeBin(as.integer(16), con, size = 2, endian = "little")  # Bits per sample

  # data chunk
  writeBin(charToRaw("data"), con)
  writeBin(as.integer(data_size), con, size = 4, endian = "little")
  writeBin(samples_16, con, size = 2, endian = "little")

  invisible(filename)
}

#' Print Method for AVAudio
#'
#' @param x AVAudio object
#' @param ... Additional arguments (unused)
#' @name print.AVAudio
#' @keywords internal
S7::method(print, AVAudio) <- function(x, ...) {
  n_samples <- length(x@samples)
  duration <- n_samples / x@channels / x@sample_rate

  cat("<AVAudio>\n")
  cat("  Sample rate:", x@sample_rate, "Hz\n")
  cat("  Channels:", x@channels, "\n")
  cat("  Samples:", n_samples, "\n")
  cat("  Duration:", round(duration, 3), "seconds\n")

  if (!is.na(x@file_path)) {
    cat("  Source:", x@file_path, "\n")
  }

  invisible(x)
}

#' Summary Method for AVAudio
#'
#' @param object AVAudio object
#' @param ... Additional arguments (unused)
#' @name summary.AVAudio
#' @keywords internal
S7::method(summary, AVAudio) <- function(object, ...) {
  n_samples <- length(object@samples)
  duration <- n_samples / object@channels / object@sample_rate

  # Calculate statistics per channel
  if (object@channels == 1) {
    samples_matrix <- matrix(object@samples, ncol = 1)
  } else {
    samples_matrix <- matrix(object@samples, ncol = object@channels, byrow = TRUE)
  }

  cat("<AVAudio Summary>\n")
  cat("  Sample rate:", object@sample_rate, "Hz\n")
  cat("  Channels:", object@channels, "\n")
  cat("  Total samples:", n_samples, "\n")
  cat("  Duration:", round(duration, 3), "seconds\n")

  if (!is.na(object@file_path)) {
    cat("  Source:", object@file_path, "\n")
  }

  cat("\n  Sample statistics (per channel):\n")
  for (ch in 1:object@channels) {
    channel_samples <- samples_matrix[, ch]
    cat("    Channel", ch, ":\n")
    cat("      Min:", min(channel_samples), "\n")
    cat("      Max:", max(channel_samples), "\n")
    cat("      Mean:", round(mean(channel_samples), 2), "\n")
    cat("      SD:", round(sd(channel_samples), 2), "\n")
  }

  invisible(object)
}

#' Package-level Constants for Signal Processing
#'
#' @description
#' This file defines package-level constants used throughout superassp for
#' signal processing, audio handling, and file format specifications.
#'
#' @name constants
#' @keywords internal
NULL

# Audio Processing Constants ====

#' @rdname constants
#' @keywords internal
PHONET_SAMPLE_RATE <- 16000L  # Phonet requires 16 kHz

#' @rdname constants
#' @keywords internal
DEFAULT_CHANNELS <- 1L  # Mono audio

#' @rdname constants
#' @keywords internal
SACC_SAMPLE_RATE <- 16000L  # SAcC pitch tracker requires 16 kHz

#' @rdname constants
#' @keywords internal
BROUHAHA_SAMPLE_RATE <- 16000L  # Brouhaha-VAD requires 16 kHz

# Time Conversion Constants ====

#' @rdname constants
#' @keywords internal
MS_TO_SECONDS <- 1000  # Milliseconds per second

#' @rdname constants
#' @keywords internal
US_TO_SECONDS <- 1000000  # Microseconds per second

# JSTF Format Constants ====

#' @rdname constants
#' @keywords internal
JSTF_FORMAT_VERSION <- "1.0"

#' @rdname constants
#' @keywords internal
JSTF_DEFAULT_COMPRESSION <- TRUE

# Sample Rate Conversion ====

#' @rdname constants
#' @keywords internal
STANDARD_SAMPLE_RATES <- c(8000L, 11025L, 16000L, 22050L, 32000L, 44100L, 48000L, 96000L)

# AVQI Constants ====

#' @rdname constants
#' @keywords internal
AVQI_MIN_SV_DURATION_MS <- 1000  # Minimum sustained vowel duration (1 second)

#' @rdname constants
#' @keywords internal
AVQI_MIN_CS_DURATION_MS <- 1000  # Minimum continuous speech duration (1 second)

# Time Conversion Helper Functions ====

#' Convert milliseconds to seconds
#'
#' @param ms Numeric value in milliseconds
#' @return Numeric value in seconds
#' @keywords internal
#' @noRd
ms_to_sec <- function(ms) {
  ms / MS_TO_SECONDS
}

#' Convert seconds to milliseconds
#'
#' @param sec Numeric value in seconds
#' @return Numeric value in milliseconds
#' @keywords internal
#' @noRd
sec_to_ms <- function(sec) {
  sec * MS_TO_SECONDS
}

#' Convert microseconds to seconds
#'
#' @param us Numeric value in microseconds
#' @return Numeric value in seconds
#' @keywords internal
#' @noRd
us_to_sec <- function(us) {
  us / US_TO_SECONDS
}

#' Convert seconds to microseconds
#'
#' @param sec Numeric value in seconds
#' @return Numeric value in microseconds
#' @keywords internal
#' @noRd
sec_to_us <- function(sec) {
  sec * US_TO_SECONDS
}

# Sample Rate Validation ====

#' Check if sample rate is standard
#'
#' @param sr Sample rate in Hz
#' @return Logical indicating if sample rate is standard
#' @keywords internal
#' @noRd
is_standard_sample_rate <- function(sr) {
  sr %in% STANDARD_SAMPLE_RATES
}

#' Find nearest standard sample rate
#'
#' @param sr Sample rate in Hz
#' @return Nearest standard sample rate
#' @keywords internal
#' @noRd
nearest_standard_sample_rate <- function(sr) {
  STANDARD_SAMPLE_RATES[which.min(abs(STANDARD_SAMPLE_RATES - sr))]
}

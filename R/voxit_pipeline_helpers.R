#' Internal WORLD vocoder pipeline
#'
#' Orchestrate Harvest, CheapTrick, D4C on pre-loaded raw PCM signal.
#' Computes normalized linear power from CheapTrick spectrogram.
#'
#' @param wave Numeric vector; raw PCM samples in [-1, 1]
#' @param fs Integer; original audio sample rate (Hz)
#' @param frame_period Numeric; frame shift in milliseconds (default 5ms)
#'
#' @return List with three elements:
#'   - f0_parameter: list(temporal_positions, f0, vuv)
#'   - spectrum_parameter: list(spectrogram, temporal_positions, fs, linPower)
#'   - source_parameter: list(aperiodicity, temporal_positions, fs)
#'
#' @keywords internal
.voxit_world_pipeline <- function(wave, fs, frame_period = 5.0) {
  f0_param <- harvest_r(wave, fs, f0_floor = 71.0, f0_ceil = 800.0,
                        frame_period = frame_period)
  spec_param <- cheap_trick_r(wave, fs, f0_param$temporal_positions,
                              f0_param$f0)
  src_param <- d4c_r(wave, fs, f0_param$temporal_positions, f0_param$f0)

  lin_power <- colSums(spec_param$spectrogram) / nrow(spec_param$spectrogram)
  lin_power <- lin_power / max(lin_power)

  list(
    f0_parameter = f0_param,
    spectrum_parameter = modifyList(spec_param, list(linPower = lin_power)),
    source_parameter = src_param
  )
}

#' Internal: Harvest F0 on raw PCM for voxit pipeline
#'
#' @param wave Numeric vector; PCM samples
#' @param fs Integer; sample rate
#' @param f0_floor Numeric; lower F0 bound (Hz)
#' @param f0_ceil Numeric; upper F0 bound (Hz)
#' @param frame_period Numeric; frame shift in ms
#' @return list(f0, temporal_positions, fs)
#' @keywords internal
harvest_r <- function(wave, fs, f0_floor = 71.0, f0_ceil = 800.0,
                      frame_period = 5.0) {
  audio_obj <- list(audio = matrix(as.double(wave), ncol = 1))
  attr(audio_obj, "sampleRate") <- as.integer(fs)
  attr(audio_obj, "origFreq")   <- as.integer(fs)
  class(audio_obj) <- "AsspDataObj"

  res <- harvest_cpp(
    audio_obj = audio_obj,
    minF = f0_floor,
    maxF = f0_ceil,
    windowShift = frame_period,
    voicing_threshold = 0.1,
    verbose = FALSE
  )
  list(f0 = as.numeric(res$f0),
       temporal_positions = as.numeric(res$times),
       fs = fs)
}

#' Internal: CheapTrick spectral envelope on raw PCM for voxit pipeline
#'
#' @param wave Numeric vector; PCM samples
#' @param fs Integer; sample rate
#' @param temporal_positions Numeric vector of frame times
#' @param f0 Numeric vector of F0 values
#' @param q1 Regularization parameter
#' @param f0_floor Lower F0 bound for FFT size
#' @return list(spectrogram, temporal_positions, fs, fft_size)
#' @keywords internal
cheap_trick_r <- function(wave, fs, temporal_positions, f0,
                          q1 = -0.15, f0_floor = 71.0) {
  audio_obj <- list(audio = matrix(as.double(wave), ncol = 1))
  attr(audio_obj, "sampleRate") <- as.integer(fs)
  attr(audio_obj, "origFreq")   <- as.integer(fs)
  class(audio_obj) <- "AsspDataObj"

  res <- cheap_trick_cpp(
    audio_obj = audio_obj,
    f0 = f0,
    temporal_positions = temporal_positions,
    q1 = q1,
    f0_floor = f0_floor,
    verbose = FALSE
  )
  list(spectrogram = res$spectrogram,
       temporal_positions = temporal_positions,
       fs = fs,
       fft_size = res$fft_size)
}

#' Internal: D4C aperiodicity on raw PCM for voxit pipeline
#'
#' @param wave Numeric vector; PCM samples
#' @param fs Integer; sample rate
#' @param temporal_positions Numeric vector of frame times
#' @param f0 Numeric vector of F0 values
#' @param threshold D4C VUV threshold
#' @return list(aperiodicity, temporal_positions, fs)
#' @keywords internal
d4c_r <- function(wave, fs, temporal_positions, f0, threshold = 0.85) {
  audio_obj <- list(audio = matrix(as.double(wave), ncol = 1))
  attr(audio_obj, "sampleRate") <- as.integer(fs)
  attr(audio_obj, "origFreq")   <- as.integer(fs)
  class(audio_obj) <- "AsspDataObj"

  res <- d4c_cpp(
    audio_obj = audio_obj,
    f0 = f0,
    temporal_positions = temporal_positions,
    threshold = threshold,
    verbose = FALSE
  )
  list(aperiodicity = res$aperiodicity,
       temporal_positions = temporal_positions,
       fs = fs)
}

#' Check if Voxit WORLD functionality is available
#'
#' Returns TRUE if the WORLD vocoder (Harvest, CheapTrick, D4C) is compiled and available.
#'
#' @return Logical; always TRUE after compilation
#'
#' @export
has_voxit_support <- function() {
  TRUE
}

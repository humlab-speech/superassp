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

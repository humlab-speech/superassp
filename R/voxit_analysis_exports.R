#' Voxit analysis statistics
#'
#' Compute basic F0 and intensity statistics for Voxit analysis.
#'
#' @name voxit-analysis-stats
NULL

#' F0 statistics
#'
#' Compute geometric mean F0, range (octaves), and voicing percentage.
#'
#' @param f0 Numeric vector of F0 values (Hz; 0 = unvoiced)
#' @param vuv Numeric vector of voicing confidence (0–1)
#'
#' @return Numeric vector: `c(geometric_mean_hz, range_octaves, voicing_percent)`
#' @keywords internal
#' @noRd
compute_f0_stats_simple_cpp <- function(f0, vuv) {
  .Call(`_superassp_compute_f0_stats_simple_cpp`, f0, vuv)
}

#' Pause frame count
#'
#' Count frames below intensity threshold and unvoiced.
#'
#' @param lin_power Linear power vector (0–1)
#' @param vuv Voicing vector (0–1)
#' @param thresh_db dB threshold below 90th percentile (default 10)
#'
#' @return Integer pause frame count
#' @keywords internal
#' @noRd
compute_pause_count_cpp <- function(lin_power, vuv, thresh_db = 10.0) {
  .Call(`_superassp_compute_pause_count_cpp`, lin_power, vuv, thresh_db)
}

#' Intensity mean
#'
#' Convert linear power to dB and compute mean.
#'
#' @param lin_power Linear power vector (0–1)
#'
#' @return Mean intensity in dB
#' @keywords internal
#' @noRd
compute_intensity_mean_cpp <- function(lin_power) {
  .Call(`_superassp_compute_intensity_mean_cpp`, lin_power)
}

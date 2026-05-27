#' DSP utility functions from Voxit
#'
#' Fast C++ implementations of Savitzky-Golay filtering, run detection,
#' and histogram binning.
#'
#' @name voxit-dsp-utils
NULL

#' Savitzky-Golay filter
#'
#' Apply a fixed Savitzky-Golay smoothing filter (degree 2, window 7).
#'
#' @param x Numeric vector to filter
#' @param order Polynomial degree (fixed at 2)
#' @param window Window size (fixed at 7)
#'
#' @return Smoothed numeric vector
#' @keywords internal
#' @noRd
sgolay_filter_cpp <- function(x, order = 2, window = 7) {
  .Call(`_superassp_sgolay_filter_cpp`, x, order, window)
}

#' Find contiguous runs
#'
#' Detect contiguous runs of a target value in a vector.
#'
#' @param x Integer vector
#' @param target_value Value to find runs of (default 1)
#'
#' @return List with elements `starts`, `stops`, `lengths` (1-indexed)
#' @keywords internal
#' @noRd
contiguous_runs_cpp <- function(x, target_value = 1L) {
  .Call(`_superassp_contiguous_runs_cpp`, as.integer(x), as.integer(target_value))
}

#' Fixed-width histogram
#'
#' Compute histogram bin counts for a fixed interval.
#'
#' @param x Numeric vector
#' @param nbins Number of bins
#' @param lo Lower bound of range
#' @param hi Upper bound of range
#'
#' @return Integer vector of bin counts
#' @keywords internal
#' @noRd
histcounts_cpp <- function(x, nbins, lo, hi) {
  .Call(`_superassp_histcounts_cpp`, x, as.integer(nbins), lo, hi)
}

#' Quantile (Type 7 linear interpolation)
#'
#' Compute quantiles matching R's Type 7 method.
#'
#' @param x Numeric vector
#' @param p Quantile (0–1)
#'
#' @return Numeric quantile value
#' @keywords internal
#' @noRd
quantile_cpp <- function(x, p) {
  .Call(`_superassp_quantile_cpp`, x, p)
}

#' Linear interpolation
#'
#' Linear interpolation at arbitrary points.
#'
#' @param xp X positions of known points
#' @param fp Y values at known points (same length as xp)
#' @param xi X positions to interpolate at
#'
#' @return Interpolated Y values at xi
#' @keywords internal
#' @noRd
interp1_linear_cpp <- function(xp, fp, xi) {
  .Call(`_superassp_interp1_linear_cpp`, xp, fp, xi)
}

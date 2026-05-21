#' Peak Slope measurement
#' Computes the PeakSlope parameter on a fixed frame basis using multi-level
#' wavelet decomposition. Translates \code{get_peakSlope.m} (John Kane, 2011).
#' @note Now uses the original Daless wavelet bank recovered from earlier git
#'   history (commit bb9b314). Pass \code{backend = "r"} to use the legacy db4
#'   approximation kept for back-compat.
#' @param s Speech signal (numeric vector)
#' @param fs Sampling frequency (Hz)
#' @return Numeric vector of PeakSlope values measured every 10 ms.
#' @references Kane, J., Gobl, C. (2011). Identifying regions of non-modal
#'   phonation using features of the wavelet transform. Proc. Interspeech.
#' @param backend "cpp" (default, Daless wavelets in Rcpp) or "r" (db4 fallback).
#' @export
.vat_peak_slope <- function(s, fs, backend = c("cpp", "r")) {
  backend <- match.arg(backend)
  s <- as.numeric(s)
  if (backend == "cpp") {
    return(vat_peakslope_cpp(s, fs))
  }

  frame_len_ms  <- 40   # 40 ms -> covers one period down to 25 Hz F0
  frame_shft_ms <- 10
  frame_len     <- round(frame_len_ms  / 1000 * fs)
  frame_shft    <- round(frame_shft_ms / 1000 * fs)

  n_frames   <- floor((length(s) - frame_len) / frame_shft)
  peak_slope <- numeric(n_frames)

  # Multi-level DWT using db4 (approximation of Daless wavelets)
  # 7 bands: i = 0:6 corresponding to fs/2, fs/4, fs/8, ..., fs/128
  n_levels <- 6L  # yields 7 sub-bands: 1 approx + 6 detail

  if (!requireNamespace("wavelets", quietly = TRUE)) {
    # Fallback: use simple octave-band filtering via successive half-band LP filters
    y <- wavelet_octave_bands_r(s, n_levels)
  } else {
    y <- wavelet_bands_wavelets_pkg(s, n_levels)
  }

  # Measure peak slope per frame
  n_bands <- length(y)
  for (m in seq_len(n_frames)) {
    start  <- (m - 1L) * frame_shft + 1L
    finish <- start + frame_len - 1L
    if (finish > length(s)) break

    maxima <- vapply(y, function(band) max(abs(band[start:finish])), numeric(1))
    # Reverse order (low freq first) and convert to log10
    maxima <- log10(rev(maxima) + 1e-10)
    t      <- seq_along(maxima)
    # Linear regression slope
    p_coef <- stats::lm.fit(cbind(1, t), maxima)$coefficients
    peak_slope[m] <- p_coef[2]
  }

  peak_slope[is.nan(peak_slope)] <- 0
  peak_slope
}

# ── Internal: wavelet band decomposition ─────────────────────────────────────

#' DWT band decomposition using 'wavelets' package
#' @keywords internal
wavelet_bands_wavelets_pkg <- function(s, n_levels) {
  wt  <- wavelets::dwt(s, filter = "d4", n.levels = n_levels,
                       boundary = "periodic")
  # Collect approximation + detail coefficients reconstructed to signal length
  n_samp <- length(s)
  bands  <- vector("list", n_levels + 1L)

  for (lvl in seq_len(n_levels)) {
    # Reconstruct single detail level
    wt_tmp <- wt
    for (l2 in seq_len(n_levels)) {
      if (l2 != lvl) wt_tmp@W[[l2]] <- numeric(length(wt@W[[l2]]))
    }
    wt_tmp@V[[n_levels]] <- numeric(length(wt@V[[n_levels]]))
    rec <- wavelets::idwt(wt_tmp)
    if (length(rec) > n_samp) rec <- rec[seq_len(n_samp)]
    else if (length(rec) < n_samp) rec <- c(rec, numeric(n_samp - length(rec)))
    bands[[lvl]] <- rec
  }
  # Approximation at coarsest level
  wt_tmp <- wt
  for (lvl in seq_len(n_levels)) wt_tmp@W[[lvl]] <- numeric(length(wt@W[[lvl]]))
  rec <- wavelets::idwt(wt_tmp)
  if (length(rec) > n_samp) rec <- rec[seq_len(n_samp)]
  else if (length(rec) < n_samp) rec <- c(rec, numeric(n_samp - length(rec)))
  bands[[n_levels + 1L]] <- rec

  bands
}

#' Octave-band decomposition fallback (successive half-band LP filters)
#' @keywords internal
wavelet_octave_bands_r <- function(s, n_levels) {
  bands  <- vector("list", n_levels + 1L)
  n_samp <- length(s)
  x_cur  <- s

  for (lvl in seq_len(n_levels)) {
    # Half-band lowpass
    fc   <- 0.5  # normalised cutoff (Nyquist of current signal)
    filt <- signal::butter(4L, 0.45)  # slight rolloff
    lo   <- signal::filtfilt(filt$b, filt$a, x_cur)
    # Detail = original minus lowpass
    hi   <- x_cur - lo
    # Ensure same length
    if (length(hi) > n_samp) hi <- hi[seq_len(n_samp)]
    else if (length(hi) < n_samp) hi <- c(hi, numeric(n_samp - length(hi)))
    bands[[lvl]] <- hi
    x_cur <- lo[seq(1, length(lo), by = 2)]  # downsample for next level
  }
  # Approximation: upsample back to original length
  approx_up <- stats::approx(seq_along(x_cur), x_cur, xout = seq_len(n_samp), rule = 2)$y
  bands[[n_levels + 1L]] <- approx_up

  bands
}

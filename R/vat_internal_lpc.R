#' Compute LP residual (overlap-add)
#' Frame-by-frame LPC analysis with overlap-add residual synthesis.
#' Translates \code{GetLPCresidual.m} (Thomas Drugman / John Kane).
#' @param wave Speech signal (numeric vector)
#' @param window_len Window length in samples (typ. 25 ms * fs)
#' @param shift Window shift in samples (typ. 5 ms * fs)
#' @param order LPC order (typ. fs/1000 + 2)
#' @return Normalized LP residual signal (same length as \code{wave})
#' @keywords internal
#' @noRd
.vat_lpc_residual <- function(wave, window_len, shift, order) {
  wave <- as.numeric(wave)
  vat_lpc_residual_cpp(wave, as.integer(round(window_len)),
                   as.integer(round(shift)),
                   as.integer(order))
}

#' Autocorrelation LPC analysis
#' Computes LPC coefficients and residual energy for a windowed frame using
#' the autocorrelation method (Levinson-Durbin).
#' @param s Signal frame (numeric vector), should already be windowed
#' @param p LPC order
#' @return List with \code{ar} (LPC coefficients, first = 1) and \code{e} (residual energy)
#' @keywords internal
lpcauto_frame <- function(s, p) {
  s <- as.numeric(s)
  n <- length(s)
  if (n <= p) {
    return(list(ar = c(1, rep(0, p)), e = sum(s^2)))
  }
  # Autocorrelation
  r <- stats::acf(s, lag.max = p, type = "covariance", plot = FALSE)$acf[, 1, 1]
  if (r[1] == 0) return(list(ar = c(1, rep(0, p)), e = 0))

  # Levinson-Durbin
  a <- numeric(p)
  e <- r[1]
  for (m in seq_len(p)) {
    k <- -(r[m + 1] + sum(a[seq_len(m - 1)] * r[m:2])) / e
    if (m > 1) {
      a_new <- a[seq_len(m - 1)] + k * rev(a[seq_len(m - 1)])
      a[seq_len(m - 1)] <- a_new
    }
    a[m] <- k
    e <- e * (1 - k^2)
    if (e <= 0) { e <- 1e-10; break }
  }
  list(ar = c(1, a), e = e)
}

#' Resonator applied to LP residual (RCVD)
#' Zero-phase resonator filter at F0mean, used to generate post-processing
#' reference signal for creak GCI removal. Translates \code{RCVD_reson_GCI.m}.
#' @param res LP residual signal
#' @param fs Sampling frequency (Hz)
#' @param f0_mean Mean F0 (Hz)
#' @return Normalized resonator output (same length as \code{res})
#' @keywords internal
#' @noRd
.vat_rcvd_reson <- function(res, fs, f0_mean) {
  Phi <- 2 * pi * f0_mean / fs
  Rho <- 0.9
  b_filt <- c(1, 0, 0)
  a_filt <- c(1, -2 * Rho * cos(Phi), Rho^2)
  rep_out <- signal::filtfilt(b_filt, a_filt, res)
  rep_out / max(abs(rep_out))
}

#' Zero-phase high-pass Butterworth filter
#' @param x Input signal
#' @param fs Sampling frequency (Hz)
#' @param f_p Pass-band frequency (Hz)
#' @param f_s Stop-band frequency (Hz)
#' @return Filtered signal
#' @keywords internal
zero_phase_hp_filt <- function(x, fs, f_p, f_s) {
  Wp <- f_p / (fs / 2)
  Ws <- f_s / (fs / 2)
  # Butterworth order estimation
  ord <- signal::buttord(Wp, Ws, Rp = 0.5, Rs = 40)
  filt <- signal::butter(ord$n, ord$Wc, type = "high")
  signal::filtfilt(filt$b, filt$a, x)
}

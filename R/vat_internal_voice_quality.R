#' Voice quality parameters: NAQ, QOQ, H1H2, HRF
#' Computes Normalized Amplitude Quotient (NAQ), Quasi-Open Quotient (QOQ),
#' H1-H2 spectral measure, and Harmonic Richness Factor (HRF) from an
#' estimated glottal source signal. Translates \code{get_NAQ_QOQ_H1H2.m}.
#' @param glot Glottal flow derivative estimate (from \code{.vat_iaif()})
#' @param fs Sampling frequency (Hz)
#' @param GCI Glottal closure instants (integer sample positions)
#' @return List with numeric vectors \code{NAQ}, \code{QOQ}, \code{H1H2}, \code{HRF},
#'   each of length \code{length(GCI)}.
#' @references \insertCite{Kane2013VoiceQuality}{superassp}
#' @param backend "cpp" (default, Rcpp) or "r" (legacy R orchestration).
#' @keywords internal
#' @noRd
.vat_voice_quality <- function(glot, fs, GCI, backend = c("cpp", "r")) {
  backend <- match.arg(backend)
  glot <- as.numeric(glot)
  GCI  <- as.integer(GCI)
  if (backend == "cpp") {
    return(vat_naq_qoq_h1h2_cpp(glot, fs, GCI))
  }
  n    <- length(GCI)

  F0min    <- 20
  F0max    <- 500
  qoq_level <- 0.5
  T0_num    <- 3L
  glot_shift <- round(0.5 / 1000 * fs)

  NAQ  <- numeric(n)
  QOQ  <- numeric(n)
  H1H2 <- numeric(n)
  HRF  <- numeric(n)

  glot_int <- integrat_r(glot, fs)

  for (i in seq_len(n)) {
    if (i == 1L) {
      start <- 1L
      stop  <- GCI[i]
      if (n > 1) T0 <- GCI[2] - GCI[1] else next
    } else {
      start <- GCI[i - 1]
      stop  <- GCI[i]
      T0    <- GCI[i] - GCI[i - 1]
    }

    if (T0 == 0) next
    F0 <- fs / T0
    if (!is.finite(F0) || F0 <= F0min || F0 >= F0max) next

    # Drift-compensated glottal integral over [start, stop]
    gi_ends <- glot_int[c(start, stop)]
    seg_len <- stop - start + 1L
    if (stop > start && seg_len >= 2L) {
      line <- stats::approx(c(1L, seg_len), gi_ends,
                            xout = seq_len(seg_len))$y
    } else {
      line <- 0
    }
    glot_int_cur  <- glot_int[start:stop]
    glot_int_comp <- glot_int_cur - line

    stop2 <- min(stop + glot_shift, length(glot))
    glot_cur <- glot[start:stop2]

    # H1-H2 frame
    half_T0  <- round(T0 * T0_num / 2)
    f_start  <- max(1L, GCI[i] - half_T0)
    f_stop   <- min(length(glot), GCI[i] + half_T0)
    f_frame  <- glot[f_start:f_stop]
    f_win    <- f_frame * signal::hamming(length(f_frame))
    # FFT at fs points to get Hz-indexed spectrum
    nfft     <- as.integer(fs)
    f_spec   <- 20 * log10(abs(stats::fft(c(f_win, numeric(nfft - length(f_win))))) + 1e-10)
    f_spec   <- f_spec[seq_len(nfft %/% 2)]  # one-sided

    # Find harmonic peaks
    h_peaks <- pracma::findpeaks(f_spec, minpeakdistance = round(F0 / 2))
    if (!is.null(h_peaks) && nrow(h_peaks) >= 5) {
      h_idx <- h_peaks[, 2]
      h_amp <- h_peaks[, 1]
      f0_idx  <- which.min(abs(h_idx - F0))
      f02_idx <- which.min(abs(h_idx - F0 * 2))
      H1H2[i] <- h_amp[f0_idx] - h_amp[f02_idx]
      HRF[i]  <- sum(h_amp[f0_idx]) / h_amp[f0_idx[1]]
    }

    # NAQ and QOQ
    d_peak  <- max(abs(glot_cur))
    f_ac    <- max(glot_int_comp)
    max_idx <- which.max(glot_int_comp)
    Amid    <- f_ac * qoq_level

    t12 <- find_amid_t(glot_int_comp, Amid, max_idx)
    T1  <- t12[1]
    T2  <- t12[2]

    if (d_peak != 0) NAQ[i] <- (f_ac / d_peak) * F0
    QOQ[i] <- (T2 - T1) / (fs / F0)
  }

  QOQ[QOQ < 0 | QOQ > 1] <- 0

  list(NAQ = NAQ, QOQ = QOQ, H1H2 = H1H2, HRF = HRF)
}

# ── Internal helpers ─────────────────────────────────────────────────────────

#' Numerical integrator (Euler, Ts = 1/fs)
#' @keywords internal
integrat_r <- function(x, fs) {
  Ts <- 1 / fs
  cumsum(x) * Ts
}

#' Find open-phase start and end times for QOQ
#' @keywords internal
find_amid_t <- function(glot_adj, Amid, Tz) {
  T1 <- 0L
  T2 <- 0L
  if (Tz == 0) return(c(T1, T2))
  n <- Tz
  while (n > 3 && glot_adj[n] > Amid) n <- n - 1L
  T1 <- n
  n  <- Tz
  lim <- length(glot_adj) - 2L
  while (n < lim && glot_adj[n] > Amid) n <- n + 1L
  T2 <- n
  c(T1, T2)
}

#' SE-VQ GCI Detection
#' Detects glottal closure instants (GCIs) using the SE-VQ algorithm
#' (Kane & Gobl 2013), an improved SEDREAMS method optimised for non-modal
#' voice qualities. Translates \code{SE_VQ.m} / \code{SE_VQ_varF0.m}.
#' @param x Speech signal (numeric vector, samples)
#' @param fs Sampling frequency (Hz)
#' @param f0 F0 contour at sample rate (Hz). If NULL, computed via \code{.vat_pitch()}.
#' @param VUV Voiced/unvoiced binary decisions (same length as x). If NULL,
#'   computed with f0.
#' @param creak Optional creaky voice binary decisions (same length as x) for
#'   post-processing.
#' @param var_f0 If TRUE, use variable-F0 window length (SE_VQ_varF0).
#' @param backend "cpp" (default, full Rcpp pipeline) or "r" (legacy R orchestration).
#' @return List with: GCI (integer), rep, res, MBS, F0mean, F0max.
#' @references \insertCite{Kane2013GCI}{superassp}
#' @keywords internal
#' @noRd
.vat_se_vq <- function(x, fs, f0 = NULL, VUV = NULL, creak = NULL,
                       var_f0 = FALSE, backend = c("cpp", "r")) {
  backend <- match.arg(backend)
  x <- as.numeric(x)
  F0min <- 20
  F0max <- 500

  if (is.null(f0) || is.null(VUV)) {
    pitch_res <- .vat_pitch(x, fs, F0min, F0max)
    f0  <- pitch_res$f0
    VUV <- pitch_res$VUV
    if (length(f0) != length(x)) {
      ts <- seq(0, by = 0.010, length.out = length(f0))
      tx <- (seq_along(x) - 1) / fs
      f0  <- stats::approx(ts, f0,  tx, rule = 2)$y
      VUV <- as.integer(stats::approx(ts, VUV, tx, method = "constant", rule = 2)$y >= 0.5)
    }
  }

  if (backend == "cpp") {
    creak_arg <- if (is.null(creak)) NULL else as.integer(creak)
    return(vat_se_vq_cpp(x, fs, f0, as.integer(VUV), creak_arg, var_f0))
  }

  # ---- legacy R-orchestration path ----

  voiced_f0 <- f0[VUV == 1 & f0 > F0min & f0 < F0max]
  F0mean <- if (length(voiced_f0) > 0) stats::median(voiced_f0) else 100
  voiced_f0_for_max <- f0[VUV == 1]
  if (length(voiced_f0_for_max) >= 13) {
    F0max_local <- max(signal::medfilt1(voiced_f0_for_max, 13))
  } else {
    F0max_local <- F0max
  }
  if (F0mean < 70) {
    message("Utterance likely to contain creak")
    F0mean <- 80
  }
  T0mean <- fs / F0mean

  # Parameters
  win_len  <- round(25 / 1000 * fs)   # 25 ms window
  win_shft <- round(5  / 1000 * fs)   # 5 ms shift
  lpc_ord  <- round(fs / 1000) + 2
  Ncand    <- 5L
  trans_wgt   <- 1.0
  rel_amp_wgt <- 0.3
  rep_num         <- 2L
  remove_thresh   <- 0.4
  search_reg      <- round(1.3 / 1000 * fs)

  # LP residual
  res <- .vat_lpc_residual(x, win_len, win_shft, lpc_ord)

  # Resonator output
  rep <- .vat_rcvd_reson(res, fs, F0mean)

  # Mean-based signal
  MBS <- get_mbs(x, fs, T0mean)

  # GCI search intervals from MBS
  interval <- get_mbs_gci_intervals(MBS, fs, T0mean, F0max_local)

  # Find Ncand residual peaks per interval
  peaks <- search_res_interval_peaks(res, interval, Ncand)
  GCI_N       <- peaks$GCI       # N x Ncand
  GCI_relAmp  <- peaks$rel_amp   # N x Ncand

  # Dynamic programming (Rcpp)
  GCI <- vat_reson_dynprog_cpp(
    t(GCI_relAmp),   # ncands x nframe
    t(GCI_N),        # ncands x nframe
    F0mean, x, fs, trans_wgt, rel_amp_wgt
  )

  # Creak post-processing
  if (!is.null(creak) && length(creak) == length(x)) {
    message("Doing post-processing in detected creaky voice regions")
    GCI <- gci_creak_postproc(GCI, creak, search_reg, rep, remove_thresh, rep_num)
  }

  list(GCI = GCI, rep = rep, res = res, MBS = MBS)
}

# ── Internal helpers ────────────────────────────────────────────────────────

#' Mean-based signal (MBS)
#' @keywords internal
get_mbs <- function(x, fs, T0mean) {
  n <- length(x)
  MBS <- numeric(n)

  halfL <- round(1.6 * T0mean[1] / 2)
  Step  <- 2^3  # 8

  idx <- seq(halfL + 1, n - halfL, by = Step)
  for (m in idx) {
    halfL_m <- if (length(T0mean) > 1) round(1.7 * T0mean[m] / 2) else round(1.7 * T0mean[1] / 2)
    start <- round(m - halfL_m)
    stop  <- round(m + halfL_m)
    if (stop > n) break
    if (start > 0) {
      vec <- x[start:stop] * signal::blackman(stop - start + 1)
      MBS[m] <- mean(vec)
    }
  }

  # Interpolate over zero-valued gaps
  nz <- which(MBS != 0)
  if (length(nz) < 2) return(rep(0, n))
  MBS <- stats::approx(nz, MBS[nz], xout = seq_len(n), rule = 2)$y
  MBS[is.na(MBS)] <- 0

  # High-pass filter at 70 Hz
  MBS <- zero_phase_hp_filt(MBS, fs, 70, 10)

  # Normalize
  mx <- max(MBS)
  if (mx != 0) MBS <- MBS / mx

  # Smooth (7-point moving average)
  MBS <- as.numeric(stats::filter(MBS, rep(1/7, 7), sides = 2))
  MBS[is.na(MBS)] <- 0

  MBS
}

#' GCI search intervals from MBS
#' @keywords internal
get_mbs_gci_intervals <- function(MBS, fs, T0mean, F0max) {
  F0max_2  <- F0max * 2
  T0max    <- round(fs / F0max_2)
  neg_peaks <- pracma::findpeaks(-MBS, minpeakdistance = T0max)

  if (is.null(neg_peaks) || nrow(neg_peaks) == 0) return(matrix(0, 0, 2))

  idx <- neg_peaks[, 2]   # peak locations
  N   <- length(idx)

  search_rate      <- 0.28
  search_left_rate <- 0.01
  n_samp           <- length(MBS)
  interval         <- matrix(0, N, 2)

  for (n in seq_len(N)) {
    if (length(T0mean) > 1) {
      t0 <- T0mean[min(idx[n], length(T0mean))]
    } else {
      t0 <- T0mean[1]
    }
    start <- idx[n] - round(t0 * search_left_rate)
    stop  <- idx[n] + round(t0 * search_rate)
    if (start < 1)       start <- 1L
    if (stop > n_samp && start < n_samp) {
      stop <- n_samp
    } else if (stop > n_samp) {
      break
    }
    interval[n, 1] <- start
    interval[n, 2] <- stop
  }
  interval
}

#' Find top-N residual peaks in each MBS interval
#' @keywords internal
search_res_interval_peaks <- function(res, interval, Ncand) {
  N   <- nrow(interval)
  GCI     <- matrix(NA_real_, N, Ncand)
  rel_amp <- matrix(0,        N, Ncand)

  for (n in seq_len(N)) {
    start <- interval[n, 1]
    stop  <- interval[n, 2]
    if (stop <= start) next
    seg <- res[start:stop]
    len <- length(seg)
    if (len < Ncand) {
      best_idx <- which.max(seg) + start - 1L
      GCI[n, ]     <- best_idx
      rel_amp[n, ] <- 0
    } else {
      ord <- order(seg, decreasing = TRUE)
      top_idx <- ord[seq_len(Ncand)] + start - 1L
      top_amp <- seg[ord[seq_len(Ncand)]]
      GCI[n, ]     <- top_idx
      rel_amp[n, ] <- 1 - (top_amp / top_amp[1])
    }
  }
  list(GCI = GCI, rel_amp = rel_amp)
}

#' Remove false GCIs in creaky regions
#' @keywords internal
gci_creak_postproc <- function(GCI, creak, search_reg, rep, remove_thresh, rep_num) {
  GCI <- as.integer(GCI)
  valid <- GCI >= 1 & GCI <= length(creak)
  GCI_creak <- GCI[valid & creak[GCI] == 1]
  GCI       <- GCI[!(valid & creak[GCI] == 1)]

  for (m in seq_len(rep_num)) {
    n <- 2L
    while (n < length(GCI_creak)) {
      pos <- GCI_creak[n]
      lo <- max(1L, pos - round(search_reg))
      hi <- min(length(rep), pos + round(search_reg))
      cur_rep_max <- abs(min(rep[lo:hi]))
      nb_mean <- mean(c(abs(rep[GCI_creak[n - 1]]),
                        abs(rep[GCI_creak[n + 1]])))
      if (nb_mean * remove_thresh > cur_rep_max) {
        GCI_creak[n] <- NA_integer_
        n <- n + 2L
      } else {
        n <- n + 1L
      }
    }
    GCI_creak <- GCI_creak[!is.na(GCI_creak)]
  }
  sort(unique(c(GCI, GCI_creak)))
}

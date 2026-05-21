#' Creaky voice feature extraction
#' Extracts acoustic features for creaky voice detection. These are the
#' Kane-Drugman features plus Ishi et al. (2008) features. Translates the
#' feature extraction part of \code{get_ALL_creak_features.m}.
#' @note The MATLAB ANN classifier weights (SystemNet_creak.mat) can be
#'   converted via \code{tools/convert_ann_weights.R} to
#'   \code{inst/extdata/creak_ann.rds} and consumed by \code{.vat_creak_detect()}.
#' @param x Speech signal (numeric vector)
#' @param fs Sampling frequency (Hz)
#' @param frame_len_ms Frame length in ms (default 25)
#' @param frame_shift_ms Frame shift in ms (default 10)
#' @return Matrix of features (n_frames x n_features) with static, delta,
#'   and delta-delta coefficients. Columns: H2H1, res_p, ZCR, F0, F0mean,
#'   enerN, pow_std, creak_F0 (+ deltas).
#' @references Kane, J., Drugman, T., Gobl, C. (2013). Improved automatic
#'   detection of creak. Computer Speech and Language, 27(4), 1028-1047.
#' @export
.vat_creak_features <- function(x, fs, frame_len_ms = 25, frame_shift_ms = 10) {
  x <- as.numeric(x)
  frame_len  <- round(frame_len_ms  / 1000 * fs)
  frame_shft <- round(frame_shift_ms / 1000 * fs)

  n_samples <- length(x)
  n_frames  <- floor((n_samples - frame_len) / frame_shft) + 1L

  # LP residual for residual-based features
  lpc_ord <- round(fs / 1000) + 2
  res <- .vat_lpc_residual(x, frame_len, frame_shft, lpc_ord)

  # Resonator output (F0mean estimate for creak resonator)
  F0min <- 20; F0max <- 500
  # Simple autocorrelation-based F0 estimate for feature extraction
  F0est <- estimate_f0_simple(x, fs, F0min, F0max)

  feat_mat <- matrix(0, n_frames, 8)
  colnames(feat_mat) <- c("H2H1", "res_p", "ZCR", "F0", "F0mean",
                           "enerN", "pow_std", "creak_F0")

  for (m in seq_len(n_frames)) {
    start <- (m - 1L) * frame_shft + 1L
    stop  <- start + frame_len - 1L
    if (stop > n_samples) break

    frame  <- x[start:stop]
    res_fr <- res[start:stop]

    # Energy normalised
    enerN    <- sum(frame^2) / frame_len

    # Power std
    pow_std  <- stats::sd(frame^2)

    # ZCR
    zcr      <- sum(abs(diff(sign(frame)))) / (2 * frame_len)

    # H2H1: spectral balance at H2 / H1
    win    <- frame * signal::hamming(frame_len)
    spec   <- abs(stats::fft(win))
    spec_db <- 20 * log10(spec[seq_len(frame_len %/% 2)] + 1e-10)
    F0_fr   <- F0est[min(start, length(F0est))]
    if (F0_fr > F0min && F0_fr < F0max) {
      h1_idx <- round(F0_fr / (fs / frame_len))
      h2_idx <- round(2 * F0_fr / (fs / frame_len))
      h1_idx <- max(1L, min(h1_idx, length(spec_db)))
      h2_idx <- max(1L, min(h2_idx, length(spec_db)))
      H2H1   <- spec_db[h2_idx] - spec_db[h1_idx]
    } else {
      H2H1   <- 0
    }

    # Residual periodicity
    res_p  <- max(abs(res_fr)) / (mean(abs(res_fr)) + 1e-10)

    feat_mat[m, ] <- c(H2H1, res_p, zcr, F0_fr, F0est[1], enerN, pow_std, 0)
  }

  # Deltas and delta-deltas
  feat_d  <- compute_delta(feat_mat)
  feat_dd <- compute_delta(feat_d)

  cbind(feat_mat, feat_d, feat_dd)
}

#' Creaky voice detection
#' Runs feature extraction + the Kane-Drugman ANN classifier (if weights are
#' bundled in \code{inst/extdata/creak_ann.rds}) to produce posterior
#' probabilities and binary decisions per 10 ms frame.
#' @param x Speech signal (numeric vector)
#' @param fs Sampling frequency (Hz)
#' @param threshold ANN decision threshold (default 0.3, matches MATLAB)
#' @return List with: \code{posterior} (probabilities per frame),
#'   \code{decision} (binary), \code{time} (frame center times in s).
#' @param backend "cpp" (default, bit-faithful Kane-Drugman + Ishi port)
#'   or "r" (legacy lighter R-side approximation).
#' @export
.vat_creak_detect <- function(x, fs, threshold = 0.3, backend = c("cpp", "r")) {
  backend <- match.arg(backend)
  ann_path <- system.file("extdata", "creak_ann.rds", package = "voiceanalysis")
  if (ann_path == "") {
    stop("ANN weights not bundled. Run tools/convert_ann_weights.R to install.")
  }
  ann <- readRDS(ann_path)
  if (backend == "cpp") {
    return(vat_creak_detect_cpp(as.numeric(x), fs,
                             ann$IW, ann$b_h, ann$LW, ann$b_o,
                             ann$mini, ann$maxi, threshold))
  }
  # Legacy R path
  feats <- .vat_creak_features(x, fs)
  if (ncol(feats) < length(ann$mini)) {
    stop(sprintf("Feature dim mismatch: %d extracted vs %d expected by ANN",
                 ncol(feats), length(ann$mini)))
  }
  X <- t(feats[, seq_len(length(ann$mini)), drop = FALSE])
  y <- vat_ann_forward_cpp(X, ann$IW, ann$b_h, ann$LW, ann$b_o, ann$mini, ann$maxi)
  list(
    posterior = y,
    decision  = as.integer(y > threshold),
    time      = seq_along(y) * 0.010
  )
}

# ── Internal helpers ─────────────────────────────────────────────────────────

#' Simple autocorrelation F0 estimate (per-sample via frame interpolation)
#' @keywords internal
estimate_f0_simple <- function(x, fs, F0min, F0max) {
  frame_len  <- round(25 / 1000 * fs)
  frame_shft <- round(10 / 1000 * fs)
  n  <- length(x)
  nf <- floor((n - frame_len) / frame_shft) + 1L

  lag_min <- round(fs / F0max)
  lag_max <- round(fs / F0min)

  f0_frames <- numeric(nf)
  times     <- numeric(nf)
  for (m in seq_len(nf)) {
    start <- (m - 1L) * frame_shft + 1L
    stop  <- start + frame_len - 1L
    if (stop > n) break
    frame <- x[start:stop] * signal::hanning(frame_len)
    acf   <- stats::acf(frame, lag.max = lag_max, plot = FALSE)$acf[, 1, 1]
    lags  <- seq(lag_min, lag_max)
    if (length(lags) == 0 || lag_max > length(acf) - 1) {
      f0_frames[m] <- 0
    } else {
      peak_lag    <- lags[which.max(acf[lags + 1])]
      f0_frames[m] <- fs / peak_lag
    }
    times[m] <- (start + stop) / 2
  }

  # Interpolate to per-sample
  if (sum(times > 0) < 2) return(rep(0, n))
  stats::approx(times[times > 0], f0_frames[times > 0],
                xout = seq_len(n), rule = 2)$y
}

#' Compute delta (first derivative) of feature matrix
#' @keywords internal
compute_delta <- function(feat) {
  n  <- nrow(feat)
  nc <- ncol(feat)
  d  <- matrix(0, n, nc)
  for (i in seq_len(n)) {
    il <- max(1L, i - 1L)
    ir <- min(n, i + 1L)
    d[i, ] <- (feat[ir, ] - feat[il, ]) / (ir - il)
  }
  d
}

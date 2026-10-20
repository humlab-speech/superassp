# Pins the vectorised / SIMD helpers of trk_covarep_vad_drugman() and
# lst_polarity() against verbatim copies of the original scalar R
# implementations they replaced (equal within rounding).

# ---- original reference implementations ---------------------------------

ref_xcorr_self <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  out <- numeric(2 * n - 1L)
  lags <- -(n - 1L):(n - 1L)
  for (i in seq_along(lags)) {
    lag <- lags[i]
    if (lag < 0) {
      out[i] <- sum(x[1:(n + lag)] * x[(1 - lag):n])
    } else {
      out[i] <- sum(x[(1 + lag):n] * x[1:(n - lag)])
    }
  }
  out
}

ref_medfilt1 <- function(x, k = 11L) {
  x <- as.numeric(x)
  half <- floor(k / 2)
  out <- numeric(length(x))
  for (i in seq_along(x)) {
    lo <- max(1L, i - half)
    hi <- min(length(x), i + half)
    out[i] <- stats::median(x[lo:hi])
  }
  out
}

ref_filtering_and_derivatives <- function(feat) {
  nfeats <- nrow(feat)
  nframes <- ncol(feat)
  for (k in seq_len(nfeats)) feat[k, ] <- ref_medfilt1(feat[k, ], 11)
  mat_feat <- matrix(0, nrow = nframes, ncol = 3 * nfeats)
  mat_feat[, 1:nfeats] <- t(feat)
  n_around <- 10
  for (k in seq_len(nfeats)) {
    for (l in (n_around + 1L):nframes) {
      vec <- feat[k, (l - n_around):l]
      mat_feat[l, k + nfeats] <- sum(vec[length(vec)] - vec[1:(length(vec) - 1L)])
    }
  }
  for (k in seq_len(nfeats)) {
    for (l in (n_around + 1L):nframes) {
      vec <- mat_feat[(l - n_around):l, k + nfeats]
      mat_feat[l, k + 2 * nfeats] <- sum(vec[length(vec)] - vec[1:(length(vec) - 1L)])
    }
  }
  mat_feat
}

ref_mfcc <- function(wave, fs) {
  win_l <- round(30 / 1000 * fs)
  hopsize <- round(10 / 1000 * fs)
  win <- as.numeric(av::hanning(win_l))
  n_mfcc <- 13; n_mels <- 40; n_fft <- 512
  start <- 1L; stop <- win_l
  mfcc_mat <- matrix(0, nrow = n_mfcc, ncol = 0)
  while (stop <= length(wave)) {
    seg <- wave[start:stop] * win
    spec <- Mod(stats::fft(seg, inverse = FALSE))[1:(n_fft / 2 + 1)]
    spec <- spec^2
    mel_spec <- numeric(n_mels)
    for (m in 1:n_mels) {
      idx <- round(m / (n_mels + 1) * n_fft / 2)
      if (idx >= 1 && idx <= length(spec)) {
        mel_spec[m] <- log(max(spec[idx], .Machine$double.eps))
      }
    }
    mfcc <- mel_spec[1:n_mfcc]
    mfcc_mat <- cbind(mfcc_mat, mfcc)
    start <- start + hopsize; stop <- stop + hopsize
  }
  mfcc_mat
}

ref_hps <- function(seg, fs) {
  xspec <- Mod(stats::fft(seg, inverse = FALSE)[1:1024])
  xspec <- 20 * log10(pmax(xspec, .Machine$double.eps))
  lmin_hz <- round(62.5 / fs * 2048)
  lmax_hz <- round(500 / fs * 2048)
  vec_tmp <- numeric(lmax_hz - lmin_hz + 1L)
  for (kk in lmin_hz:lmax_hz) {
    s <- 0
    for (hrm in 1:8) {
      idx <- hrm * kk
      if (idx <= length(xspec)) s <- s + xspec[idx]
    }
    vec_tmp[kk - lmin_hz + 1L] <- s
  }
  max(vec_tmp)
}

ref_polarity_lpc_autocorr <- function(x, order) {
  n <- length(x)
  r <- numeric(order + 1L)
  for (k in 0:order) {
    idx <- 1:(n - k)
    if (length(idx) > 0) r[k + 1L] <- sum(x[idx] * x[(k + 1):n])
  }
  r
}

# ---- tests ----------------------------------------------------------------

test_that(".vad_acf_lags / .vad_acf_fft equal the lags >= 0 of the original xcorr", {
  set.seed(42)
  for (n in c(1L, 5L, 50L, 480L)) {
    x <- rnorm(n)
    full <- ref_xcorr_self(x)
    for (ml in unique(pmin(c(0L, 3L, 255L), n - 1L))) {
      expect_identical(superassp:::.vad_acf_lags(x, ml), full[n:(n + ml)])
      expect_equal(drop(superassp:::.vad_acf_fft(matrix(x), ml, stats::nextn(n + ml + 1L))),
                   full[n:(n + ml)])
    }
  }
  w <- as.numeric(av::hanning(480))
  m <- cbind(rnorm(480) * w, sin(2 * pi * 200 * (0:479) / 16000) * w, 1e-6 * rnorm(480) * w)
  got <- superassp:::.vad_acf_fft(m, 255L, stats::nextn(480 + 256))
  for (j in seq_len(ncol(m))) expect_equal(got[, j], ref_xcorr_self(m[, j])[480:735])
})

test_that(".vad_frame_features matches the original per-frame loop", {
  set.seed(21)
  fs <- 16000
  wave <- c(sin(2 * pi * 140 * (1:8000) / fs), rnorm(9000) * 0.1)
  win_l <- 480L; hop <- 160L
  win <- as.numeric(av::hanning(win_l))
  n_frames <- as.integer(floor((length(wave) - win_l - 1) / hop) + 1)
  got <- superassp:::.vad_frame_features(wave, win, hop, n_frames, 32, 256, fs, 12L, chunk = 37L)
  corr_win <- ref_xcorr_self(win)
  for (i in seq_len(n_frames)) {
    seg <- wave[(i - 1) * hop + 1:win_l] * win
    rxx <- (ref_xcorr_self(seg) / corr_win)[win_l:(2 * win_l - 1)]
    vec <- rxx[32:256]; mx <- max(vec)
    expect_equal(got$harm[i], log10(abs(mx / (rxx[1] - mx + .Machine$double.eps))))
    dxx <- 0.8 * sqrt(2 * abs(rxx[1] - rxx))[32:256]
    expect_equal(got$clarity[i], 1 - min(dxx) / (max(dxx) + .Machine$double.eps))
    expect_identical(got$hps[i], ref_hps(seg, fs))
  }
  lpc_err <- superassp:::.vad_yw_varpred(got$yw_acf, win_l)
  want <- vapply(seq_len(n_frames), function(i) {
    stats::ar.yw(wave[(i - 1) * hop + 1:win_l] * win, order.max = 12, aic = FALSE)$var.pred
  }, numeric(1))
  expect_equal(lpc_err, want)
})

test_that(".vad_yw_varpred reproduces stats::ar.yw var.pred", {
  set.seed(7)
  win <- as.numeric(av::hanning(480))
  segs <- cbind(
    rnorm(480) * win,
    sin(2 * pi * 150 * (0:479) / 16000) * win + rnorm(480, sd = 1e-3) * win,
    stats::filter(rnorm(480), 0.9, method = "recursive") * win,
    (rnorm(480) + 3) * win
  )
  acf_mat <- apply(segs, 2, function(s) {
    acf <- stats::acf(s, lag.max = 12, type = "covariance", plot = FALSE)$acf
    drop(acf)
  })
  got <- superassp:::.vad_yw_varpred(acf_mat, 480)
  want <- apply(segs, 2, function(s) stats::ar.yw(s, order.max = 12, aic = FALSE)$var.pred)
  expect_equal(got, want)
})

test_that("zero-variance frames still raise ar.yw's error", {
  wave <- c(rep(0, 2000), rnorm(4000))
  expect_error(superassp:::.vad_extract_features(wave, 16000), "zero-variance")
})

test_that(".vad_medfilt1 is identical to the per-element median filter", {
  set.seed(3)
  inputs <- list(rnorm(500), rnorm(11), rnorm(10), rnorm(3), 1, numeric(0),
                 rep(0.3, 30), round(rnorm(200)), c(rnorm(20), NA, rnorm(20)),
                 c(rep(0, 40), rnorm(40), rep(1, 40)))
  for (x in inputs) {
    expect_identical(superassp:::.vad_medfilt1(x, 11), ref_medfilt1(x, 11))
  }
  x <- rnorm(50)
  expect_identical(superassp:::.vad_medfilt1(x, 4), ref_medfilt1(x, 4))
})

test_that(".vad_apply_filtering_and_derivatives matches the scalar loops", {
  set.seed(11)
  for (dims in list(c(13, 400), c(4, 11), c(3, 12), c(3, 200))) {
    m <- matrix(rnorm(prod(dims)), dims[1])
    expect_equal(superassp:::.vad_apply_filtering_and_derivatives(m),
                 ref_filtering_and_derivatives(m))
  }
  m <- matrix(rnorm(3 * 8), 3)
  expect_error(ref_filtering_and_derivatives(m))
  expect_error(superassp:::.vad_apply_filtering_and_derivatives(m))
})

test_that(".vad_compute_mfcc matches the cbind loop", {
  set.seed(5)
  for (n in c(479L, 480L, 481L, 640L, 16000L)) {
    w <- rnorm(n)
    expect_identical(superassp:::.vad_compute_mfcc(w, 16000), ref_mfcc(w, 16000))
  }
})

test_that("HPS stays NA-propagating (known fft(seg)[1:1024] quirk preserved)", {
  set.seed(9)
  seg <- rnorm(480)
  expect_true(is.na(ref_hps(seg, 16000)))
  feat <- suppressWarnings(superassp:::.vad_extract_features(rnorm(16000), 16000))$mat_feat
  # hps row (row 17) is zeroed by the !is.finite() sweep
  expect_true(all(feat[17, ] == 0))
})

test_that(".polarity_frame_autocorr / .polarity_lpc match the scalar loop", {
  set.seed(13)
  order <- 18L
  fl <- 400L; sh <- 80L
  x <- rnorm(16000)
  w <- 0.5 * (1 - cos(2 * pi * (0:(fl - 1L)) / (fl - 1L)))
  n_ok <- floor((length(x) - fl) / sh) + 1L
  r_mat <- superassp:::.polarity_frame_autocorr(x, w, n_ok, fl, sh, order, chunk = 23L)
  for (i in seq_len(n_ok)) {
    fr <- x[(i - 1L) * sh + 1:fl] * w
    r_ref <- ref_polarity_lpc_autocorr(fr, order)
    expect_equal(r_mat[, i], r_ref, tolerance = 0)
    expect_identical(superassp:::.polarity_levinson(r_mat[, i], order),
                     superassp:::.polarity_lpc(fr, order))
  }
  expect_null(superassp:::.polarity_frame_autocorr(x, w[1:10], 3L, 10L, 5L, 18L))
})

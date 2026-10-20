# Pins the vectorised trk_peakslope() internals to the per-sample R loop and
# per-frame lm() fit they replaced (COVAREP peakslope port).

peakslope_filt_loop <- function(s, h_i) {
  s_padded <- c(s, rep(0, length(h_i) - 1L))
  y <- numeric(length(s_padded))
  for (i in seq_along(s_padded)) {
    start_idx <- max(1L, i - length(h_i) + 1L)
    h_idx_start <- length(h_i) - (i - start_idx)
    if (h_idx_start > 0) {
      y[i] <- sum(h_i[h_idx_start:length(h_i)] * s_padded[start_idx:i])
    }
  }
  half_len <- ceiling(length(h_i) / 2)
  y[half_len:(half_len + length(s) - 1L)]
}

peakslope_analysis_loop <- function(s, fs) {
  frame_len <- round(40 / 1000 * fs)
  frame_shift <- round(10 / 1000 * fs)
  times <- numeric()
  slopes <- numeric()
  i <- 0:6
  y <- matrix(0, nrow = length(i), ncol = length(s))
  for (n in seq_along(i)) {
    y[n, ] <- peakslope_filt_loop(s, superassp:::.peakslope_daless_mw(i, n, fs))
  }
  start <- 1L
  finish <- start + frame_len - 1L
  while (finish <= length(s)) {
    maxima <- apply(abs(y[, start:finish, drop = FALSE]), 1, max)
    maxima <- log10(pmax(rev(maxima), .Machine$double.eps))
    p <- stats::lm(maxima ~ seq_along(maxima))
    slopes <- c(slopes, unname(stats::coef(p)[2]))
    times <- c(times, (round((start + finish) / 2) - 1) / fs)
    start <- start + frame_shift
    finish <- start + frame_len - 1L
  }
  ps <- cbind(times, slopes)
  ps[is.nan(ps[, 2]), 2] <- 0
  unname(ps)
}

test_that("wavelet filter matches the per-sample loop at every scale", {
  set.seed(42)
  s <- rnorm(3000) * 1000
  for (n in 1:7) {
    h <- superassp:::.peakslope_daless_mw(0:6, n, 16000)
    expect_equal(superassp:::.peakslope_do_daless_filt(s, h),
                 peakslope_filt_loop(s, h), info = paste("scale", n))
  }
})

test_that("peakslope analysis matches the lm()-based frame loop", {
  set.seed(7)
  s <- c(rep(0, 400), rnorm(2600) * 500)   # includes digital silence
  expect_equal(superassp:::.peakslope_analysis(s, 16000),
               peakslope_analysis_loop(s, 16000))
  # Shorter than one frame: empty result, same shape
  expect_equal(dim(superassp:::.peakslope_analysis(rnorm(100), 16000)), c(0L, 2L))
})

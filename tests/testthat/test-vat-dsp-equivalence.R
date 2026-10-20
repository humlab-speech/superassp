# Pins the restructured VAT DSP kernels (perf rewrites) to small R
# re-implementations of the ORIGINAL algorithms they replaced.

ns <- asNamespace("superassp")

# Original O(P^2) findpeaks: strict local maxima >= height, sorted by height
# descending, greedy suppression within min distance, survivors ascending.
ref_findpeaks <- function(x, mph, mpd) {
  N <- length(x)
  raw <- which(x[2:(N - 1)] > x[1:(N - 2)] & x[2:(N - 1)] > x[3:N] &
                 x[2:(N - 1)] >= mph) + 1L
  raw <- raw[order(-x[raw])]  # stable; C++ std::sort only differs on exact ties
  kept <- rep(TRUE, length(raw))
  for (i in seq_along(raw)) {
    if (!kept[i]) next
    j <- seq_along(raw) > i & kept & abs(raw - raw[i]) < mpd
    kept[j] <- FALSE
  }
  sort(raw[kept])
}

test_that("vat findpeaks matches the original O(P^2) suppression", {
  set.seed(11)
  for (rep in 1:3) {
    x <- as.numeric(stats::filter(rnorm(3000), rep(1, 5) / 5, sides = 2))
    x[is.na(x)] <- 0
    for (d in c(1L, 2L, 7L, 40L, 300L)) {
      got <- ns$vat_findpeaks_cpp(x, -1e300, d)
      expect_identical(got$locs, ref_findpeaks(x, -Inf, d))
      expect_identical(as.numeric(got$pks), x[got$locs])
    }
    got <- ns$vat_findpeaks_cpp(x, 0.2, 9L)
    expect_identical(got$locs, ref_findpeaks(x, 0.2, 9L))
  }
})

# Original resample loop, including the k in (-p, 0) quirk where integer
# division truncates toward zero and adds h(k - lo) * x(0).
ref_resample <- function(x, p, q, beta = 5) {
  g <- function(a, b) if (b == 0) a else g(b, a %% b)
  d <- g(p, q); p <- p / d; q <- q / d
  L <- 2 * 10 * max(p, q) + 1
  h <- ns$vat_fir1_cpp(L - 1, 1 / max(p, q), "low") * p * ns$vat_kaiser_cpp(L, beta)
  half <- (L - 1) / 2
  Ly <- ceiling(length(x) * p / q)
  y <- numeric(Ly)
  for (n in 0:(Ly - 1)) {
    lo <- n * q - half
    acc <- 0
    for (k in lo:(n * q + half)) {
      if (k %% p != 0 && k >= 0) next
      xi <- trunc(k / p)
      if (xi >= 0 && xi < length(x)) acc <- acc + h[k - lo + 1] * x[xi + 1]
    }
    y[n + 1] <- acc
  }
  y
}

test_that("vat resample matches the original full-scan loop (incl. quirk)", {
  set.seed(12)
  x <- rnorm(120)
  for (pq in list(c(3, 2), c(2, 3), c(5, 7), c(4, 1), c(1, 3))) {
    expect_equal(ns$vat_resample_cpp(x, pq[1], pq[2]),
                 ref_resample(x, pq[1], pq[2]), tolerance = 1e-13)
  }
})

test_that("reused FFT plan is bit-identical to a fresh arma::fft", {
  set.seed(13)
  frames <- lapply(c(300, 17, 1200, 0, 5, 4410, 5000), rnorm)
  expect_true(ns$vat_fft_plan_identical_cpp(frames, 4410L))
  expect_true(ns$vat_fft_plan_identical_cpp(frames, 44100L))
  expect_true(ns$vat_fft_plan_identical_cpp(frames[1:3], 97L))
})

test_that("blocked FIR path of vat filter matches the direct sum", {
  set.seed(14)
  x <- rnorm(500)
  for (M in c(2, 9, 64, 705)) {
    b <- rnorm(M)
    ref <- vapply(seq_along(x), function(n) {
      k <- 0:min(n - 1, M - 1)
      sum(b[k + 1] * x[n - k])
    }, numeric(1))
    expect_equal(ns$vat_filter_cpp(b, 1, x), ref, tolerance = 1e-12)
  }
})

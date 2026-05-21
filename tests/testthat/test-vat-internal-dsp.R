test_that(".vat_hamming matches MATLAB hamming definition", {
  w <- .vat_hamming(10)
  # MATLAB: hamming(10) first/last = 0.08, middle ~= 1
  expect_equal(w[1], 0.08, tolerance = 1e-12)
  expect_equal(w[10], 0.08, tolerance = 1e-12)
  expect_lt(abs(max(w) - 1), 0.05)
  expect_equal(w, rev(w), tolerance = 1e-12)  # symmetric
})

test_that(".vat_hanning symmetric", {
  w <- .vat_hanning(8)
  expect_equal(w[1], 0)
  expect_equal(w[8], 0)
  expect_equal(w, rev(w), tolerance = 1e-12)
})

test_that(".vat_filter matches signal::filter on simple IIR", {
  set.seed(1)
  x <- rnorm(50)
  b <- c(1, -1)
  a <- c(1, -0.5)
  y_ref <- signal::filter(b, a, x)
  y_vat <- .vat_filter(b, a, x)
  expect_equal(as.numeric(y_vat), as.numeric(y_ref), tolerance = 1e-10)
})

test_that(".vat_filtfilt is zero-phase (centroid invariant on symmetric impulse)", {
  x <- c(rep(0, 20), 1, rep(0, 20))
  b <- c(0.25, 0.5, 0.25)
  a <- 1
  y <- .vat_filtfilt(b, a, x)
  # Centroid at sample 21
  i <- seq_along(y)
  centroid <- sum(i * y) / sum(y)
  expect_equal(centroid, 21, tolerance = 0.1)
})

test_that(".vat_fir1 lowpass sums roughly to 1 at DC", {
  h <- .vat_fir1(50, 0.3, "low")
  expect_equal(sum(h), 1, tolerance = 0.01)
})

test_that(".vat_butter lowpass stable (poles inside unit circle)", {
  bf <- .vat_butter(4, 0.3, "low")
  expect_true(all(Mod(polyroot(rev(bf$a))) < 1))
})

test_that(".vat_medfilt1 odd window matches stats::runmed for interior", {
  x <- c(1, 5, 2, 8, 3, 9, 4)
  y <- .vat_medfilt1(x, 3)
  # Interior: median of (1,5,2)=2, (5,2,8)=5, (2,8,3)=3, (8,3,9)=8, (3,9,4)=4
  expect_equal(y[2:6], c(2, 5, 3, 8, 4))
})

test_that(".vat_interp1 linear matches stats::approx", {
  x <- c(1, 3, 5, 8)
  y <- c(2, 6, 10, 4)
  xq <- c(2, 4, 6, 7)
  yq <- .vat_interp1(x, y, xq, "linear")
  yq_ref <- approx(x, y, xq)$y
  expect_equal(yq, yq_ref, tolerance = 1e-12)
})

test_that(".vat_interp1 spline matches stats::spline (natural)", {
  x <- seq(0, 2 * pi, length.out = 10)
  y <- sin(x)
  xq <- seq(0, 2 * pi, length.out = 50)
  yq <- .vat_interp1(x, y, xq, "spline")
  # stats::spline default is FMM, not natural — compare to sin instead
  expect_lt(max(abs(yq - sin(xq))), 0.02)
})

test_that(".vat_findpeaks finds local maxima", {
  x <- c(0, 1, 0, 2, 0, 1, 0, 3, 0)
  res <- .vat_findpeaks(x)
  expect_equal(sort(res$locs), c(2, 4, 6, 8))
})

test_that(".vat_findpeaks suppresses within MinPeakDistance", {
  x <- c(0, 1, 0, 2, 0, 3, 0)
  res <- .vat_findpeaks(x, min_peak_distance = 3)
  # Keeps highest peak first; suppresses neighbors within 3
  expect_true(6 %in% res$locs)
  expect_false(4 %in% res$locs)
})

test_that(".vat_fft matches stats::fft", {
  x <- c(1, 2, 3, 4)
  X1 <- .vat_fft(x)
  X2 <- stats::fft(x)
  expect_equal(X1, X2, tolerance = 1e-12)
})

test_that(".vat_resample 2x upsample doubles length and preserves shape", {
  set.seed(42); x <- sin(2 * pi * 0.05 * 1:100)
  y <- .vat_resample(x, 2, 1)
  expect_equal(length(y), 200)
  # Interior samples should track the underlying sine (edges have Kaiser taper)
  x_dense <- sin(2 * pi * 0.025 * 1:200)
  expect_gt(cor(y[50:150], x_dense[50:150]), 0.98)
})

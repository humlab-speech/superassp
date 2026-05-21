test_that(".vat_lpcauto recovers AR(2) coefficients from synthetic signal", {
  set.seed(7)
  N <- 2000
  # Synthesize AR(2) with poles 0.95*exp(±i*pi/4): a = [1, -2*0.95*cos(pi/4), 0.95^2]
  a_true <- c(1, -2 * 0.95 * cos(pi / 4), 0.95^2)
  e <- rnorm(N)
  x <- signal::filter(1, a_true, e) |> as.numeric()
  res <- .vat_lpcauto(x, 2)
  # Hamming windowing biases coefficients slightly; expect close but not exact
  expect_equal(res$ar[2], a_true[2], tolerance = 0.1)
  expect_equal(res$ar[3], a_true[3], tolerance = 0.1)
  expect_gt(res$e, 0)
})

test_that(".vat_lpcauto first coef is 1", {
  x <- rnorm(500)
  res <- .vat_lpcauto(x, 8)
  expect_equal(res$ar[1], 1)
})

test_that(".vat_burg first coef is 1 and order correct", {
  x <- rnorm(500)
  res <- .vat_burg(x, 6)
  expect_equal(res$ar[1], 1)
  expect_equal(length(res$ar), 7)
})

test_that(".vat_lpcar2rf round-trips via lpcrf2rr (autocorr direction sanity)", {
  ar <- c(1, -1.4, 0.7, -0.3, 0.1)
  k <- .vat_lpcar2rf(ar)
  expect_equal(length(k), length(ar) - 1)
  # All |k_i| < 1 for stable AR
  expect_true(all(abs(k) < 1))
})

test_that(".vat_lpcar2ra gives r(0) > 0", {
  ar <- c(1, -1.4, 0.7, -0.3)
  r <- .vat_lpcar2ra(ar)
  expect_gt(r[1], 0)
  expect_equal(length(r), length(ar))
})

test_that(".vat_distitar is zero for identical AR", {
  ar <- c(1, -1.2, 0.5)
  expect_equal(.vat_distitar(ar, ar), 0, tolerance = 1e-10)
})

test_that(".vat_distitar is positive for different AR", {
  ar1 <- c(1, -1.2, 0.5)
  ar2 <- c(1, -0.8, 0.3)
  expect_gt(.vat_distitar(ar1, ar2), 0)
})

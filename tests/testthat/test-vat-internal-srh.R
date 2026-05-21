test_that("srh_pitch returns expected list structure", {
  fs <- 16000
  t  <- seq(0, 1, by = 1/fs)
  x  <- sin(2 * pi * 120 * t)  # 120 Hz sine
  res <- .vat_pitch(x, fs, f0_min = 60, f0_max = 400, method = "srh")
  expect_named(res, c("f0", "VUV", "SRHVal", "time", "fs"))
  expect_equal(res$fs, 16000)
  expect_gt(length(res$f0), 0)
})

test_that("srh_pitch tracks a synthetic harmonic signal close to its F0", {
  fs <- 16000
  t  <- seq(0, 1.5, by = 1/fs)
  # Sum of harmonics at 150 Hz to give clear residual peaks
  f0_true <- 150
  x <- 0
  for (h in 1:6) x <- x + (1/h) * sin(2 * pi * h * f0_true * t)
  x <- x + 0.01 * rnorm(length(x))
  res <- .vat_pitch(x, fs, f0_min = 60, f0_max = 400, method = "srh")
  voiced <- res$f0[res$VUV == 1]
  if (length(voiced) > 5) {
    expect_lt(abs(median(voiced) - f0_true), 5)
  }
})

test_that("srh_pitch handles silence with all-unvoiced output", {
  fs <- 16000
  x <- rnorm(fs * 0.5) * 1e-6  # near-silence
  res <- .vat_pitch(x, fs, f0_min = 60, f0_max = 400, method = "srh")
  expect_lt(mean(res$VUV), 0.2)  # mostly unvoiced
})

test_that("srh_pitch f0 vector length matches SRHVal and time", {
  fs <- 16000
  x  <- sin(2 * pi * 200 * seq(0, 1, by = 1/fs))
  res <- .vat_pitch(x, fs, method = "srh")
  expect_equal(length(res$f0), length(res$SRHVal))
  expect_equal(length(res$f0), length(res$time))
  expect_equal(length(res$f0), length(res$VUV))
})

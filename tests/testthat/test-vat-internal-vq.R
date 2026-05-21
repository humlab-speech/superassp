test_that(".vat_voice_quality returns expected output on synthetic glottal flow", {
  fs <- 16000
  f0_true <- 120
  T0 <- round(fs / f0_true)
  n_cycles <- 50
  N <- T0 * n_cycles
  # Synthesize a triangular glottal flow per cycle
  glot <- numeric(N)
  for (k in 0:(n_cycles - 1)) {
    a <- k * T0 + 1
    b <- a + T0 - 1
    open_phase <- round(0.6 * T0)
    pulse <- c(seq(0, 1, length.out = open_phase),
               seq(1, 0, length.out = T0 - open_phase))
    glot[a:b] <- pulse[1:T0]
  }
  # Derivative
  dg <- c(0, diff(glot)) * fs
  GCI <- seq(T0, n_cycles * T0, by = T0)
  res <- .vat_voice_quality(dg, fs, GCI, backend = "cpp")
  expect_named(res, c("NAQ", "QOQ", "H1H2", "HRF"))
  expect_equal(length(res$NAQ), length(GCI))
  # NAQ should be in plausible range (0.05-0.5)
  naq_valid <- res$NAQ[res$NAQ != 0]
  if (length(naq_valid) > 5) {
    expect_gt(median(naq_valid), 0)
    expect_lt(median(naq_valid), 1)
  }
})

test_that(".vat_peak_slope returns a numeric vector of expected length", {
  fs <- 16000
  x <- sin(2 * pi * 120 * seq(0, 1, by = 1/fs)) +
       0.3 * sin(2 * pi * 800 * seq(0, 1, by = 1/fs))
  ps <- .vat_peak_slope(x, fs, backend = "cpp")
  expect_true(is.numeric(ps))
  expect_gt(length(ps), 50)
  expect_true(all(is.finite(ps)))
})

test_that(".vat_mdq returns one value per GCI", {
  fs <- 16000
  res <- rnorm(fs)
  GCI <- seq(100, length(res) - 100, by = 130)
  mdq <- .vat_mdq(res, fs, GCI)
  expect_equal(length(mdq), length(GCI))
  expect_true(all(is.finite(mdq)))
  # MDQ values should typically be in [0, 1]
  expect_true(all(mdq >= 0 & mdq <= 2))
})

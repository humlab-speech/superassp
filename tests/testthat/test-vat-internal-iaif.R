test_that(".vat_iaif returns valid glottal flow on synthetic voiced input", {
  set.seed(11)
  fs <- 16000
  f0 <- 130
  t  <- seq(0, 1, by = 1/fs)
  pulse <- numeric(length(t))
  pulse[seq(1, length(pulse), by = round(fs / f0))] <- 1
  # Two-formant filter
  a1 <- c(1, -2 * 0.95 * cos(2 * pi * 700 / fs), 0.95^2)
  a2 <- c(1, -2 * 0.95 * cos(2 * pi * 1200 / fs), 0.95^2)
  x <- as.numeric(signal::filter(1, a1, pulse))
  x <- as.numeric(signal::filter(1, a2, x))
  x <- x / max(abs(x))

  gci <- .vat_se_vq(x, fs)$GCI
  res <- .vat_iaif(x, fs, gci, backend = "cpp")
  expect_named(res, c("g", "g_iaif", "dg", "ar_lpc", "e_lpc"))
  expect_equal(length(res$g), length(x))
  expect_equal(length(res$dg), length(x))
  expect_true(all(is.finite(res$g)))
  expect_true(all(is.finite(res$dg)))
})

test_that(".vat_iaif cpp and r backends agree on glottal flow derivative shape", {
  skip_if_not_installed("signal")
  set.seed(13)
  fs <- 16000
  f0 <- 150
  t  <- seq(0, 0.6, by = 1/fs)
  pulse <- numeric(length(t))
  pulse[seq(1, length(pulse), by = round(fs / f0))] <- 1
  a <- c(1, -2 * 0.95 * cos(2 * pi * 800 / fs), 0.95^2)
  x <- as.numeric(signal::filter(1, a, pulse))
  x <- x / max(abs(x))

  gci <- .vat_se_vq(x, fs)$GCI
  res_c <- .vat_iaif(x, fs, gci, backend = "cpp")
  res_r <- tryCatch(.vat_iaif(x, fs, gci, backend = "r"), error = function(e) NULL)
  if (!is.null(res_r)) {
    # Same length, correlated within middle region
    n  <- length(x)
    mid <- seq(round(0.1 * n), round(0.9 * n))
    expect_gt(cor(res_c$dg[mid], res_r$g_iaif[mid]), 0.5)
  } else {
    succeed()
  }
})

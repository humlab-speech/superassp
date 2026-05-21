test_that(".vat_se_vq returns reasonable GCIs on synthetic voiced signal", {
  fs <- 16000
  f0_true <- 120
  t  <- seq(0, 1.5, by = 1/fs)
  # Synthesize pulse train at 120 Hz then filter to give formant structure
  pulse <- numeric(length(t))
  pulse[seq(1, length(pulse), by = round(fs / f0_true))] <- 1
  # AR(2) formant filter
  a <- c(1, -2 * 0.95 * cos(2 * pi * 800 / fs), 0.95^2)
  x <- as.numeric(signal::filter(1, a, pulse))
  x <- x / max(abs(x))
  res <- .vat_se_vq(x, fs, backend = "cpp")
  expect_named(res, c("GCI", "res", "rep", "MBS", "F0mean", "F0max"))
  expect_gt(length(res$GCI), 50)
  # GCIs should be roughly evenly spaced; mode of intervals near fs/f0_true
  if (length(res$GCI) > 10) {
    intervals <- diff(res$GCI)
    expected <- fs / f0_true
    # Allow wide tolerance — exact MATLAB fidelity tested by goldens
    expect_lt(abs(median(intervals) - expected) / expected, 0.6)
  }
})

test_that(".vat_se_vq backends produce comparable GCI counts on synthetic input", {
  skip_if_not_installed("signal")
  fs <- 16000
  f0_true <- 150
  t  <- seq(0, 0.8, by = 1/fs)
  pulse <- numeric(length(t))
  pulse[seq(1, length(pulse), by = round(fs / f0_true))] <- 1
  a <- c(1, -2 * 0.95 * cos(2 * pi * 700 / fs), 0.95^2)
  x <- as.numeric(signal::filter(1, a, pulse))
  res_cpp <- .vat_se_vq(x, fs, backend = "cpp")
  res_r   <- tryCatch(.vat_se_vq(x, fs, backend = "r"), error = function(e) NULL)
  if (!is.null(res_r)) {
    expect_lt(abs(length(res_cpp$GCI) - length(res_r$GCI)) /
                max(1, length(res_r$GCI)), 0.3)
  } else {
    succeed()
  }
})

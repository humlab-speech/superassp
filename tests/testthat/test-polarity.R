test_that("lst_polarity handles single file", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_polarity(test_wav, verbose = FALSE)

  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 2)
  expect_true("file" %in% names(result))
  expect_true("polarity" %in% names(result))
  expect_equal(nrow(result), 1)
  expect_true(result$polarity[1] %in% c(-1L, 1L))
})

test_that("lst_polarity handles multiple files", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  files <- c(test_wav, test_wav)
  result <- lst_polarity(files, verbose = FALSE)

  expect_equal(nrow(result), 2)
  expect_true(all(result$polarity %in% c(-1L, 1L)))
})

test_that("lst_polarity returns consistent polarity", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Run twice on same file, should get same result (deterministic)
  result1 <- lst_polarity(test_wav, verbose = FALSE)
  result2 <- lst_polarity(test_wav, verbose = FALSE)

  expect_equal(result1$polarity[1], result2$polarity[1])
})

test_that("lst_polarity respects time windowing", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Full file
  result_full <- lst_polarity(test_wav, verbose = FALSE)

  # Partial file
  result_partial <- lst_polarity(test_wav, beginTime = 0, endTime = 0.5, verbose = FALSE)

  # Both should have valid polarity (even if different due to different data)
  expect_true(result_full$polarity[1] %in% c(-1L, 1L))
  expect_true(result_partial$polarity[1] %in% c(-1L, 1L))
})

test_that("lst_polarity rejects toFile parameter", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  expect_error(lst_polarity(test_wav, toFile = TRUE, verbose = FALSE))
})

test_that(".polarity_lpc_residual_two_signals matches direct-append reference", {
  set.seed(42)
  filter_signal <- rnorm(2000)
  analysis_signal <- rnorm(1900)

  # Reference: the pre-fix O(n^2) direct-append implementation.
  reference_impl <- function(filter_signal, analysis_signal, frame_length, frame_shift, order) {
    n_frames <- floor((length(filter_signal) - frame_length) / frame_shift) + 1L
    residuals <- numeric()
    for (i in seq_len(n_frames)) {
      start_idx <- (i - 1L) * frame_shift + 1L
      end_idx <- start_idx + frame_length - 1L
      if (end_idx > length(filter_signal) || end_idx > length(analysis_signal)) break
      frame_filt <- filter_signal[start_idx:end_idx]
      w <- 0.5 * (1 - cos(2 * pi * (0:(frame_length - 1L)) / (frame_length - 1L)))
      frame_filt_windowed <- frame_filt * w
      a <- superassp:::.polarity_lpc(frame_filt_windowed, order)
      frame_ana <- analysis_signal[start_idx:end_idx]
      filter_b <- if (length(a) > 1) c(1, -a[-1]) else 1
      res_frame <- stats::filter(frame_ana, filter_b, method = "convolution", sides = 1)
      residuals <- c(residuals, res_frame[!is.na(res_frame)])
    }
    residuals
  }

  expected <- reference_impl(filter_signal, analysis_signal, 400L, 100L, 12L)
  actual <- superassp:::.polarity_lpc_residual_two_signals(
    filter_signal, analysis_signal, 400L, 100L, 12L
  )
  expect_equal(actual, expected)
})

test_that(".polarity_lpc_residual_two_signals applies filter_b as a causal FIR filter to analysis_signal", {
  set.seed(7)
  filter_signal <- rnorm(600)
  analysis_signal <- rnorm(550)

  # Independent reference: manual causal FIR convolution, no stats::filter call at
  # all, so it cannot share the production code's argument-order/method bug.
  reference_convolution <- function(filter_signal, analysis_signal, frame_length, frame_shift, order) {
    n_frames <- floor((length(filter_signal) - frame_length) / frame_shift) + 1L
    residuals <- numeric()
    for (i in seq_len(n_frames)) {
      start_idx <- (i - 1L) * frame_shift + 1L
      end_idx <- start_idx + frame_length - 1L
      if (end_idx > length(filter_signal) || end_idx > length(analysis_signal)) break
      frame_filt <- filter_signal[start_idx:end_idx]
      w <- 0.5 * (1 - cos(2 * pi * (0:(frame_length - 1L)) / (frame_length - 1L)))
      frame_filt_windowed <- frame_filt * w
      a <- superassp:::.polarity_lpc(frame_filt_windowed, order)
      frame_ana <- analysis_signal[start_idx:end_idx]
      filter_b <- if (length(a) > 1) c(1, -a[-1]) else 1
      p <- length(filter_b)
      n <- length(frame_ana)
      y <- rep(NA_real_, n)
      for (t in seq_len(n)) {
        if (t >= p) y[t] <- sum(filter_b * rev(frame_ana[(t - p + 1):t]))
      }
      residuals <- c(residuals, y[!is.na(y)])
    }
    residuals
  }

  expected <- reference_convolution(filter_signal, analysis_signal, 200L, 50L, 10L)
  actual <- superassp:::.polarity_lpc_residual_two_signals(
    filter_signal, analysis_signal, 200L, 50L, 10L
  )
  expect_equal(actual, expected, tolerance = 1e-10)
})

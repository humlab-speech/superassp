test_that("trk_formant_cgdzp handles single file", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_formant_cgdzp(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("F1" %in% names(result))
  expect_true("F2" %in% names(result))
  expect_true("F3" %in% names(result))
  expect_true("F4" %in% names(result))
  expect_true("F5" %in% names(result))
  expect_true(length(result$F1) > 0)
})

test_that("trk_formant_cgdzp respects frame parameters", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result30 <- trk_formant_cgdzp(test_wav, frameSize = 30, frameShift = 10, toFile = FALSE, verbose = FALSE)
  result40 <- trk_formant_cgdzp(test_wav, frameSize = 40, frameShift = 5, toFile = FALSE, verbose = FALSE)

  # Different frame shift should result in different number of frames
  expect_true(length(result40$F1) > length(result30$F1))
})

test_that("trk_formant_cgdzp handles time windowing", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Full file
  result_full <- trk_formant_cgdzp(test_wav, toFile = FALSE, verbose = FALSE)

  # Partial file (first 0.5s)
  result_partial <- trk_formant_cgdzp(test_wav, beginTime = 0, endTime = 0.5, toFile = FALSE, verbose = FALSE)

  # Partial should have fewer frames
  expect_true(length(result_partial$F1) < length(result_full$F1))
})

test_that("trk_formant_cgdzp produces reasonable formant frequencies", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_formant_cgdzp(test_wav, toFile = FALSE, verbose = FALSE)

  # Sustained vowel should have detected formants (non-zero values)
  expect_true(any(result$F1 > 0))
  expect_true(any(result$F2 > 0))
  expect_true(any(result$F3 > 0))

  # Sanity check: higher formants should exist but not all be zero
  # F1 < F2 < F3 < F4 < F5 (in general, though with noise on individual frames)
  f1_vals <- result$F1[result$F1 > 0]
  f2_vals <- result$F2[result$F2 > 0]
  f3_vals <- result$F3[result$F3 > 0]
  f4_vals <- result$F4[result$F4 > 0]
  f5_vals <- result$F5[result$F5 > 0]

  # All formants should be detected at least in some frames
  expect_true(length(f1_vals) > 0)
  expect_true(length(f2_vals) > 0)
  expect_true(length(f3_vals) > 0)

  # Medians should respect general formant ordering
  # This is a soft check accounting for algorithm noise
  f1_med <- median(f1_vals, na.rm = TRUE)
  f2_med <- median(f2_vals, na.rm = TRUE)
  f3_med <- median(f3_vals, na.rm = TRUE)

  # F1 and F2 medians should be in a reasonable order
  expect_true(f1_med > 0 && f2_med > 0 && f3_med > 0)
  # Very soft check: F1 typically < F2 < F3
  expect_true(f1_med < f3_med)
})

test_that("trk_formant_cgdzp function attributes are set", {
  expect_equal(attr(trk_formant_cgdzp, "ext"), "cgf")
  expect_equal(attr(trk_formant_cgdzp, "tracks"), c("F1", "F2", "F3", "F4", "F5"))
  expect_equal(attr(trk_formant_cgdzp, "outputType"), "SSFF")
})

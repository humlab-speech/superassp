test_that("trk_peakslope handles single file", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_peakslope(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("peakslope" %in% names(result))
  expect_true(length(result$peakslope) > 0)
})

test_that("trk_peakslope produces numeric slopes", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_peakslope(test_wav, toFile = FALSE, verbose = FALSE)

  # Should have valid numeric values
  expect_true(is.numeric(result$peakslope))
  expect_true(all(is.finite(result$peakslope)))
})

test_that("trk_peakslope handles multiple files", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  files <- c(test_wav, test_wav)
  result <- trk_peakslope(files, toFile = FALSE, verbose = FALSE)

  expect_is(result, "list")
  expect_equal(length(result), 2)
  expect_s3_class(result[[1]], "AsspDataObj")
  expect_s3_class(result[[2]], "AsspDataObj")
})

test_that("trk_peakslope respects time windowing", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Full file
  result_full <- trk_peakslope(test_wav, toFile = FALSE, verbose = FALSE)

  # Partial file (first 0.5s)
  result_partial <- trk_peakslope(test_wav, beginTime = 0, endTime = 0.5, toFile = FALSE, verbose = FALSE)

  # Partial should have fewer frames
  expect_true(length(result_partial$peakslope) < length(result_full$peakslope))
})

test_that("trk_peakslope produces reasonable slope values", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_peakslope(test_wav, toFile = FALSE, verbose = FALSE)

  # Spectral slopes should be finite and in reasonable range
  # Typically between -1 and 1 for natural speech
  slopes <- result$peakslope[is.finite(result$peakslope)]
  expect_true(length(slopes) > 0)

  # Most slopes should be in reasonable range (allow some outliers)
  # Negative slope = bright (high frequency energy)
  # Positive slope = dull (low frequency energy)
  expect_true(mean(slopes) > -2 && mean(slopes) < 2)
})

test_that("trk_peakslope function attributes are set", {
  expect_equal(attr(trk_peakslope, "ext"), "psl")
  expect_equal(attr(trk_peakslope, "tracks"), c("peakslope"))
  expect_equal(attr(trk_peakslope, "outputType"), "SSFF")
})

test_that("trk_peakslope returns consistent results", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Run twice on same file
  result1 <- trk_peakslope(test_wav, toFile = FALSE, verbose = FALSE)
  result2 <- trk_peakslope(test_wav, toFile = FALSE, verbose = FALSE)

  # Results should be identical (deterministic algorithm)
  expect_equal(result1$peakslope, result2$peakslope)
})

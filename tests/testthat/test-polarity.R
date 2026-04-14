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

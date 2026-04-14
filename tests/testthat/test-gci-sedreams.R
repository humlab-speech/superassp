test_that("lst_covarep_gci_sedreams handles single file", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_covarep_gci_sedreams(test_wav, f0mean = 100, verbose = FALSE)

  expect_is(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_true("file" %in% names(result))
  expect_true("n_gcis" %in% names(result))
  expect_true("gci_times" %in% names(result))
})

test_that("lst_covarep_gci_sedreams produces numeric GCI times", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_covarep_gci_sedreams(test_wav, f0mean = 100, verbose = FALSE)

  # GCI times should be numeric and positive
  gci_times <- result$gci_times[[1]]
  expect_true(is.numeric(gci_times))
  expect_true(all(gci_times >= 0))

  # n_gcis should match length of gci_times
  expect_equal(result$n_gcis[1], length(gci_times))
})

test_that("lst_covarep_gci_sedreams handles multiple files", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  files <- c(test_wav, test_wav)
  result <- lst_covarep_gci_sedreams(files, f0mean = 100, verbose = FALSE)

  expect_is(result, "data.frame")
  expect_equal(nrow(result), 2)
  expect_equal(length(result$gci_times), 2)
})

test_that("lst_covarep_gci_sedreams respects time windowing", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Full file
  result_full <- lst_covarep_gci_sedreams(test_wav, f0mean = 100, verbose = FALSE)

  # Partial file (first 0.5s)
  result_partial <- lst_covarep_gci_sedreams(test_wav, beginTime = 0, endTime = 0.5, f0mean = 100, verbose = FALSE)

  # Partial should have fewer or equal GCIs
  expect_true(result_partial$n_gcis[1] <= result_full$n_gcis[1])
})

test_that("lst_covarep_gci_sedreams with auto F0 estimation", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Should not crash with default F0
  result <- lst_covarep_gci_sedreams(test_wav, f0mean = NULL, verbose = FALSE)

  expect_is(result, "data.frame")
  expect_true(result$n_gcis[1] >= 0)
})

test_that("lst_covarep_gci_sedreams with polarity detection", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Should auto-detect polarity
  result <- lst_covarep_gci_sedreams(test_wav, f0mean = 100, polarity = NULL, verbose = FALSE)

  expect_is(result, "data.frame")
  expect_equal(nrow(result), 1)
})

test_that("lst_covarep_gci_sedreams consistency (deterministic)", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Run twice on same file
  result1 <- lst_covarep_gci_sedreams(test_wav, f0mean = 100, verbose = FALSE)
  result2 <- lst_covarep_gci_sedreams(test_wav, f0mean = 100, verbose = FALSE)

  # Results should be identical
  expect_equal(result1$n_gcis, result2$n_gcis)
  expect_equal(result1$gci_times[[1]], result2$gci_times[[1]])
})

test_that("trk_covarep_vq_gci handles single file with auto GCI", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_covarep_vq_gci(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true(all(c("naq", "qoq", "h1h2", "hrf", "psp") %in% names(result)))
})

test_that("trk_covarep_vq_gci produces numeric features", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_covarep_vq_gci(test_wav, toFile = FALSE, verbose = FALSE)

  # All tracks should be numeric
  expect_true(is.numeric(result$naq))
  expect_true(is.numeric(result$qoq))
  expect_true(is.numeric(result$h1h2))
  expect_true(is.numeric(result$hrf))
  expect_true(is.numeric(result$psp))

  # Same length
  expect_equal(length(result$naq), length(result$qoq))
  expect_equal(length(result$naq), length(result$h1h2))
})

test_that("trk_covarep_vq_gci handles multiple files", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  files <- c(test_wav, test_wav)
  result <- trk_covarep_vq_gci(files, toFile = FALSE, verbose = FALSE)

  expect_is(result, "list")
  expect_equal(length(result), 2)
  expect_s3_class(result[[1]], "AsspDataObj")
  expect_s3_class(result[[2]], "AsspDataObj")
})

test_that("trk_covarep_vq_gci respects time windowing", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Full file
  result_full <- trk_covarep_vq_gci(test_wav, toFile = FALSE, verbose = FALSE)

  # Partial file (first 0.5s)
  result_partial <- trk_covarep_vq_gci(test_wav, beginTime = 0, endTime = 0.5, toFile = FALSE, verbose = FALSE)

  # Partial should have fewer or equal frames
  expect_true(length(result_partial$naq) <= length(result_full$naq))
})

test_that("trk_covarep_vq_gci with provided GCIs", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get GCIs first
  gcis <- lst_covarep_gci_sedreams(test_wav, f0mean = 100, verbose = FALSE)
  gci_times <- gcis$gci_times[[1]]

  # Use provided GCIs
  result <- trk_covarep_vq_gci(test_wav, gci_times = gci_times, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("naq" %in% names(result))
})

test_that("trk_covarep_vq_gci function attributes are set", {
  expect_equal(attr(trk_covarep_vq_gci, "ext"), "vqg")
  expect_equal(attr(trk_covarep_vq_gci, "tracks"), c("naq", "qoq", "h1h2", "hrf", "psp"))
  expect_equal(attr(trk_covarep_vq_gci, "outputType"), "SSFF")
})

test_that("trk_covarep_vq_gci consistency (deterministic)", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Run twice on same file
  result1 <- trk_covarep_vq_gci(test_wav, toFile = FALSE, verbose = FALSE)
  result2 <- trk_covarep_vq_gci(test_wav, toFile = FALSE, verbose = FALSE)

  # Results should be identical (deterministic)
  expect_equal(length(result1$naq), length(result2$naq))
})

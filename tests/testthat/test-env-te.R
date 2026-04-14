test_that("trk_covarep_env_te handles single file", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_covarep_env_te(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true(all(c("env_te", "env_cc") %in% names(result)))
})

test_that("trk_covarep_env_te produces numeric spectra", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_covarep_env_te(test_wav, toFile = FALSE, verbose = FALSE)

  expect_true(is.numeric(result$env_te))
  expect_true(is.numeric(result$env_cc))
  expect_true(nrow(result$env_te) > 0)
})

test_that("trk_covarep_env_te handles multiple files", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  files <- c(test_wav, test_wav)
  result <- trk_covarep_env_te(files, toFile = FALSE, verbose = FALSE)

  expect_is(result, "list")
  expect_equal(length(result), 2)
  expect_s3_class(result[[1]], "AsspDataObj")
})

test_that("trk_covarep_env_te respects frame parameters", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result1 <- trk_covarep_env_te(test_wav, frameShift = 5, toFile = FALSE, verbose = FALSE)
  result2 <- trk_covarep_env_te(test_wav, frameShift = 10, toFile = FALSE, verbose = FALSE)

  # Smaller frame shift should give more frames
  expect_true(nrow(result1$env_te) >= nrow(result2$env_te))
})

test_that("trk_covarep_env_te function attributes are set", {
  expect_equal(attr(trk_covarep_env_te, "ext"), "ete")
  expect_equal(attr(trk_covarep_env_te, "tracks"), c("env_te", "env_cc"))
  expect_equal(attr(trk_covarep_env_te, "outputType"), "SSFF")
})

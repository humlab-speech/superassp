test_that("trk_hmpd handles single file", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_hmpd(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true(all(c("ae", "pdm", "pdd") %in% names(result)))
  expect_true(nrow(result$ae) > 0)
})

test_that("trk_hmpd produces numeric features", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_hmpd(test_wav, toFile = FALSE, verbose = FALSE)

  # All should be numeric and finite
  expect_true(is.numeric(result$ae))
  expect_true(is.numeric(result$pdm))
  expect_true(is.numeric(result$pdd))

  # Check dimensions: compressed representation
  expect_equal(ncol(result$ae), 25)  # 1 + 24
  expect_equal(ncol(result$pdm), 25) # 1 + 24
  expect_equal(ncol(result$pdd), 13) # 1 + 12

  # Mostly finite values
  expect_true(sum(is.finite(result$ae)) / length(result$ae) > 0.8)
})

test_that("trk_hmpd handles multiple files", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  files <- c(test_wav, test_wav)
  result <- trk_hmpd(files, toFile = FALSE, verbose = FALSE)

  expect_is(result, "list")
  expect_equal(length(result), 2)
  expect_s3_class(result[[1]], "AsspDataObj")
  expect_s3_class(result[[2]], "AsspDataObj")
})

test_that("trk_hmpd respects time windowing", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Full file
  result_full <- trk_hmpd(test_wav, toFile = FALSE, verbose = FALSE)

  # Partial file (first 0.5s)
  result_partial <- trk_hmpd(test_wav, beginTime = 0, endTime = 0.5, toFile = FALSE, verbose = FALSE)

  # Partial should have fewer frames
  expect_true(nrow(result_partial$ae) < nrow(result_full$ae))
})

test_that("trk_hmpd function attributes are set", {
  expect_equal(attr(trk_hmpd, "ext"), "hpd")
  expect_equal(attr(trk_hmpd, "tracks"), c("ae", "pdm", "pdd"))
  expect_equal(attr(trk_hmpd, "outputType"), "SSFF")
})

test_that("trk_hmpd with external F0s", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Compute F0 externally
  f0_res <- trk_srh_variant(test_wav, toFile = FALSE, verbose = FALSE)
  f0s <- cbind(
    seq(0, (nrow(f0_res$f0) - 1) * 0.01, by = 0.01),
    f0_res$f0[, 1]
  )

  # Use external F0s
  result <- trk_hmpd(test_wav, f0s = f0s, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true(all(c("ae", "pdm", "pdd") %in% names(result)))
})

test_that("trk_hmpd consistency (deterministic)", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Run twice on same file
  result1 <- trk_hmpd(test_wav, toFile = FALSE, verbose = FALSE)
  result2 <- trk_hmpd(test_wav, toFile = FALSE, verbose = FALSE)

  # Results should be identical (except for possible numeric rounding)
  expect_equal(nrow(result1$ae), nrow(result2$ae))
  expect_equal(ncol(result1$ae), ncol(result2$ae))
})

test_that("trk_hmpd produces non-zero features", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_hmpd(test_wav, toFile = FALSE, verbose = FALSE)

  # Amplitude envelope should have non-zero values
  ae_nonzero <- sum(result$ae[is.finite(result$ae)] != 0)
  expect_true(ae_nonzero > nrow(result$ae) * 0.5)  # At least 50% non-zero

  # PDM should have some non-zero values
  pdm_nonzero <- sum(result$pdm[is.finite(result$pdm)] != 0)
  expect_true(pdm_nonzero > 0)
})

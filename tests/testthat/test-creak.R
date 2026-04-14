test_that("trk_covarep_creak handles single file", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_covarep_creak(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true(all(c("creak_pp", "creak_bin") %in% names(result)))
})

test_that("trk_covarep_creak produces numeric creakiness values", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_covarep_creak(test_wav, toFile = FALSE, verbose = FALSE)

  expect_true(is.numeric(result$creak_pp))
  expect_true(is.numeric(result$creak_bin))
  expect_true(all(result$creak_pp >= 0 & result$creak_pp <= 1))
  expect_true(all(result$creak_bin %in% c(0, 1)))
})

test_that("trk_covarep_creak handles multiple files", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  files <- c(test_wav, test_wav)
  result <- trk_covarep_creak(files, toFile = FALSE, verbose = FALSE)

  expect_is(result, "list")
  expect_equal(length(result), 2)
  expect_s3_class(result[[1]], "AsspDataObj")
})

test_that("trk_covarep_creak function attributes are set", {
  expect_equal(attr(trk_covarep_creak, "ext"), "crk")
  expect_equal(attr(trk_covarep_creak, "tracks"), c("creak_pp", "creak_bin"))
  expect_equal(attr(trk_covarep_creak, "outputType"), "SSFF")
})

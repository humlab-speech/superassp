test_that("trk_covarep_vad_drugman handles single file", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_covarep_vad_drugman(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true(all(c("vad_final", "vad_mfcc", "vad_sadjadi", "vad_new") %in% names(result)))
})

test_that("trk_covarep_vad_drugman produces numeric posteriors", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_covarep_vad_drugman(test_wav, toFile = FALSE, verbose = FALSE)

  expect_true(is.numeric(result$vad_final))
  expect_true(is.numeric(result$vad_mfcc))
  expect_true(all(result$vad_final >= 0 & result$vad_final <= 1))
  expect_true(all(result$vad_mfcc >= 0 & result$vad_mfcc <= 1))
})

test_that("trk_covarep_vad_drugman handles multiple files", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  files <- c(test_wav, test_wav)
  result <- trk_covarep_vad_drugman(files, toFile = FALSE, verbose = FALSE)

  expect_is(result, "list")
  expect_equal(length(result), 2)
  expect_s3_class(result[[1]], "AsspDataObj")
})

test_that("trk_covarep_vad_drugman function attributes are set", {
  expect_equal(attr(trk_covarep_vad_drugman, "ext"), "cvd")
  expect_equal(attr(trk_covarep_vad_drugman, "tracks"),
               c("vad_final", "vad_mfcc", "vad_sadjadi", "vad_new"))
  expect_equal(attr(trk_covarep_vad_drugman, "outputType"), "SSFF")
})

test_that("trk_covarep_vad_drugman respects time windowing", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result_full <- trk_covarep_vad_drugman(test_wav, toFile = FALSE, verbose = FALSE)
  result_partial <- trk_covarep_vad_drugman(test_wav, beginTime = 0, endTime = 0.5,
                                            toFile = FALSE, verbose = FALSE)

  expect_true(length(result_partial$vad_final) < length(result_full$vad_final))
})

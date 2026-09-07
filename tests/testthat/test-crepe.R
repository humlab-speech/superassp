test_that("trk_pitch_crepe returns valid AsspDataObj for sustained /a/", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  skip_if(!superassp:::ort_available_cpp(), "ONNX Runtime not available")
  skip_on_cran()
  skip_if_not_installed("huggingfaceR")
  skip_if_offline()

  result <- trk_pitch_crepe(test_wav, model = "tiny", toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true(all(c("f0", "periodicity") %in% names(result)))
  expect_true(nrow(result$f0) > 0L)

  voiced <- result$f0[result$f0 > 0]
  expect_true(length(voiced) > 0)
  expect_gt(median(voiced), 50)
  expect_lt(median(voiced), 550)

  expect_true(all(result$periodicity >= 0, na.rm = TRUE))
  expect_true(all(result$periodicity <= 1, na.rm = TRUE))

  # Frame rate 100 Hz (default 10 ms windowShift)
  expect_equal(attr(result, "sampleRate"), 100.0)
})

test_that("trk_pitch_crepe writes SSFF file with .crp extension", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  skip_if(!superassp:::ort_available_cpp(), "ONNX Runtime not available")
  skip_on_cran()
  skip_if_not_installed("huggingfaceR")
  skip_if_offline()

  out_dir <- tempdir()
  n <- trk_pitch_crepe(test_wav, model = "tiny", toFile = TRUE,
                        outputDirectory = out_dir, verbose = FALSE)
  expect_equal(n, 1L)
  out_file <- file.path(out_dir, "a1.crp")
  expect_true(file.exists(out_file))

  reread <- read_ssff(out_file)
  expect_true(all(c("f0", "periodicity") %in% names(reread)))
})

test_that("trk_pitch_crepe function attributes are correct", {
  expect_equal(attr(trk_pitch_crepe, "ext"),        "crp")
  expect_equal(attr(trk_pitch_crepe, "tracks"),     c("f0", "periodicity"))
  expect_equal(attr(trk_pitch_crepe, "outputType"), "SSFF")
})

test_that("trk_pitch_crepe respects minF/maxF range clamping", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  skip_if(!superassp:::ort_available_cpp(), "ONNX Runtime not available")
  skip_on_cran()
  skip_if_not_installed("huggingfaceR")
  skip_if_offline()

  result <- trk_pitch_crepe(test_wav, model = "tiny", minF = 100, maxF = 200,
                             toFile = FALSE, verbose = FALSE)
  voiced <- result$f0[result$f0 > 0]
  if (length(voiced) > 0) {
    expect_true(all(voiced >= 100 & voiced <= 200))
  }
})

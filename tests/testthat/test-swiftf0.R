test_that("trk_pitch_swiftf0 returns valid AsspDataObj for sustained /a/", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  skip_if(!superassp:::ort_available_cpp(), "ONNX Runtime not available")
  skip_on_cran()
  skip_if_not_installed("huggingfaceR")
  skip_if_offline()

  result <- trk_pitch_swiftf0(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true(all(c("f0", "confidence") %in% names(result)))
  expect_true(nrow(result$f0) > 0L)

  voiced <- result$f0[result$f0 > 0]
  expect_true(length(voiced) > 0)
  expect_gt(median(voiced), 75)
  expect_lt(median(voiced), 400)

  expect_true(all(result$confidence >= 0, na.rm = TRUE))
  # Allow float32 rounding noise at the ceiling (model output, not our code)
  expect_true(all(result$confidence <= 1 + 1e-6, na.rm = TRUE))

  # Fixed frame rate 62.5 Hz
  expect_equal(attr(result, "sampleRate"), 62.5)
})

test_that("trk_pitch_swiftf0 writes SSFF file with .sf0 extension", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  skip_if(!superassp:::ort_available_cpp(), "ONNX Runtime not available")
  skip_on_cran()
  skip_if_not_installed("huggingfaceR")
  skip_if_offline()

  out_dir <- tempdir()
  n <- trk_pitch_swiftf0(test_wav, toFile = TRUE,
                          outputDirectory = out_dir, verbose = FALSE)
  expect_equal(n, 1L)
  out_file <- file.path(out_dir, "a1.sf0")
  expect_true(file.exists(out_file))

  reread <- read_ssff(out_file)
  expect_true(all(c("f0", "confidence") %in% names(reread)))
})

test_that("trk_pitch_swiftf0 function attributes are correct", {
  expect_equal(attr(trk_pitch_swiftf0, "ext"),        "sf0")
  expect_equal(attr(trk_pitch_swiftf0, "tracks"),     c("f0", "confidence"))
  expect_equal(attr(trk_pitch_swiftf0, "outputType"), "SSFF")
})

test_that("trk_pitch_swiftf0 respects minF/maxF voicing range", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  skip_if(!superassp:::ort_available_cpp(), "ONNX Runtime not available")
  skip_on_cran()
  skip_if_not_installed("huggingfaceR")
  skip_if_offline()

  result <- trk_pitch_swiftf0(test_wav, minF = 100, maxF = 200,
                               toFile = FALSE, verbose = FALSE)
  voiced <- result$f0[result$f0 > 0]
  if (length(voiced) > 0) {
    expect_true(all(voiced >= 100 & voiced <= 200))
  }
})

test_that("trk_pitch_swiftf0 rejects out-of-range minF/maxF", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  expect_error(trk_pitch_swiftf0(test_wav, minF = 10, toFile = FALSE, verbose = FALSE))
  expect_error(trk_pitch_swiftf0(test_wav, maxF = 3000, toFile = FALSE, verbose = FALSE))
  expect_error(trk_pitch_swiftf0(test_wav, minF = 300, maxF = 100, toFile = FALSE, verbose = FALSE))
})

test_that("trk_pitch_swiftf0 handles time windowing", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  skip_if(!superassp:::ort_available_cpp(), "ONNX Runtime not available")
  skip_on_cran()
  skip_if_not_installed("huggingfaceR")
  skip_if_offline()

  full <- trk_pitch_swiftf0(test_wav, beginTime = 0, endTime = 0,
                             toFile = FALSE, verbose = FALSE)
  win  <- trk_pitch_swiftf0(test_wav, beginTime = 0, endTime = 0.2,
                             toFile = FALSE, verbose = FALSE)
  expect_gt(nrow(full$f0), nrow(win$f0))
})

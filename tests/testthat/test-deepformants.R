test_that("trk_formant_deepformants returns valid AsspDataObj for sustained /a/", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  skip_if(!superassp:::ort_available_cpp(), "ONNX Runtime not available")

  result <- trk_formant_deepformants(test_wav, numFormants = 3L,
                                      toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("fm" %in% names(result))
  expect_false("bw" %in% names(result))
  expect_equal(ncol(result$fm), 3L)
  expect_true(nrow(result$fm) > 0L)

  # Sustained /a/ vowel: F1 ~ 700–900 Hz
  f1_vals <- result$fm[, 1]
  f2_vals <- result$fm[, 2]
  expect_gt(median(f1_vals), 400)
  expect_lt(median(f1_vals), 1500)

  # Formant ordering: F1 < F2 < F3 (medians)
  expect_lt(median(result$fm[, 1]), median(result$fm[, 2]))
  expect_lt(median(result$fm[, 2]), median(result$fm[, 3]))

  # Frame rate 100 Hz
  expect_equal(attr(result, "sampleRate"), 100.0)
})

test_that("trk_formant_deepformants writes SSFF file with .dff extension", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  skip_if(!superassp:::ort_available_cpp(), "ONNX Runtime not available")

  out_dir <- tempdir()
  n <- trk_formant_deepformants(test_wav, toFile = TRUE,
                                 outputDirectory = out_dir, verbose = FALSE)
  expect_equal(n, 1L)
  out_file <- file.path(out_dir, "a1.dff")
  expect_true(file.exists(out_file))

  reread <- wrassp::read.AsspDataObj(out_file)
  expect_true("fm" %in% names(reread))
  expect_false("bw" %in% names(reread))
})

test_that("trk_formant_deepformants function attributes are correct", {
  expect_equal(attr(trk_formant_deepformants, "ext"),        "dff")
  expect_equal(attr(trk_formant_deepformants, "tracks"),     "fm")
  expect_equal(attr(trk_formant_deepformants, "outputType"), "SSFF")
})

test_that("trk_formant_deepformants respects numFormants", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  skip_if(!superassp:::ort_available_cpp(), "ONNX Runtime not available")

  r1 <- trk_formant_deepformants(test_wav, numFormants = 1L,
                                  toFile = FALSE, verbose = FALSE)
  r4 <- trk_formant_deepformants(test_wav, numFormants = 4L,
                                  toFile = FALSE, verbose = FALSE)
  expect_equal(ncol(r1$fm), 1L)
  expect_equal(ncol(r4$fm), 4L)
})

test_that("trk_formant_deepformants handles time windowing", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  skip_if(!superassp:::ort_available_cpp(), "ONNX Runtime not available")

  full <- trk_formant_deepformants(test_wav, beginTime = 0, endTime = 0,
                                    toFile = FALSE, verbose = FALSE)
  win  <- trk_formant_deepformants(test_wav, beginTime = 0, endTime = 0.2,
                                    toFile = FALSE, verbose = FALSE)
  expect_gt(nrow(full$fm), nrow(win$fm))
})

test_that("trk_formant_deepformants respects windowShift", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  skip_if(!superassp:::ort_available_cpp(), "ONNX Runtime not available")

  r10 <- trk_formant_deepformants(test_wav, windowShift = 10.0,
                                   toFile = FALSE, verbose = FALSE)
  r5  <- trk_formant_deepformants(test_wav, windowShift = 5.0,
                                   toFile = FALSE, verbose = FALSE)
  expect_equal(attr(r10, "sampleRate"), 100.0)
  expect_equal(attr(r5,  "sampleRate"), 200.0)
  expect_gt(nrow(r5$fm), nrow(r10$fm))
})

test_that("covarep module availability check works", {
  # Should return logical
  result <- covarep_available()
  expect_type(result, "logical")
})

test_that("covarep_info returns proper structure", {
  skip_if_not(covarep_available(), "COVAREP not installed")

  info <- covarep_info()

  expect_type(info, "list")
  expect_named(info, c("available", "numba", "optimization_level", "performance"))
  expect_type(info$available, "logical")
  expect_type(info$numba, "logical")
  expect_type(info$optimization_level, "character")
  expect_type(info$performance, "character")
  expect_true(info$available)
})

test_that("trk_covarep_srh works with single file", {
  skip_if_not(covarep_available(), "COVAREP not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_covarep_srh(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("F0[Hz]" %in% names(result))
  expect_true("VUV" %in% names(result))
  expect_true("SRH" %in% names(result))

  # Check dimensions
  expect_true(nrow(result$`F0[Hz]`) > 0)
  expect_equal(ncol(result$`F0[Hz]`), 1)

  # Check attributes
  expect_true(!is.null(attr(result, "sampleRate")))
  expect_true(!is.null(attr(result, "startTime")))
  expect_equal(attr(result, "startTime"), 0.0)
})

test_that("trk_covarep_srh respects F0 range parameters", {
  skip_if_not(covarep_available(), "COVAREP not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_covarep_srh(test_wav,
                            minF = 75,
                            maxF = 300,
                            toFile = FALSE,
                            verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")

  # F0 values should be within range (for voiced frames)
  f0_values <- result$`F0[Hz]`[, 1]
  vuv_values <- result$VUV[, 1]
  voiced_f0 <- f0_values[vuv_values == 1]

  if (length(voiced_f0) > 0) {
    expect_true(all(voiced_f0 >= 75 & voiced_f0 <= 300))
  }
})

test_that("trk_covarep_srh works with custom window shift", {
  skip_if_not(covarep_available(), "COVAREP not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result_5ms <- trk_covarep_srh(test_wav, windowShift = 5.0,
                                toFile = FALSE, verbose = FALSE)
  result_10ms <- trk_covarep_srh(test_wav, windowShift = 10.0,
                                 toFile = FALSE, verbose = FALSE)

  expect_s3_class(result_5ms, "AsspDataObj")
  expect_s3_class(result_10ms, "AsspDataObj")

  # 10ms shift should give approximately half the frames
  expect_true(nrow(result_10ms$`F0[Hz]`) < nrow(result_5ms$`F0[Hz]`))
})

test_that("trk_covarep_srh supports time windowing", {
  skip_if_not(covarep_available(), "COVAREP not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result_full <- trk_covarep_srh(test_wav, toFile = FALSE, verbose = FALSE)
  result_window <- trk_covarep_srh(test_wav,
                                   beginTime = 0.5,
                                   endTime = 1.5,
                                   toFile = FALSE,
                                   verbose = FALSE)

  expect_s3_class(result_window, "AsspDataObj")
  expect_true(nrow(result_window$`F0[Hz]`) < nrow(result_full$`F0[Hz]`))
  expect_equal(attr(result_window, "startTime"), 0.5)
})

test_that("trk_covarep_srh handles batch processing", {
  skip_if_not(covarep_available(), "COVAREP not installed")

  test_files <- list.files(
    system.file("samples", "sustained", package = "superassp"),
    pattern = "\\.wav$",
    full.names = TRUE
  )

  skip_if(length(test_files) < 2, "Need at least 2 test files")
  test_files <- test_files[1:2]

  results <- trk_covarep_srh(test_files, toFile = FALSE, verbose = FALSE)

  expect_type(results, "list")
  expect_length(results, 2)
  expect_s3_class(results[[1]], "AsspDataObj")
  expect_s3_class(results[[2]], "AsspDataObj")
})

test_that("trk_covarep_srh writes files when toFile=TRUE", {
  skip_if_not(covarep_available(), "COVAREP not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  tmp_dir <- tempdir()
  result <- trk_covarep_srh(test_wav,
                            toFile = TRUE,
                            outputDirectory = tmp_dir,
                            verbose = FALSE)

  expect_type(result, "character")
  expect_true(file.exists(result))
  expect_match(result, "\\.f0$")

  # Clean up
  unlink(result)
})

test_that("trk_covarep_srh handles non-WAV formats", {
  skip_if_not(covarep_available(), "COVAREP not installed")

  # Find MP3 test file if available
  test_mp3 <- system.file("samples", "test.mp3", package = "superassp")

  skip_if(test_mp3 == "", "MP3 test file not available")

  result <- trk_covarep_srh(test_mp3, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("F0[Hz]" %in% names(result))
})

test_that("trk_covarep_srh handles missing files gracefully", {
  skip_if_not(covarep_available(), "COVAREP not installed")

  expect_warning(
    result <- trk_covarep_srh("nonexistent.wav", toFile = FALSE, verbose = FALSE),
    "File not found"
  )

  expect_null(result)
})

test_that("trk_covarep_srh validates parameters", {
  skip_if_not(covarep_available(), "COVAREP not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Invalid window shift
  expect_error(
    trk_covarep_srh(test_wav, windowShift = -5, toFile = FALSE, verbose = FALSE),
    "windowShift must be positive"
  )

  # Invalid F0 range
  expect_error(
    trk_covarep_srh(test_wav, minF = 500, maxF = 100, toFile = FALSE, verbose = FALSE),
    "F0 range invalid"
  )
})

test_that("trk_covarep_srh provides consistent results", {
  skip_if_not(covarep_available(), "COVAREP not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result1 <- trk_covarep_srh(test_wav, toFile = FALSE, verbose = FALSE)
  result2 <- trk_covarep_srh(test_wav, toFile = FALSE, verbose = FALSE)

  expect_equal(nrow(result1$`F0[Hz]`), nrow(result2$`F0[Hz]`))
  expect_equal(result1$`F0[Hz]`, result2$`F0[Hz]`)
  expect_equal(result1$VUV, result2$VUV)
})

test_that("trk_covarep_srh function attributes are set", {
  expect_equal(attr(trk_covarep_srh, "ext"), "f0")
  expect_equal(attr(trk_covarep_srh, "tracks"), c("F0[Hz]", "VUV", "SRH"))
  expect_equal(attr(trk_covarep_srh, "outputType"), "SSFF")
})

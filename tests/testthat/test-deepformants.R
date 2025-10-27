test_that("DeepFormants availability check works", {
  # Test deepformants_available function
  is_available <- deepformants_available()
  expect_type(is_available, "logical")
  expect_length(is_available, 1)
})

test_that("trk_deepformants requires Python modules", {
  skip_if(deepformants_available(), "DeepFormants dependencies are installed, skipping error test")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  expect_error(
    trk_deepformants(test_wav, toFile = FALSE, verbose = FALSE),
    "torch|numpy|scipy|pandas|numba"
  )
})

test_that("trk_deepformants works with single file", {
  skip_if_not(deepformants_available(), "DeepFormants dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_deepformants(test_wav, toFile = FALSE, verbose = FALSE)

  # Check result structure
  expect_s3_class(result, "AsspDataObj")
  expect_true("F1" %in% names(result))
  expect_true("F2" %in% names(result))
  expect_true("F3" %in% names(result))
  expect_true("F4" %in% names(result))

  # Check formant tracks are numeric
  expect_type(result$F1, "double")
  expect_type(result$F2, "double")
  expect_type(result$F3, "double")
  expect_type(result$F4, "double")

  # Check all tracks have same length
  expect_equal(length(result$F1), length(result$F2))
  expect_equal(length(result$F2), length(result$F3))
  expect_equal(length(result$F3), length(result$F4))

  # Check attributes
  expect_true(!is.null(attr(result, "sampleRate")))
  expect_true(!is.null(attr(result, "startTime")))
  expect_true(!is.null(attr(result, "origFreq")))

  # Check sample rate is approximately 100 Hz (10ms frames)
  sample_rate <- attr(result, "sampleRate")
  expect_equal(sample_rate, 100, tolerance = 1)
})

test_that("trk_deepformants time windowing", {
  skip_if_not(deepformants_available(), "DeepFormants dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get file info
  info <- av::av_media_info(test_wav)
  duration <- as.numeric(info$duration)

  # Process subset
  if (duration > 1.0) {
    result <- trk_deepformants(test_wav,
                              beginTime = 0.5,
                              endTime = 1.5,
                              toFile = FALSE,
                              verbose = FALSE)

    expect_s3_class(result, "AsspDataObj")
    expect_true("F1" %in% names(result))

    # Check duration is approximately 1 second
    n_frames <- length(result$F1)
    sample_rate <- attr(result, "sampleRate")
    duration_processed <- n_frames / sample_rate
    expect_true(duration_processed >= 0.8 && duration_processed <= 1.2)
  }
})

test_that("trk_deepformants batch processing", {
  skip_if_not(deepformants_available(), "DeepFormants dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Create temporary directory
  temp_dir <- tempdir()

  # Process multiple files (same file twice for testing)
  files <- rep(test_wav, 2)

  n_processed <- trk_deepformants(files,
                                 toFile = TRUE,
                                 outputDirectory = temp_dir,
                                 verbose = FALSE)

  expect_equal(n_processed, 2)

  # Check output files exist
  output_file <- file.path(temp_dir, sub("wav$", "dfm", basename(test_wav)))
  expect_true(file.exists(output_file))

  # Clean up
  unlink(output_file)
})

test_that("trk_deepformants toFile modes", {
  skip_if_not(deepformants_available(), "DeepFormants dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # toFile = FALSE returns AsspDataObj
  result_obj <- trk_deepformants(test_wav, toFile = FALSE, verbose = FALSE)
  expect_s3_class(result_obj, "AsspDataObj")

  # toFile = TRUE returns count
  temp_dir <- tempdir()
  result_count <- trk_deepformants(test_wav,
                                  toFile = TRUE,
                                  outputDirectory = temp_dir,
                                  verbose = FALSE)
  expect_equal(result_count, 1)

  # Check file was written
  output_file <- file.path(temp_dir, sub("wav$", "dfm", basename(test_wav)))
  expect_true(file.exists(output_file))

  # Clean up
  unlink(output_file)
})

test_that("trk_deepformants validates inputs", {
  skip_if_not(deepformants_available(), "DeepFormants dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Multiple files with toFile=FALSE should error
  expect_error(
    trk_deepformants(c(test_wav, test_wav), toFile = FALSE, verbose = FALSE),
    "toFile=FALSE"
  )

  # Missing file should error
  expect_error(
    trk_deepformants("nonexistent.wav", toFile = FALSE, verbose = FALSE),
    "Unable to find"
  )
})

test_that("trk_deepformants output format", {
  skip_if_not(deepformants_available(), "DeepFormants dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_deepformants(test_wav, toFile = FALSE, verbose = FALSE)

  # Check track formats match specification
  track_formats <- attr(result, "trackFormats")
  expect_equal(track_formats, c("REAL32", "REAL32", "REAL32", "REAL32"))

  # Check SSFF file format
  expect_equal(AsspFileFormat(result), "SSFF")
  expect_equal(AsspDataFormat(result), 2)  # binary
})

test_that("trk_deepformants formant values are reasonable", {
  skip_if_not(deepformants_available(), "DeepFormants dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_deepformants(test_wav, toFile = FALSE, verbose = FALSE)

  # Check that formant values are in reasonable ranges
  # F1: typically 200-1000 Hz
  # F2: typically 800-2500 Hz
  # F3: typically 2000-3500 Hz
  # F4: typically 3000-4500 Hz

  # Allow some outliers but most should be in range
  f1_in_range <- sum(result$F1 >= 100 & result$F1 <= 1200, na.rm = TRUE)
  f2_in_range <- sum(result$F2 >= 600 & result$F2 <= 3000, na.rm = TRUE)
  f3_in_range <- sum(result$F3 >= 1500 & result$F3 <= 4000, na.rm = TRUE)
  f4_in_range <- sum(result$F4 >= 2500 & result$F4 <= 5000, na.rm = TRUE)

  total_frames <- length(result$F1)

  # At least 70% should be in reasonable range
  expect_true(f1_in_range / total_frames > 0.7)
  expect_true(f2_in_range / total_frames > 0.7)
  expect_true(f3_in_range / total_frames > 0.7)
  expect_true(f4_in_range / total_frames > 0.7)
})

# ===== lst_deepformants tests =====

test_that("lst_deepformants requires time windows", {
  skip_if_not(deepformants_available(), "DeepFormants dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Missing beginTime/endTime should error
  expect_error(
    lst_deepformants(test_wav, verbose = FALSE),
    "beginTime and endTime"
  )
})

test_that("lst_deepformants works with single file", {
  skip_if_not(deepformants_available(), "DeepFormants dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_deepformants(test_wav,
                             beginTime = 0.5,
                             endTime = 0.6,
                             verbose = FALSE)

  # Check result structure
  expect_type(result, "list")
  expect_true("F1" %in% names(result))
  expect_true("F2" %in% names(result))
  expect_true("F3" %in% names(result))
  expect_true("F4" %in% names(result))
  expect_true("file" %in% names(result))
  expect_true("beginTime" %in% names(result))
  expect_true("endTime" %in% names(result))

  # Check formants are single numeric values
  expect_type(result$F1, "double")
  expect_length(result$F1, 1)
  expect_type(result$F2, "double")
  expect_length(result$F2, 1)
  expect_type(result$F3, "double")
  expect_length(result$F3, 1)
  expect_type(result$F4, "double")
  expect_length(result$F4, 1)
})

test_that("lst_deepformants formant values are reasonable", {
  skip_if_not(deepformants_available(), "DeepFormants dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_deepformants(test_wav,
                             beginTime = 0.5,
                             endTime = 0.6,
                             verbose = FALSE)

  # Check formant values are in reasonable ranges
  expect_true(result$F1 >= 100 && result$F1 <= 1200)
  expect_true(result$F2 >= 600 && result$F2 <= 3000)
  expect_true(result$F3 >= 1500 && result$F3 <= 4000)
  expect_true(result$F4 >= 2500 && result$F4 <= 5000)
})

test_that("lst_deepformants batch processing", {
  skip_if_not(deepformants_available(), "DeepFormants dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Process same file with different time windows
  files <- rep(test_wav, 3)
  begins <- c(0.3, 0.5, 0.7)
  ends <- c(0.4, 0.6, 0.8)

  results <- lst_deepformants(files, begins, ends, verbose = FALSE)

  # Should return list of lists
  expect_type(results, "list")
  expect_length(results, 3)

  # Each result should have formant values
  for (i in 1:3) {
    expect_true("F1" %in% names(results[[i]]))
    expect_true("F2" %in% names(results[[i]]))
    expect_true("F3" %in% names(results[[i]]))
    expect_true("F4" %in% names(results[[i]]))
  }
})

test_that("deepformants_info returns module information", {
  skip_if_not(deepformants_available(), "DeepFormants dependencies not installed")

  info <- deepformants_info()

  expect_type(info, "list")
  expect_true("python_path" %in% names(info))
  expect_true("python_version" %in% names(info))
  expect_true("torch_version" %in% names(info))
  expect_true("numpy_version" %in% names(info))
  expect_true("scipy_version" %in% names(info))
  expect_true("pandas_version" %in% names(info))
  expect_true("numba_version" %in% names(info))
  expect_true("deepformants_path" %in% names(info))
  expect_true("model_files" %in% names(info))

  # Check that model files are present
  expect_true(length(info$model_files) > 0)
})

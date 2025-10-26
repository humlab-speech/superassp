test_that("SAcC availability check works", {
  # Test sacc_available function
  is_available <- sacc_available()
  expect_type(is_available, "logical")
  expect_length(is_available, 1)
})

test_that("trk_sacc requires Python modules", {
  skip_if(sacc_available(), "SAcC dependencies are installed, skipping error test")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  expect_error(
    trk_sacc(test_wav, toFile = FALSE, verbose = FALSE),
    "numpy|scipy|soundfile"
  )
})

test_that("trk_sacc works with single file", {
  skip_if_not(sacc_available(), "SAcC dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_sacc(test_wav, toFile = FALSE, verbose = FALSE)

  # Check result structure
  expect_s3_class(result, "AsspDataObj")
  expect_true("F0" %in% names(result))
  expect_true("prob_voiced" %in% names(result))

  # Check F0 track
  expect_type(result$F0, "double")
  expect_true(length(result$F0) > 0)

  # Check prob_voiced track
  expect_type(result$prob_voiced, "double")
  expect_true(all(result$prob_voiced >= 0 & result$prob_voiced <= 1))

  # Check attributes
  expect_true(!is.null(attr(result, "sampleRate")))
  expect_true(!is.null(attr(result, "startTime")))
  expect_true(!is.null(attr(result, "origFreq")))

  # Check sample rate is approximately 100 Hz (10ms frames)
  sample_rate <- attr(result, "sampleRate")
  expect_equal(sample_rate, 100, tolerance = 1)
})

test_that("trk_sacc custom hmm_vp parameter", {
  skip_if_not(sacc_available(), "SAcC dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # High voicing penalty (more conservative)
  result_high <- trk_sacc(test_wav,
                         hmm_vp = 0.95,
                         toFile = FALSE,
                         verbose = FALSE)

  # Low voicing penalty (more permissive)
  result_low <- trk_sacc(test_wav,
                        hmm_vp = 0.7,
                        toFile = FALSE,
                        verbose = FALSE)

  expect_s3_class(result_high, "AsspDataObj")
  expect_s3_class(result_low, "AsspDataObj")

  # Both should have same number of frames
  expect_equal(length(result_high$F0), length(result_low$F0))
})

test_that("trk_sacc dithering parameter", {
  skip_if_not(sacc_available(), "SAcC dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # With dithering (default)
  result_dither <- trk_sacc(test_wav,
                           dither_level = 1e-3,
                           toFile = FALSE,
                           verbose = FALSE)

  # Without dithering
  result_no_dither <- trk_sacc(test_wav,
                              dither_level = 0,
                              toFile = FALSE,
                              verbose = FALSE)

  expect_s3_class(result_dither, "AsspDataObj")
  expect_s3_class(result_no_dither, "AsspDataObj")
})

test_that("trk_sacc time windowing", {
  skip_if_not(sacc_available(), "SAcC dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get file info
  info <- av::av_media_info(test_wav)
  duration <- as.numeric(info$duration)

  # Process subset
  if (duration > 1.0) {
    result <- trk_sacc(test_wav,
                      beginTime = 0.5,
                      endTime = 1.5,
                      toFile = FALSE,
                      verbose = FALSE)

    expect_s3_class(result, "AsspDataObj")
    expect_true("F0" %in% names(result))

    # Check duration is approximately 1 second
    n_frames <- length(result$F0)
    sample_rate <- attr(result, "sampleRate")
    duration_processed <- n_frames / sample_rate
    expect_true(duration_processed >= 0.8 && duration_processed <= 1.2)
  }
})

test_that("trk_sacc batch processing", {
  skip_if_not(sacc_available(), "SAcC dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Create temporary directory
  temp_dir <- tempdir()

  # Process multiple files (same file twice for testing)
  files <- rep(test_wav, 2)

  n_processed <- trk_sacc(files,
                         toFile = TRUE,
                         outputDirectory = temp_dir,
                         verbose = FALSE)

  expect_equal(n_processed, 2)

  # Check output files exist
  output_file <- file.path(temp_dir, sub("wav$", "sacc", basename(test_wav)))
  expect_true(file.exists(output_file))

  # Clean up
  unlink(output_file)
})

test_that("trk_sacc toFile modes", {
  skip_if_not(sacc_available(), "SAcC dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # toFile = FALSE returns AsspDataObj
  result_obj <- trk_sacc(test_wav, toFile = FALSE, verbose = FALSE)
  expect_s3_class(result_obj, "AsspDataObj")

  # toFile = TRUE returns count
  temp_dir <- tempdir()
  result_count <- trk_sacc(test_wav,
                          toFile = TRUE,
                          outputDirectory = temp_dir,
                          verbose = FALSE)
  expect_equal(result_count, 1)

  # Check file was written
  output_file <- file.path(temp_dir, sub("wav$", "sacc", basename(test_wav)))
  expect_true(file.exists(output_file))

  # Clean up
  unlink(output_file)
})

test_that("trk_sacc validates inputs", {
  skip_if_not(sacc_available(), "SAcC dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Multiple files with toFile=FALSE should error
  expect_error(
    trk_sacc(c(test_wav, test_wav), toFile = FALSE, verbose = FALSE),
    "toFile=FALSE"
  )

  # Missing file should error
  expect_error(
    trk_sacc("nonexistent.wav", toFile = FALSE, verbose = FALSE),
    "Unable to find"
  )
})

test_that("trk_sacc handles non-WAV formats", {
  skip_if_not(sacc_available(), "SAcC dependencies not installed")

  # Find any audio file in test samples
  test_samples <- system.file("samples", package = "superassp")
  audio_files <- list.files(test_samples, pattern = "\\.(wav|mp3|flac|ogg)$",
                           recursive = TRUE, full.names = TRUE)

  skip_if(length(audio_files) == 0, "No audio files found")

  test_file <- audio_files[1]

  result <- trk_sacc(test_file, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("F0" %in% names(result))
})

test_that("trk_sacc output format matches other pitch trackers", {
  skip_if_not(sacc_available(), "SAcC dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_sacc(test_wav, toFile = FALSE, verbose = FALSE)

  # Check track formats match specification
  track_formats <- attr(result, "trackFormats")
  expect_equal(track_formats, c("REAL32", "REAL32"))

  # Check SSFF file format
  expect_equal(AsspFileFormat(result), "SSFF")
  expect_equal(AsspDataFormat(result), 2)  # binary
})

test_that("sacc_info returns module information", {
  skip_if_not(sacc_available(), "SAcC dependencies not installed")

  info <- sacc_info()

  expect_type(info, "list")
  expect_true("python_path" %in% names(info))
  expect_true("python_version" %in% names(info))
  expect_true("numpy_version" %in% names(info))
  expect_true("scipy_version" %in% names(info))
  expect_true("sacc_module_path" %in% names(info))
  expect_true("config_files" %in% names(info))

  # Check that config files are present
  expect_true(length(info$config_files) > 0)
})

test_that("trk_sacc processes audio at 8kHz", {
  skip_if_not(sacc_available(), "SAcC dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_sacc(test_wav, toFile = FALSE, verbose = FALSE)

  # SAcC processes at 8kHz
  orig_freq <- attr(result, "origFreq")
  expect_equal(orig_freq, 8000, tolerance = 1)
})

test_that("trk_sacc F0 values are reasonable for speech", {
  skip_if_not(sacc_available(), "SAcC dependencies not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_sacc(test_wav, toFile = FALSE, verbose = FALSE)

  # Check that voiced F0 values are in typical speech range (80-500 Hz for SAcC)
  voiced_f0 <- result$F0[result$F0 > 0]
  if (length(voiced_f0) > 0) {
    expect_true(all(voiced_f0 >= 60))   # Allow some margin
    expect_true(all(voiced_f0 <= 600))  # Allow some margin
  }
})

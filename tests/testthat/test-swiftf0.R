test_that("Swift-F0 availability check works", {
  # Test swiftf0_available function
  is_available <- swiftf0_available()
  expect_type(is_available, "logical")
  expect_length(is_available, 1)
})

test_that("trk_swiftf0 requires swift-f0 module", {
  skip_if(swiftf0_available(), "Swift-F0 is installed, skipping error test")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  expect_error(
    trk_swiftf0(test_wav, toFile = FALSE, verbose = FALSE),
    "swift-f0 Python module"
  )
})

test_that("trk_swiftf0 works with single file", {
  skip_if_not(swiftf0_available(), "Swift-F0 not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_swiftf0(test_wav, toFile = FALSE, verbose = FALSE)

  # Check result structure
  expect_s3_class(result, "AsspDataObj")
  expect_true("F0" %in% names(result))
  expect_true("confidence" %in% names(result))
  expect_true("voicing" %in% names(result))

  # Check F0 track
  expect_type(result$F0, "integer")
  expect_true(length(result$F0) > 0)

  # Check confidence track
  expect_type(result$confidence, "double")
  expect_true(all(result$confidence >= 0 & result$confidence <= 1))

  # Check voicing track
  expect_type(result$voicing, "integer")
  expect_true(all(result$voicing %in% c(0, 1)))

  # Check attributes
  expect_true(!is.null(attr(result, "sampleRate")))
  expect_true(!is.null(attr(result, "startTime")))
  expect_true(!is.null(attr(result, "origFreq")))
})

test_that("trk_swiftf0 custom F0 range (speech)", {
  skip_if_not(swiftf0_available(), "Swift-F0 not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Speech parameters (75-400 Hz)
  result <- trk_swiftf0(test_wav,
                       minF = 75,
                       maxF = 400,
                       toFile = FALSE,
                       verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("F0" %in% names(result))

  # Check that voiced frames are within range (allowing for 0 = unvoiced)
  voiced_f0 <- result$F0[result$F0 > 0]
  if (length(voiced_f0) > 0) {
    expect_true(all(voiced_f0 >= 75))
    expect_true(all(voiced_f0 <= 400))
  }
})

test_that("trk_swiftf0 custom F0 range (music)", {
  skip_if_not(swiftf0_available(), "Swift-F0 not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Music parameters (full range)
  result <- trk_swiftf0(test_wav,
                       minF = 46.875,
                       maxF = 2093.75,
                       toFile = FALSE,
                       verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("F0" %in% names(result))
})

test_that("trk_swiftf0 confidence threshold parameter", {
  skip_if_not(swiftf0_available(), "Swift-F0 not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # High threshold (more conservative voicing)
  result_high <- trk_swiftf0(test_wav,
                            confidence_threshold = 0.95,
                            toFile = FALSE,
                            verbose = FALSE)

  # Low threshold (more permissive voicing)
  result_low <- trk_swiftf0(test_wav,
                           confidence_threshold = 0.7,
                           toFile = FALSE,
                           verbose = FALSE)

  expect_s3_class(result_high, "AsspDataObj")
  expect_s3_class(result_low, "AsspDataObj")

  # Lower threshold should generally result in more voiced frames
  # (but this depends on the signal, so we just check both work)
  expect_true(length(result_high$F0) == length(result_low$F0))
})

test_that("trk_swiftf0 time windowing", {
  skip_if_not(swiftf0_available(), "Swift-F0 not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get file info
  info <- av::av_media_info(test_wav)
  duration <- as.numeric(info$duration)

  # Process subset
  if (duration > 1.0) {
    result <- trk_swiftf0(test_wav,
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

test_that("trk_swiftf0 batch processing", {
  skip_if_not(swiftf0_available(), "Swift-F0 not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Create temporary directory
  temp_dir <- tempdir()

  # Process multiple files (same file twice for testing)
  files <- rep(test_wav, 2)

  n_processed <- trk_swiftf0(files,
                            toFile = TRUE,
                            outputDirectory = temp_dir,
                            verbose = FALSE)

  expect_equal(n_processed, 2)

  # Check output files exist
  output_file <- file.path(temp_dir, sub("wav$", "sf0", basename(test_wav)))
  expect_true(file.exists(output_file))

  # Clean up
  unlink(output_file)
})

test_that("trk_swiftf0 toFile modes", {
  skip_if_not(swiftf0_available(), "Swift-F0 not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # toFile = FALSE returns AsspDataObj
  result_obj <- trk_swiftf0(test_wav, toFile = FALSE, verbose = FALSE)
  expect_s3_class(result_obj, "AsspDataObj")

  # toFile = TRUE returns count
  temp_dir <- tempdir()
  result_count <- trk_swiftf0(test_wav,
                             toFile = TRUE,
                             outputDirectory = temp_dir,
                             verbose = FALSE)
  expect_equal(result_count, 1)

  # Check file was written
  output_file <- file.path(temp_dir, sub("wav$", "sf0", basename(test_wav)))
  expect_true(file.exists(output_file))

  # Clean up
  unlink(output_file)
})

test_that("trk_swiftf0 validates inputs", {
  skip_if_not(swiftf0_available(), "Swift-F0 not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Multiple files with toFile=FALSE should error
  expect_error(
    trk_swiftf0(c(test_wav, test_wav), toFile = FALSE, verbose = FALSE),
    "toFile=FALSE"
  )

  # Missing file should error
  expect_error(
    trk_swiftf0("nonexistent.wav", toFile = FALSE, verbose = FALSE),
    "Unable to find"
  )
})

test_that("trk_swiftf0 handles non-WAV formats", {
  skip_if_not(swiftf0_available(), "Swift-F0 not installed")

  # Find any audio file in test samples
  test_samples <- system.file("samples", package = "superassp")
  audio_files <- list.files(test_samples, pattern = "\\.(wav|mp3|flac|ogg)$",
                           recursive = TRUE, full.names = TRUE)

  skip_if(length(audio_files) == 0, "No audio files found")

  test_file <- audio_files[1]

  result <- trk_swiftf0(test_file, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("F0" %in% names(result))
})

test_that("trk_swiftf0 S7 dispatch with AVAudio", {
  skip_if_not(swiftf0_available(), "Swift-F0 not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load audio into AVAudio object
  audio <- read_avaudio(test_wav)

  # Process with AVAudio (automatic S7 dispatch)
  result <- trk_swiftf0(audio, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("F0" %in% names(result))
  expect_true("confidence" %in% names(result))
  expect_true("voicing" %in% names(result))
})

test_that("trk_swiftf0 S7 dispatch consistency", {
  skip_if_not(swiftf0_available(), "Swift-F0 not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Process with file path
  result_path <- trk_swiftf0(test_wav, toFile = FALSE, verbose = FALSE)

  # Process with AVAudio
  audio <- read_avaudio(test_wav)
  result_avaudio <- trk_swiftf0(audio, toFile = FALSE, verbose = FALSE)

  # Results should be similar (may not be identical due to file I/O)
  expect_equal(length(result_path$F0), length(result_avaudio$F0))
  expect_equal(attr(result_path, "sampleRate"), attr(result_avaudio, "sampleRate"),
               tolerance = 0.1)
})

test_that("trk_swiftf0 AVAudio preprocessing", {
  skip_if_not(swiftf0_available(), "Swift-F0 not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load with preprocessing (16kHz mono - optimal for Swift-F0)
  audio <- read_avaudio(test_wav,
                       sample_rate = 16000,
                       channels = 1)

  result <- trk_swiftf0(audio, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("F0" %in% names(result))

  # Check original frequency was 16kHz
  expect_equal(attr(result, "origFreq"), 16000, tolerance = 1)
})

test_that("trk_swiftf0 output format matches other pitch trackers", {
  skip_if_not(swiftf0_available(), "Swift-F0 not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_swiftf0(test_wav, toFile = FALSE, verbose = FALSE)

  # Check track formats match specification
  track_formats <- attr(result, "trackFormats")
  expect_equal(track_formats, c("INT16", "REAL32", "INT16"))

  # Check SSFF file format
  expect_equal(AsspFileFormat(result), "SSFF")
  expect_equal(AsspDataFormat(result), 2)  # binary
})

test_that("swiftf0_info returns module information", {
  skip_if_not(swiftf0_available(), "Swift-F0 not installed")

  info <- swiftf0_info()

  expect_type(info, "list")
  expect_true("version" %in% names(info))
  expect_true("target_sample_rate" %in% names(info))
  expect_true("hop_length" %in% names(info))
  expect_true("model_fmin" %in% names(info))
  expect_true("model_fmax" %in% names(info))

  # Check values are reasonable
  expect_equal(info$target_sample_rate, 16000)
  expect_equal(info$hop_length, 256)
  expect_equal(info$model_fmin, 46.875)
  expect_equal(info$model_fmax, 2093.75)
})

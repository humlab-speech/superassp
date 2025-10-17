# Tests for SPTK C++ Pitch Tracking Functions

test_that("rapt_cpp works with default parameters", {
  skip_if_not_installed("superassp")

  # Get test file
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load audio
  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Run RAPT pitch extraction
  result <- superassp::rapt_cpp(audio_obj)

  # Check result structure
  expect_type(result, "list")
  expect_named(result, c("f0", "times", "sample_rate", "n_frames"))

  # Check data types
  expect_true(is.matrix(result$f0))
  expect_type(result$times, "double")
  expect_type(result$sample_rate, "double")
  expect_type(result$n_frames, "double")

  # Check reasonable values
  expect_true(all(result$times >= 0))
  expect_true(result$sample_rate > 0)
  expect_true(result$n_frames > 0)
  expect_equal(length(result$times), result$n_frames)
  expect_equal(nrow(result$f0), result$n_frames)

  # Check F0 values are in reasonable speech range (or 0 for unvoiced)
  f0_values <- result$f0[result$f0 > 0]
  if (length(f0_values) > 0) {
    expect_true(all(f0_values >= 50))   # At least 50 Hz
    expect_true(all(f0_values <= 500))  # At most 500 Hz
  }
})

test_that("swipe_cpp works with default parameters", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)
  result <- superassp::swipe_cpp(audio_obj)

  # Check result structure
  expect_type(result, "list")
  expect_named(result, c("f0", "times", "sample_rate", "n_frames"))

  # Check data types
  expect_true(is.matrix(result$f0))
  expect_type(result$times, "double")
  expect_type(result$sample_rate, "double")
  expect_type(result$n_frames, "double")

  # Check reasonable values
  expect_true(result$n_frames > 0)
  expect_equal(nrow(result$f0), result$n_frames)

  # Check F0 values
  f0_values <- result$f0[result$f0 > 0]
  if (length(f0_values) > 0) {
    expect_true(all(f0_values >= 50))
    expect_true(all(f0_values <= 500))
  }
})

test_that("reaper_cpp works with default parameters", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)
  result <- superassp::reaper_cpp(audio_obj)

  # Check result structure - REAPER includes additional epoch information
  expect_type(result, "list")
  expect_named(result, c("f0", "times", "sample_rate", "n_frames", "epochs", "n_epochs", "polarity"))

  # Check data types
  expect_true(is.matrix(result$f0))
  expect_type(result$times, "double")
  expect_type(result$sample_rate, "double")
  expect_type(result$n_frames, "double")
  expect_type(result$epochs, "double")
  expect_type(result$n_epochs, "double")
  expect_type(result$polarity, "character")

  # Check reasonable values
  expect_true(result$n_frames > 0)
  expect_equal(nrow(result$f0), result$n_frames)

  # Check epochs
  expect_true(result$n_epochs >= 0)
  expect_equal(length(result$epochs), result$n_epochs)
  if (result$n_epochs > 0) {
    expect_true(all(result$epochs >= 0))
  }

  # Check polarity
  expect_true(result$polarity %in% c("positive", "negative", "unknown"))

  # Check F0 values
  f0_values <- result$f0[result$f0 > 0]
  if (length(f0_values) > 0) {
    expect_true(all(f0_values >= 50))
    expect_true(all(f0_values <= 500))
  }
})

test_that("dio_cpp works with default parameters", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)
  result <- superassp::dio_cpp(audio_obj)

  # Check result structure
  expect_type(result, "list")
  expect_named(result, c("f0", "times", "sample_rate", "n_frames"))

  # Check data types
  expect_true(is.matrix(result$f0))
  expect_type(result$times, "double")
  expect_type(result$sample_rate, "double")
  expect_type(result$n_frames, "double")

  # Check reasonable values
  expect_true(result$n_frames > 0)
  expect_equal(nrow(result$f0), result$n_frames)

  # Check F0 values
  f0_values <- result$f0[result$f0 > 0]
  if (length(f0_values) > 0) {
    expect_true(all(f0_values >= 50))
    expect_true(all(f0_values <= 500))
  }
})

test_that("SPTK C++ functions work with custom F0 range", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Test with female voice range
  result_female <- superassp::rapt_cpp(audio_obj, minF = 100, maxF = 500)
  expect_true(result_female$n_frames > 0)

  # Test with male voice range
  result_male <- superassp::swipe_cpp(audio_obj, minF = 50, maxF = 250)
  expect_true(result_male$n_frames > 0)

  # Test with wide range
  result_wide <- superassp::dio_cpp(audio_obj, minF = 40, maxF = 600)
  expect_true(result_wide$n_frames > 0)
})

test_that("SPTK C++ functions work with custom windowShift", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Test with shorter frame shift (more frames)
  result_5ms <- superassp::rapt_cpp(audio_obj, windowShift = 5)

  # Test with longer frame shift (fewer frames)
  result_20ms <- superassp::rapt_cpp(audio_obj, windowShift = 20)

  # More frequent sampling should produce more frames
  expect_true(result_5ms$n_frames > result_20ms$n_frames)
})

test_that("SPTK C++ functions work with custom voicing threshold", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Test RAPT with different voicing thresholds
  result_low <- superassp::rapt_cpp(audio_obj, voicing_threshold = 0.3)
  result_high <- superassp::rapt_cpp(audio_obj, voicing_threshold = 0.95)

  # Both should produce results
  expect_true(result_low$n_frames > 0)
  expect_true(result_high$n_frames > 0)

  # Test SWIPE with different voicing thresholds
  result_swipe_low <- superassp::swipe_cpp(audio_obj, voicing_threshold = 0.1)
  result_swipe_high <- superassp::swipe_cpp(audio_obj, voicing_threshold = 0.5)

  expect_true(result_swipe_low$n_frames > 0)
  expect_true(result_swipe_high$n_frames > 0)
})

test_that("SPTK C++ functions error handling", {
  skip_if_not_installed("superassp")

  # Test with invalid input
  expect_error(
    superassp::rapt_cpp("not an audio object"),
    "must be an AsspDataObj"
  )

  expect_error(
    superassp::swipe_cpp("not an audio object"),
    "must be an AsspDataObj"
  )

  expect_error(
    superassp::reaper_cpp("not an audio object"),
    "must be an AsspDataObj"
  )

  expect_error(
    superassp::dio_cpp("not an audio object"),
    "must be an AsspDataObj"
  )

  # Test with object missing audio track
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Remove audio track
  audio_obj_bad <- audio_obj
  audio_obj_bad$audio <- NULL

  expect_error(
    superassp::rapt_cpp(audio_obj_bad),
    "must contain 'audio' track"
  )

  expect_error(
    superassp::swipe_cpp(audio_obj_bad),
    "must contain 'audio' track"
  )
})

test_that("SPTK C++ functions produce consistent results", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Run twice with same parameters for each algorithm
  rapt1 <- superassp::rapt_cpp(audio_obj, minF = 60, maxF = 400, windowShift = 10)
  rapt2 <- superassp::rapt_cpp(audio_obj, minF = 60, maxF = 400, windowShift = 10)

  swipe1 <- superassp::swipe_cpp(audio_obj, minF = 60, maxF = 400, windowShift = 10)
  swipe2 <- superassp::swipe_cpp(audio_obj, minF = 60, maxF = 400, windowShift = 10)

  # Results should be identical
  expect_equal(rapt1$n_frames, rapt2$n_frames)
  expect_equal(rapt1$f0, rapt2$f0)
  expect_equal(rapt1$times, rapt2$times)

  expect_equal(swipe1$n_frames, swipe2$n_frames)
  expect_equal(swipe1$f0, swipe2$f0)
  expect_equal(swipe1$times, swipe2$times)
})

test_that("SPTK C++ verbose output works", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Capture output for each algorithm
  output_rapt <- capture.output({
    result <- superassp::rapt_cpp(audio_obj, verbose = TRUE)
  })

  output_swipe <- capture.output({
    result <- superassp::swipe_cpp(audio_obj, verbose = TRUE)
  })

  output_reaper <- capture.output({
    result <- superassp::reaper_cpp(audio_obj, verbose = TRUE)
  })

  output_dio <- capture.output({
    result <- superassp::dio_cpp(audio_obj, verbose = TRUE)
  })

  # Should have some output
  expect_true(length(output_rapt) > 0)
  expect_true(any(grepl("Processing", output_rapt)))

  expect_true(length(output_swipe) > 0)
  expect_true(any(grepl("Processing", output_swipe)))

  expect_true(length(output_reaper) > 0)
  expect_true(any(grepl("Processing", output_reaper)))

  expect_true(length(output_dio) > 0)
  expect_true(any(grepl("Processing", output_dio)))
})

test_that("SPTK C++ functions handle short audio", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load and truncate audio to 0.1 seconds
  audio_obj <- superassp::av_to_asspDataObj(test_wav)
  n_samples <- round(0.1 * attr(audio_obj, "sampleRate"))
  audio_obj$audio <- audio_obj$audio[1:n_samples, , drop = FALSE]

  # Should still work with short audio
  result_rapt <- superassp::rapt_cpp(audio_obj)
  result_swipe <- superassp::swipe_cpp(audio_obj)
  result_reaper <- superassp::reaper_cpp(audio_obj)
  result_dio <- superassp::dio_cpp(audio_obj)

  # All should produce results
  expect_true(result_rapt$n_frames > 0)
  expect_true(result_swipe$n_frames > 0)
  expect_true(result_reaper$n_frames > 0)
  expect_true(result_dio$n_frames > 0)
})

test_that("REAPER epochs are properly ordered", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)
  result <- superassp::reaper_cpp(audio_obj)

  # Check that epochs are monotonically increasing
  if (result$n_epochs > 1) {
    expect_true(all(diff(result$epochs) > 0))
  }

  # Check that epoch times are reasonable (within signal duration)
  if (result$n_epochs > 0) {
    duration <- max(result$times)
    expect_true(all(result$epochs <= duration + 0.1))  # Allow small tolerance
  }
})

# =============================================================================
# Tests for R Wrapper Functions (rapt, swipe, reaper, dio)
# =============================================================================

test_that("rapt() wrapper works with single file", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Test with toFile = FALSE (returns AsspDataObj)
  result <- superassp::rapt(test_wav, toFile = FALSE, verbose = FALSE)

  # Check result structure
  expect_s3_class(result, "AsspDataObj")
  expect_true("f0" %in% names(result))

  # Check result has proper attributes
  expect_true(!is.null(attr(result, "sampleRate")))
  expect_true(!is.null(attr(result, "startTime")))
  expect_true(!is.null(attr(result, "startRecord")))
  expect_true(!is.null(attr(result, "endRecord")))

  # Check F0 values are reasonable
  f0_values <- result$f0[result$f0 > 0]
  if (length(f0_values) > 0) {
    expect_true(all(f0_values >= 50))
    expect_true(all(f0_values <= 500))
  }
})

test_that("swipe() wrapper works with single file", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- superassp::swipe(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("f0" %in% names(result))

  f0_values <- result$f0[result$f0 > 0]
  if (length(f0_values) > 0) {
    expect_true(all(f0_values >= 50))
    expect_true(all(f0_values <= 500))
  }
})

test_that("reaper() wrapper works with single file", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- superassp::reaper(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("f0" %in% names(result))

  # Check that epoch attributes exist
  expect_true(!is.null(attr(result, "epochs")))
  expect_true(!is.null(attr(result, "n_epochs")))
  expect_true(!is.null(attr(result, "polarity")))

  f0_values <- result$f0[result$f0 > 0]
  if (length(f0_values) > 0) {
    expect_true(all(f0_values >= 50))
    expect_true(all(f0_values <= 500))
  }
})

test_that("dio() wrapper works with single file", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- superassp::dio(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("f0" %in% names(result))

  f0_values <- result$f0[result$f0 > 0]
  if (length(f0_values) > 0) {
    expect_true(all(f0_values >= 50))
    expect_true(all(f0_values <= 500))
  }
})

test_that("R wrappers work with custom F0 range", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Test with female voice range
  result_rapt <- superassp::rapt(test_wav, minF = 100, maxF = 500,
                                  toFile = FALSE, verbose = FALSE)
  expect_s3_class(result_rapt, "AsspDataObj")

  result_swipe <- superassp::swipe(test_wav, minF = 100, maxF = 500,
                                    toFile = FALSE, verbose = FALSE)
  expect_s3_class(result_swipe, "AsspDataObj")
})

test_that("R wrappers work with time windowing", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get duration
  info <- av::av_media_info(test_wav)
  duration <- info$duration

  # Test with time window
  result <- superassp::rapt(test_wav, beginTime = 0.1,
                           endTime = min(0.5, duration - 0.1),
                           toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("f0" %in% names(result))
})

test_that("R wrappers work with multiple files", {
  skip_if_not_installed("superassp")

  # Get multiple test files
  test_files <- system.file("samples", "sustained",
                           c("a1.wav", "a32b.wav"),
                           package = "superassp")
  test_files <- test_files[file.exists(test_files)]

  skip_if(length(test_files) < 2, "Not enough test files found")

  # Test rapt with multiple files
  results <- superassp::rapt(test_files, toFile = FALSE, verbose = FALSE)

  expect_type(results, "list")
  expect_equal(length(results), length(test_files))

  for (result in results) {
    expect_s3_class(result, "AsspDataObj")
    expect_true("f0" %in% names(result))
  }
})

test_that("R wrappers can write to file", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Create temp directory
  temp_dir <- tempdir()

  # Test rapt with toFile = TRUE
  n_written <- superassp::rapt(test_wav, toFile = TRUE,
                               outputDirectory = temp_dir,
                               explicitExt = "f0",
                               verbose = FALSE)

  expect_equal(n_written, 1)

  # Check that file was created
  output_file <- file.path(temp_dir, "a1.f0")
  expect_true(file.exists(output_file))

  # Clean up
  unlink(output_file)
})

test_that("R wrappers handle non-WAV files via av package", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("av")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Convert to MP3 for testing
  temp_mp3 <- tempfile(fileext = ".mp3")
  on.exit(unlink(temp_mp3), add = TRUE)

  tryCatch({
    av::av_audio_convert(test_wav, temp_mp3, format = "mp3")
  }, error = function(e) {
    skip("Could not create MP3 test file")
  })

  skip_if(!file.exists(temp_mp3), "MP3 file not created")

  # Test that rapt can process MP3
  result <- superassp::rapt(temp_mp3, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("f0" %in% names(result))
})

test_that("R wrappers have consistent output with C++ functions", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load audio
  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Compare R wrapper with direct C++ call for RAPT
  result_wrapper <- superassp::rapt(test_wav, minF = 60, maxF = 400,
                                    windowShift = 10.0, voicing_threshold = 0.9,
                                    toFile = FALSE, verbose = FALSE)

  result_cpp <- superassp::rapt_cpp(audio_obj, minF = 60, maxF = 400,
                                    windowShift = 10.0, voicing_threshold = 0.9)

  # Both should produce similar number of frames
  n_frames_wrapper <- attr(result_wrapper, "endRecord") - attr(result_wrapper, "startRecord") + 1
  expect_equal(n_frames_wrapper, result_cpp$n_frames, tolerance = 1)

  # F0 values should be comparable
  expect_equal(dim(result_wrapper$f0), dim(result_cpp$f0))
})

test_that("R wrapper error handling works correctly", {
  skip_if_not_installed("superassp")

  # Test with non-existent file
  expect_error(
    superassp::rapt("/nonexistent/file.wav", toFile = FALSE),
    "do not exist"
  )

  # Test with empty file list
  expect_error(
    superassp::rapt(NULL, toFile = FALSE),
    "No input files specified"
  )

  expect_error(
    superassp::swipe(character(0), toFile = FALSE),
    "No input files specified"
  )
})

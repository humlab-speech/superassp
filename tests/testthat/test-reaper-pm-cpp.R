# Tests for REAPER Pitch Mark C++ Implementation (trk_reaper_pm)
# New in v0.9.0 - C++ implementation replacing Python reaper_pm()

test_that("trk_reaper_pm works with default parameters", {
  skip_if_not_installed("superassp")

  # Get test file
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Run REAPER pitch mark extraction with toFile = FALSE
  result <- superassp::trk_reaper_pm(test_wav, toFile = FALSE, verbose = FALSE)

  # Check result structure
  expect_s3_class(result, "AsspDataObj")
  expect_true("pm" %in% names(result))

  # Check data type (should be INT16 matrix)
  expect_true(is.matrix(result$pm))
  expect_type(result$pm, "integer")

  # Check dimensions
  expect_equal(ncol(result$pm), 1)  # Single column
  expect_true(nrow(result$pm) > 0)  # At least some frames

  # Check metadata
  expect_true(!is.null(attr(result, "sampleRate")))
  expect_true(!is.null(attr(result, "startTime")))
  expect_true(!is.null(attr(result, "origFreq")))

  # Check frame rate matches windowShift (default 10ms = 100 Hz)
  expect_equal(attr(result, "sampleRate"), 100)
})

test_that("trk_reaper_pm returns correct binary grid format", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- superassp::trk_reaper_pm(test_wav, toFile = FALSE, verbose = FALSE)

  # Binary grid should contain only 0 or 1
  pm_values <- as.vector(result$pm)
  expect_true(all(pm_values %in% c(0L, 1L)))

  # Should have at least some pitch marks (1s) in voiced speech
  n_pitch_marks <- sum(pm_values == 1L)
  expect_true(n_pitch_marks > 0,
              info = "Expected at least some pitch marks in sustained vowel")
})

test_that("trk_reaper_pm epoch attributes are valid", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- superassp::trk_reaper_pm(test_wav, toFile = FALSE, verbose = FALSE)

  # Check epoch attributes exist
  expect_true(!is.null(attr(result, "epoch_times")))
  expect_true(!is.null(attr(result, "n_epochs")))
  expect_true(!is.null(attr(result, "polarity")))

  # Validate epoch data
  epoch_times <- attr(result, "epoch_times")
  n_epochs <- attr(result, "n_epochs")
  polarity <- attr(result, "polarity")

  # n_epochs should match length of epoch_times
  expect_equal(length(epoch_times), n_epochs)

  # Epoch times should be numeric and non-negative
  expect_type(epoch_times, "double")
  if (n_epochs > 0) {
    expect_true(all(epoch_times >= 0))
  }

  # Polarity should be valid
  expect_true(polarity %in% c("positive", "negative", "unknown"))

  # Epoch times should be monotonically increasing
  if (n_epochs > 1) {
    expect_true(all(diff(epoch_times) > 0))
  }
})

test_that("trk_reaper_pm works with custom F0 range", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Test with female voice range (100-500 Hz)
  result_female <- superassp::trk_reaper_pm(
    test_wav,
    minF = 100,
    maxF = 500,
    toFile = FALSE,
    verbose = FALSE
  )
  expect_s3_class(result_female, "AsspDataObj")
  expect_true("pm" %in% names(result_female))

  # Test with male voice range (50-250 Hz)
  result_male <- superassp::trk_reaper_pm(
    test_wav,
    minF = 50,
    maxF = 250,
    toFile = FALSE,
    verbose = FALSE
  )
  expect_s3_class(result_male, "AsspDataObj")
  expect_true("pm" %in% names(result_male))

  # Test with wide range (40-600 Hz)
  result_wide <- superassp::trk_reaper_pm(
    test_wav,
    minF = 40,
    maxF = 600,
    toFile = FALSE,
    verbose = FALSE
  )
  expect_s3_class(result_wide, "AsspDataObj")
  expect_true("pm" %in% names(result_wide))
})

test_that("trk_reaper_pm works with custom windowShift", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Test with 5ms windowShift (200 Hz)
  result_5ms <- superassp::trk_reaper_pm(
    test_wav,
    windowShift = 5.0,
    toFile = FALSE,
    verbose = FALSE
  )
  expect_equal(attr(result_5ms, "sampleRate"), 200)

  # Test with 20ms windowShift (50 Hz)
  result_20ms <- superassp::trk_reaper_pm(
    test_wav,
    windowShift = 20.0,
    toFile = FALSE,
    verbose = FALSE
  )
  expect_equal(attr(result_20ms, "sampleRate"), 50)

  # Shorter windowShift should produce more frames
  expect_true(nrow(result_5ms$pm) > nrow(result_20ms$pm))
})

test_that("trk_reaper_pm works with custom voicing threshold", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Test with low threshold (more permissive)
  result_low <- superassp::trk_reaper_pm(
    test_wav,
    voicing_threshold = 0.5,
    toFile = FALSE,
    verbose = FALSE
  )
  expect_s3_class(result_low, "AsspDataObj")

  # Test with high threshold (more strict)
  result_high <- superassp::trk_reaper_pm(
    test_wav,
    voicing_threshold = 0.95,
    toFile = FALSE,
    verbose = FALSE
  )
  expect_s3_class(result_high, "AsspDataObj")

  # Both should produce results
  expect_true(nrow(result_low$pm) > 0)
  expect_true(nrow(result_high$pm) > 0)
})

test_that("trk_reaper_pm handles time windowing", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get audio duration
  info <- av::av_media_info(test_wav)
  duration <- as.numeric(info$duration)

  skip_if(duration < 1.0, "Audio file too short for windowing test")

  # Extract full file
  result_full <- superassp::trk_reaper_pm(
    test_wav,
    toFile = FALSE,
    verbose = FALSE
  )

  # Extract middle section (0.2 to 0.8 seconds)
  result_windowed <- superassp::trk_reaper_pm(
    test_wav,
    beginTime = 0.2,
    endTime = 0.8,
    toFile = FALSE,
    verbose = FALSE
  )

  # Windowed result should have fewer frames
  expect_true(nrow(result_windowed$pm) < nrow(result_full$pm))

  # Windowed result should have approximately 60 frames (0.6 sec * 100 Hz)
  # Allow tolerance for edge effects
  expect_gt(nrow(result_windowed$pm), 50)
  expect_lt(nrow(result_windowed$pm), 70)
})

test_that("trk_reaper_pm writes SSFF files correctly", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Create temporary output directory
  temp_dir <- tempdir()

  # Process with toFile = TRUE
  n_written <- superassp::trk_reaper_pm(
    test_wav,
    toFile = TRUE,
    outputDirectory = temp_dir,
    explicitExt = "rpm",
    verbose = FALSE
  )

  # Check return value
  expect_equal(n_written, 1)

  # Check that .rpm file was created
  expected_file <- file.path(temp_dir, "a1.rpm")
  expect_true(file.exists(expected_file))

  # Read back the SSFF file
  if (file.exists(expected_file)) {
    result <- wrassp::read.AsspDataObj(expected_file)
    expect_s3_class(result, "AsspDataObj")
    expect_true("pm" %in% names(result))
    expect_true(is.matrix(result$pm))

    # Verify data type
    expect_type(result$pm, "integer")

    # Verify binary values
    expect_true(all(result$pm %in% c(0L, 1L)))

    # Clean up
    unlink(expected_file)
  }
})

test_that("trk_reaper_pm processes multiple files", {
  skip_if_not_installed("superassp")

  # Get multiple test files
  test_files <- system.file("samples", "sustained",
                           c("a1.wav", "a32b.wav"),
                           package = "superassp")
  test_files <- test_files[file.exists(test_files)]

  skip_if(length(test_files) < 2, "Not enough test files found")

  # Process multiple files with toFile = FALSE
  results <- superassp::trk_reaper_pm(
    test_files,
    toFile = FALSE,
    verbose = FALSE
  )

  # Check results structure
  expect_type(results, "list")
  expect_equal(length(results), length(test_files))

  # Check each result
  for (result in results) {
    expect_s3_class(result, "AsspDataObj")
    expect_true("pm" %in% names(result))
    expect_true(is.matrix(result$pm))
    expect_type(result$pm, "integer")
  }
})

test_that("trk_reaper_pm matches reaper_cpp epochs", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load audio
  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Get epochs from reaper_cpp
  reaper_result <- superassp::reaper_cpp(audio_obj, windowShift = 10.0)

  # Get pitch marks from trk_reaper_pm
  pm_result <- superassp::trk_reaper_pm(
    test_wav,
    windowShift = 10.0,
    toFile = FALSE,
    verbose = FALSE
  )

  # Extract epoch data
  epochs_from_reaper <- reaper_result$epochs
  epochs_from_pm <- attr(pm_result, "epoch_times")

  # Should have similar number of epochs
  expect_equal(
    length(epochs_from_reaper),
    length(epochs_from_pm),
    tolerance = 2,
    info = "trk_reaper_pm should extract same epochs as reaper_cpp"
  )

  # Epoch times should match (they come from same C++ function)
  if (length(epochs_from_reaper) > 0 && length(epochs_from_pm) > 0) {
    min_len <- min(length(epochs_from_reaper), length(epochs_from_pm))
    expect_equal(
      epochs_from_reaper[1:min_len],
      epochs_from_pm[1:min_len],
      tolerance = 0.001  # 1ms tolerance
    )
  }
})

test_that("trk_reaper_pm binary grid conversion is accurate", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get result
  result <- superassp::trk_reaper_pm(
    test_wav,
    windowShift = 10.0,
    toFile = FALSE,
    verbose = FALSE
  )

  # Extract data
  pm_grid <- result$pm
  epoch_times <- attr(result, "epoch_times")
  n_epochs <- attr(result, "n_epochs")
  frame_rate <- attr(result, "sampleRate")  # 100 Hz for 10ms shift

  # Count 1s in binary grid
  n_marks_in_grid <- sum(pm_grid == 1L)

  # Should be approximately equal to n_epochs
  # (may differ slightly due to grid quantization)
  expect_equal(n_marks_in_grid, n_epochs, tolerance = n_epochs * 0.1)

  # Verify frame indices for each epoch
  if (n_epochs > 0) {
    frame_shift_sec <- 1.0 / frame_rate

    for (i in seq_along(epoch_times)) {
      # Calculate expected frame index
      expected_frame <- floor(epoch_times[i] / frame_shift_sec) + 1

      # Frame should be within grid bounds
      expect_true(expected_frame >= 1)
      expect_true(expected_frame <= nrow(pm_grid))
    }
  }
})

test_that("trk_reaper_pm handles short audio files", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load and truncate to 0.1 seconds
  audio_obj <- superassp::av_to_asspDataObj(test_wav)
  n_samples <- round(0.1 * attr(audio_obj, "sampleRate"))
  audio_obj_short <- audio_obj
  audio_obj_short$audio <- audio_obj$audio[1:n_samples, , drop = FALSE]

  # Write short audio to temp file
  temp_wav <- tempfile(fileext = ".wav")
  wrassp::write.AsspDataObj(audio_obj_short, temp_wav)
  on.exit(unlink(temp_wav), add = TRUE)

  # Should still process short audio
  result <- superassp::trk_reaper_pm(
    temp_wav,
    toFile = FALSE,
    verbose = FALSE
  )

  expect_s3_class(result, "AsspDataObj")
  expect_true("pm" %in% names(result))

  # Should have approximately 10 frames (0.1 sec * 100 Hz)
  expect_gt(nrow(result$pm), 5)
  expect_lt(nrow(result$pm), 15)
})

test_that("trk_reaper_pm error handling works", {
  skip_if_not_installed("superassp")

  # Test with non-existent file
  expect_error(
    superassp::trk_reaper_pm("/nonexistent/file.wav", toFile = FALSE),
    "do not exist"
  )

  # Test with empty file list
  expect_error(
    superassp::trk_reaper_pm(NULL, toFile = FALSE),
    "No input files specified"
  )

  expect_error(
    superassp::trk_reaper_pm(character(0), toFile = FALSE),
    "No input files specified"
  )

  # Test with invalid parameters
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Invalid F0 range (minF > maxF)
  expect_error(
    superassp::trk_reaper_pm(test_wav, minF = 500, maxF = 100, toFile = FALSE),
    "minF.*maxF"
  )
})

test_that("trk_reaper_pm produces consistent results", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Run twice with same parameters
  result1 <- superassp::trk_reaper_pm(
    test_wav,
    minF = 60,
    maxF = 400,
    windowShift = 10.0,
    voicing_threshold = 0.9,
    toFile = FALSE,
    verbose = FALSE
  )

  result2 <- superassp::trk_reaper_pm(
    test_wav,
    minF = 60,
    maxF = 400,
    windowShift = 10.0,
    voicing_threshold = 0.9,
    toFile = FALSE,
    verbose = FALSE
  )

  # Results should be identical
  expect_equal(nrow(result1$pm), nrow(result2$pm))
  expect_equal(result1$pm, result2$pm)
  expect_equal(attr(result1, "n_epochs"), attr(result2, "n_epochs"))
  expect_equal(attr(result1, "epoch_times"), attr(result2, "epoch_times"))
  expect_equal(attr(result1, "polarity"), attr(result2, "polarity"))
})

test_that("trk_reaper_pm handles non-WAV files via av package", {
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

  # Test that trk_reaper_pm can process MP3
  result <- superassp::trk_reaper_pm(
    temp_mp3,
    toFile = FALSE,
    verbose = FALSE
  )

  expect_s3_class(result, "AsspDataObj")
  expect_true("pm" %in% names(result))
  expect_type(result$pm, "integer")
})

test_that("trk_reaper_pm epoch times are within signal duration", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get audio duration
  info <- av::av_media_info(test_wav)
  duration <- as.numeric(info$duration)

  # Get pitch marks
  result <- superassp::trk_reaper_pm(
    test_wav,
    toFile = FALSE,
    verbose = FALSE
  )

  # Extract epoch times
  epoch_times <- attr(result, "epoch_times")

  # All epochs should be within signal duration (with small tolerance)
  if (length(epoch_times) > 0) {
    expect_true(all(epoch_times >= 0))
    expect_true(all(epoch_times <= duration + 0.1))  # Allow 100ms tolerance
  }
})

test_that("trk_reaper_pm verbose output works", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Capture output with verbose = TRUE
  output <- capture.output({
    result <- superassp::trk_reaper_pm(
      test_wav,
      toFile = FALSE,
      verbose = TRUE
    )
  })

  # Should have some output
  expect_true(length(output) > 0)

  # Output should mention processing
  expect_true(any(grepl("Processing", output, ignore.case = TRUE)))
})

test_that("trk_reaper_pm handles files with no voiced regions", {
  skip_if_not_installed("superassp")

  # This test would need a silent or whispered audio file
  # Skip if not available
  skip("Test requires audio file with no voiced regions")

  # If we had such a file:
  # silent_wav <- system.file("samples", "silent.wav", package = "superassp")
  # result <- superassp::trk_reaper_pm(silent_wav, toFile = FALSE)
  # expect_equal(attr(result, "n_epochs"), 0)
  # expect_true(all(result$pm == 0L))
})

# =============================================================================
# Summary Test
# =============================================================================

test_that("trk_reaper_pm comprehensive functionality check", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Run comprehensive test
  result <- superassp::trk_reaper_pm(
    test_wav,
    minF = 60,
    maxF = 400,
    windowShift = 10.0,
    voicing_threshold = 0.9,
    toFile = FALSE,
    verbose = FALSE
  )

  # Comprehensive checks
  checks <- list(
    is_asspobj = inherits(result, "AsspDataObj"),
    has_pm_track = "pm" %in% names(result),
    pm_is_matrix = is.matrix(result$pm),
    pm_is_integer = is.integer(result$pm),
    pm_is_binary = all(result$pm %in% c(0L, 1L)),
    has_epochs = !is.null(attr(result, "epoch_times")),
    has_n_epochs = !is.null(attr(result, "n_epochs")),
    has_polarity = !is.null(attr(result, "polarity")),
    frame_rate_correct = attr(result, "sampleRate") == 100,
    has_some_marks = sum(result$pm == 1L) > 0
  )

  # All checks should pass
  expect_true(all(unlist(checks)),
              info = paste("Failed checks:",
                          paste(names(checks)[!unlist(checks)], collapse = ", ")))
})

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

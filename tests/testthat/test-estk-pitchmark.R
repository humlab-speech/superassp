# Tests for ESTK Pitchmark (Glottal Closure Instant Detection)

test_that("estk_pitchmark_cpp works with default parameters", {
  skip_if_not_installed("superassp")

  # Get test file
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load audio
  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Run pitchmark detection
  result <- superassp::estk_pitchmark_cpp(audio_obj)

  # Check result structure
  expect_type(result, "list")
  expect_named(result, c("pitchmarks", "n_pitchmarks", "sample_rate", "duration"))

  # Check data types
  expect_type(result$pitchmarks, "double")
  expect_type(result$n_pitchmarks, "double")  # Can be double or integer
  expect_type(result$sample_rate, "double")
  expect_type(result$duration, "double")

  # Check reasonable values
  expect_true(all(result$pitchmarks >= 0))  # Times should be non-negative
  expect_true(all(result$pitchmarks <= result$duration))  # Within signal duration
  expect_true(result$sample_rate > 0)
  expect_true(result$n_pitchmarks >= 0)
  expect_equal(length(result$pitchmarks), result$n_pitchmarks)
})

test_that("estk_pitchmark_cpp filter parameters work correctly", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Test with different low-pass filter settings
  result_lp1 <- superassp::estk_pitchmark_cpp(audio_obj, lx_low_frequency = 400, lx_low_order = 19)
  result_lp2 <- superassp::estk_pitchmark_cpp(audio_obj, lx_low_frequency = 600, lx_low_order = 19)

  # Both should produce results
  expect_type(result_lp1$pitchmarks, "double")
  expect_type(result_lp2$pitchmarks, "double")

  # Test with different high-pass filter settings
  result_hp1 <- superassp::estk_pitchmark_cpp(audio_obj, lx_high_frequency = 40, lx_high_order = 19)
  result_hp2 <- superassp::estk_pitchmark_cpp(audio_obj, lx_high_frequency = 60, lx_high_order = 19)

  # Both should produce results
  expect_type(result_hp1$pitchmarks, "double")
  expect_type(result_hp2$pitchmarks, "double")

  # Test with different delta filter settings
  result_df1 <- superassp::estk_pitchmark_cpp(audio_obj, df_low_frequency = 1000, df_low_order = 19)
  result_df2 <- superassp::estk_pitchmark_cpp(audio_obj, df_low_frequency = 1500, df_low_order = 19)

  # Both should produce results
  expect_type(result_df1$pitchmarks, "double")
  expect_type(result_df2$pitchmarks, "double")
})

test_that("estk_pitchmark_cpp median_order parameter works", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Test with no smoothing
  result_no_smooth <- superassp::estk_pitchmark_cpp(audio_obj, median_order = 0)

  # Test with smoothing
  result_smooth <- superassp::estk_pitchmark_cpp(audio_obj, median_order = 19)

  # Both should produce results
  expect_type(result_no_smooth$pitchmarks, "double")
  expect_type(result_smooth$pitchmarks, "double")

  # Smoothing should typically reduce jitter (fewer pitchmarks)
  # but we're just checking both work
  expect_true(result_no_smooth$n_pitchmarks >= 0)
  expect_true(result_smooth$n_pitchmarks >= 0)
})

test_that("estk_pitchmark_cpp fill parameter works correctly", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Without filling
  result_no_fill <- superassp::estk_pitchmark_cpp(audio_obj, fill = FALSE)

  # With filling
  result_fill <- superassp::estk_pitchmark_cpp(audio_obj,
                                     fill = TRUE,
                                     min_period = 0.003,  # 3 ms (333 Hz)
                                     max_period = 0.02,   # 20 ms (50 Hz)
                                     def_period = 0.01)   # 10 ms (100 Hz)

  # Both should produce results
  expect_type(result_no_fill$pitchmarks, "double")
  expect_type(result_fill$pitchmarks, "double")

  # Filling should typically produce more regular pitchmarks
  # Check that filled pitchmarks respect period constraints
  if (result_fill$n_pitchmarks > 1) {
    periods <- diff(result_fill$pitchmarks)
    # Most periods should be within min/max range (allowing some edge cases)
    expect_true(mean(periods >= 0.003 & periods <= 0.02) > 0.8)
  }
})

test_that("estk_pitchmark_cpp period constraints work correctly", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Test with narrow period range (female voice)
  result_female <- superassp::estk_pitchmark_cpp(audio_obj,
                                       fill = TRUE,
                                       min_period = 0.004,  # 4 ms (250 Hz)
                                       max_period = 0.0067, # 6.7 ms (150 Hz)
                                       def_period = 0.005)  # 5 ms (200 Hz)

  # Test with wide period range (male voice)
  result_male <- superassp::estk_pitchmark_cpp(audio_obj,
                                     fill = TRUE,
                                     min_period = 0.0025,  # 2.5 ms (400 Hz)
                                     max_period = 0.0133,  # 13.3 ms (75 Hz)
                                     def_period = 0.008)   # 8 ms (125 Hz)

  # Both should work
  expect_type(result_female$pitchmarks, "double")
  expect_type(result_male$pitchmarks, "double")
})

test_that("estk_pitchmark_cpp invert parameter works", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Without inversion
  result_normal <- superassp::estk_pitchmark_cpp(audio_obj, invert = FALSE)

  # With inversion
  result_inverted <- superassp::estk_pitchmark_cpp(audio_obj, invert = TRUE)

  # Both should produce results
  expect_type(result_normal$pitchmarks, "double")
  expect_type(result_inverted$pitchmarks, "double")

  # Results will differ based on signal polarity
  # Just verify both work
  expect_true(result_normal$n_pitchmarks >= 0)
  expect_true(result_inverted$n_pitchmarks >= 0)
})

test_that("estk_pitchmark_cpp to_f0 conversion works", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Without F0 conversion
  result_no_f0 <- superassp::estk_pitchmark_cpp(audio_obj, to_f0 = FALSE)

  # With F0 conversion
  result_f0 <- superassp::estk_pitchmark_cpp(audio_obj, to_f0 = TRUE)

  # Check that to_f0 = FALSE doesn't include f0
  expect_false("f0" %in% names(result_no_f0))

  # Check that to_f0 = TRUE includes f0
  expect_true("f0" %in% names(result_f0))
  expect_type(result_f0$f0, "double")

  # F0 should be derived from pitchmark intervals
  if (result_f0$n_pitchmarks > 1) {
    # F0 values should be in reasonable range for speech
    f0_matrix <- result_f0$f0
    f0_values <- as.vector(f0_matrix[f0_matrix > 0])
    if (length(f0_values) > 0) {
      expect_true(all(f0_values >= 30))    # At least 30 Hz
      expect_true(all(f0_values <= 2000))  # Allow up to 2000 Hz (may have harmonics)
    }
  }
})

test_that("estk_pitchmark_cpp handles short audio files", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load and truncate audio
  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Take only first 0.1 seconds
  n_samples <- round(0.1 * attr(audio_obj, "sampleRate"))
  audio_obj$audio <- audio_obj$audio[1:n_samples, , drop = FALSE]

  # Should still work with short audio
  result <- superassp::estk_pitchmark_cpp(audio_obj)

  expect_type(result$pitchmarks, "double")
  expect_true(result$n_pitchmarks >= 0)  # May have very few pitchmarks
  expect_equal(result$duration, 0.1, tolerance = 0.01)
})

test_that("estk_pitchmark_cpp error handling", {
  skip_if_not_installed("superassp")

  # Test with invalid input
  expect_error(
    superassp::estk_pitchmark_cpp("not an audio object"),
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
    superassp::estk_pitchmark_cpp(audio_obj_bad),
    "must contain 'audio' track"
  )
})

test_that("estk_pitchmark_cpp verbose output works", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Capture output
  output <- capture.output({
    result <- superassp::estk_pitchmark_cpp(audio_obj, verbose = TRUE)
  })

  # Should have some output
  expect_true(length(output) > 0)
  expect_true(any(grepl("Processing", output)))
  expect_true(any(grepl("Applying", output)))
  expect_true(any(grepl("Found", output)))
})

test_that("estk_pitchmark_cpp produces consistent results", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Run twice with same parameters
  result1 <- superassp::estk_pitchmark_cpp(audio_obj, lx_low_frequency = 400, lx_high_frequency = 40)
  result2 <- superassp::estk_pitchmark_cpp(audio_obj, lx_low_frequency = 400, lx_high_frequency = 40)

  # Results should be identical
  expect_equal(result1$n_pitchmarks, result2$n_pitchmarks)
  expect_equal(result1$pitchmarks, result2$pitchmarks)
  expect_equal(result1$sample_rate, result2$sample_rate)
  expect_equal(result1$duration, result2$duration)
})

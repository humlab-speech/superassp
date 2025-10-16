# Tests for ESTK PDA (Pitch Detection Algorithm)

test_that("estk_pda_cpp works with default parameters", {
  skip_if_not_installed("superassp")

  # Get test file
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load audio
  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Run PDA
  result <- superassp::estk_pda_cpp(audio_obj)

  # Check result structure
  expect_type(result, "list")
  expect_named(result, c("f0", "times", "is_voiced", "sample_rate", "windowShift", "n_frames"))

  # Check data types
  expect_type(result$f0, "double")
  expect_type(result$times, "double")
  expect_type(result$is_voiced, "logical")
  expect_type(result$sample_rate, "double")
  expect_type(result$windowShift, "double")
  expect_type(result$n_frames, "integer")

  # Check data consistency
  expect_equal(length(result$f0), result$n_frames)
  expect_equal(length(result$times), result$n_frames)
  expect_equal(length(result$is_voiced), result$n_frames)

  # Check reasonable values
  expect_true(all(result$f0 >= 0))  # F0 should be non-negative
  expect_true(all(result$times >= 0))  # Times should be non-negative
  expect_true(result$sample_rate > 0)
  expect_true(result$n_frames > 0)
})

test_that("estk_pda_cpp respects minF and maxF parameters", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Test with male voice range
  result_male <- superassp::estk_pda_cpp(audio_obj, minF = 75, maxF = 200)

  # Check that voiced F0 values are within range (allowing some tolerance)
  voiced_f0 <- result_male$f0[result_male$is_voiced & result_male$f0 > 0]
  if (length(voiced_f0) > 0) {
    expect_true(all(voiced_f0 >= 70))  # Allow 5 Hz tolerance
    expect_true(all(voiced_f0 <= 205))
  }

  # Test with female voice range
  result_female <- superassp::estk_pda_cpp(audio_obj, minF = 150, maxF = 400)

  # F0 range should affect detection
  expect_type(result_female$f0, "double")
})

test_that("estk_pda_cpp windowShift parameter works correctly", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Test with 5 ms window shift
  result_5ms <- superassp::estk_pda_cpp(audio_obj, windowShift = 5.0)

  # Test with 10 ms window shift
  result_10ms <- superassp::estk_pda_cpp(audio_obj, windowShift = 10.0)

  # 5 ms should produce approximately 2x more frames than 10 ms
  ratio <- result_5ms$n_frames / result_10ms$n_frames
  expect_true(ratio >= 1.8 && ratio <= 2.2)  # Allow some tolerance

  # Check time spacing
  if (result_5ms$n_frames > 1) {
    time_diff_5ms <- result_5ms$times[2] - result_5ms$times[1]
    expect_equal(time_diff_5ms, 0.005, tolerance = 0.001)
  }

  if (result_10ms$n_frames > 1) {
    time_diff_10ms <- result_10ms$times[2] - result_10ms$times[1]
    expect_equal(time_diff_10ms, 0.010, tolerance = 0.001)
  }
})

test_that("estk_pda_cpp decimation parameter works", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Test with different decimation factors
  result_dec2 <- superassp::estk_pda_cpp(audio_obj, decimation = 2L)
  result_dec4 <- superassp::estk_pda_cpp(audio_obj, decimation = 4L)

  # Both should work
  expect_type(result_dec2$f0, "double")
  expect_type(result_dec4$f0, "double")

  # Dec 2 is more accurate but slower (we're just checking it works)
  expect_true(result_dec2$n_frames > 0)
  expect_true(result_dec4$n_frames > 0)
})

test_that("estk_pda_cpp peak_tracking parameter works", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Without peak tracking
  result_no_track <- superassp::estk_pda_cpp(audio_obj, peak_tracking = FALSE)

  # With peak tracking
  result_track <- superassp::estk_pda_cpp(audio_obj, peak_tracking = TRUE)

  # Both should produce results
  expect_type(result_no_track$f0, "double")
  expect_type(result_track$f0, "double")

  # Peak tracking should reduce F0 jumps in sustained vowels
  # (we're just checking it works, not enforcing specific behavior)
  expect_true(result_track$n_frames > 0)
})

test_that("estk_pda_cpp handles short audio files", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Load and truncate audio
  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Take only first 0.1 seconds
  n_samples <- round(0.1 * attr(audio_obj, "sampleRate"))
  audio_obj$audio <- audio_obj$audio[1:n_samples, , drop = FALSE]

  # Should still work with short audio
  result <- superassp::estk_pda_cpp(audio_obj)

  expect_type(result$f0, "double")
  expect_true(result$n_frames >= 0)  # May have very few frames
})

test_that("estk_pda_cpp error handling", {
  skip_if_not_installed("superassp")

  # Test with invalid input
  expect_error(
    superassp::estk_pda_cpp("not an audio object"),
    "must be an AsspDataObj"
  )

  # Test with invalid F0 range
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # maxF must be greater than minF
  expect_error(
    superassp::estk_pda_cpp(audio_obj, minF = 400, maxF = 60)
  )
})

test_that("estk_pda_cpp verbose output works", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- superassp::av_to_asspDataObj(test_wav)

  # Capture output
  output <- capture.output({
    result <- superassp::estk_pda_cpp(audio_obj, verbose = TRUE)
  })

  # Should have some output
  expect_true(length(output) > 0)
  expect_true(any(grepl("Processing", output)))
})

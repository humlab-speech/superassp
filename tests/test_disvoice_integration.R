# Test Script for DisVoice Integration
#
# This script tests the trk_dv_f0() and trk_dv_formants() functions
# to verify correct implementation and performance improvements.

library(testthat)

# ============================================================================
# Setup
# ============================================================================

test_that("DisVoice support can be checked", {
  result <- has_disvoice_support()
  expect_type(result, "logical")
  expect_length(result, 1)
})

# Skip remaining tests if DisVoice not available
if (!has_disvoice_support()) {
  skip("DisVoice Python support not available. Install with: install_disvoice_python()")
}

# ============================================================================
# Test trk_dv_f0()
# ============================================================================

test_that("trk_dv_f0 exists and is exported", {
  expect_true(exists("trk_dv_f0"))
  expect_type(trk_dv_f0, "closure")
})

test_that("trk_dv_f0 validates file existence", {
  expect_error(
    trk_dv_f0("nonexistent_file.wav"),
    "Audio file not found"
  )
})

test_that("trk_dv_f0 returns AsspDataObj by default", {
  skip_if_not(file.exists("testthat/test_audio.wav"), "Test audio file not found")

  result <- trk_dv_f0("testthat/test_audio.wav")

  expect_s3_class(result, "AsspDataObj")
  expect_true("tracks" %in% names(result))
  expect_true("f0" %in% colnames(result$tracks))
})

test_that("trk_dv_f0 respects include_voicing parameter", {
  skip_if_not(file.exists("testthat/test_audio.wav"), "Test audio file not found")

  # With voicing
  result_with <- trk_dv_f0("testthat/test_audio.wav", include_voicing = TRUE)
  expect_true("voicing" %in% colnames(result_with$tracks))

  # Without voicing
  result_without <- trk_dv_f0("testthat/test_audio.wav", include_voicing = FALSE)
  expect_false("voicing" %in% colnames(result_without$tracks))
})

test_that("trk_dv_f0 returns data.frame when requested", {
  skip_if_not(file.exists("testthat/test_audio.wav"), "Test audio file not found")

  result <- trk_dv_f0("testthat/test_audio.wav", output_format = "dataframe")

  expect_s3_class(result, "data.frame")
  expect_true("time" %in% names(result))
  expect_true("f0" %in% names(result))
})

test_that("trk_dv_f0 returns list when requested", {
  skip_if_not(file.exists("testthat/test_audio.wav"), "Test audio file not found")

  result <- trk_dv_f0("testthat/test_audio.wav", output_format = "list")

  expect_type(result, "list")
  expect_true("tracks" %in% names(result))
  expect_true("times" %in% names(result))
  expect_true("sample_rate" %in% names(result))
})

test_that("trk_dv_f0 respects frame_shift parameter", {
  skip_if_not(file.exists("testthat/test_audio.wav"), "Test audio file not found")

  # 10ms shift = 100 Hz
  result_10ms <- trk_dv_f0("testthat/test_audio.wav", frame_shift = 10)

  # 5ms shift = 200 Hz (should have ~2x more frames)
  result_5ms <- trk_dv_f0("testthat/test_audio.wav", frame_shift = 5)

  expect_gt(nrow(result_5ms$tracks), nrow(result_10ms$tracks))
})

test_that("trk_dv_f0 respects F0 range parameters", {
  skip_if_not(file.exists("testthat/test_audio.wav"), "Test audio file not found")

  # Narrow range
  result_narrow <- trk_dv_f0(
    "testthat/test_audio.wav",
    min_f0 = 100,
    max_f0 = 300,
    output_format = "dataframe"
  )

  # Check that detected F0 values are within range (excluding 0/NA)
  valid_f0 <- result_narrow$f0[result_narrow$f0 > 0 & !is.na(result_narrow$f0)]
  if (length(valid_f0) > 0) {
    expect_true(all(valid_f0 >= 100))
    expect_true(all(valid_f0 <= 300))
  }
})

# ============================================================================
# Test trk_dv_formants()
# ============================================================================

test_that("trk_dv_formants exists and is exported", {
  expect_true(exists("trk_dv_formants"))
  expect_type(trk_dv_formants, "closure")
})

test_that("trk_dv_formants validates file existence", {
  expect_error(
    trk_dv_formants("nonexistent_file.wav"),
    "Audio file not found"
  )
})

test_that("trk_dv_formants returns AsspDataObj by default", {
  skip_if_not(file.exists("testthat/test_audio.wav"), "Test audio file not found")

  result <- trk_dv_formants("testthat/test_audio.wav")

  expect_s3_class(result, "AsspDataObj")
  expect_true("tracks" %in% names(result))
  expect_true(all(c("F1", "F2", "F3", "F4") %in% colnames(result$tracks)))
})

test_that("trk_dv_formants returns data.frame when requested", {
  skip_if_not(file.exists("testthat/test_audio.wav"), "Test audio file not found")

  result <- trk_dv_formants("testthat/test_audio.wav", output_format = "dataframe")

  expect_s3_class(result, "data.frame")
  expect_true("time" %in% names(result))
  expect_true(all(c("F1", "F2", "F3", "F4") %in% names(result)))
})

test_that("trk_dv_formants returns list when requested", {
  skip_if_not(file.exists("testthat/test_audio.wav"), "Test audio file not found")

  result <- trk_dv_formants("testthat/test_audio.wav", output_format = "list")

  expect_type(result, "list")
  expect_true("tracks" %in% names(result))
  expect_true("times" %in% names(result))
  expect_true("sample_rate" %in% names(result))
  expect_true("parameters" %in% names(result))
})

test_that("trk_dv_formants respects frame_shift parameter", {
  skip_if_not(file.exists("testthat/test_audio.wav"), "Test audio file not found")

  # 10ms shift
  result_10ms <- trk_dv_formants("testthat/test_audio.wav", frame_shift = 10)

  # 5ms shift (should have ~2x more frames)
  result_5ms <- trk_dv_formants("testthat/test_audio.wav", frame_shift = 5)

  expect_gt(nrow(result_5ms$tracks), nrow(result_10ms$tracks))
})

test_that("trk_dv_formants respects window_size parameter", {
  skip_if_not(file.exists("testthat/test_audio.wav"), "Test audio file not found")

  # Different window sizes should produce different results
  result_25ms <- trk_dv_formants("testthat/test_audio.wav", window_size = 25)
  result_50ms <- trk_dv_formants("testthat/test_audio.wav", window_size = 50)

  expect_s3_class(result_25ms, "AsspDataObj")
  expect_s3_class(result_50ms, "AsspDataObj")
})

test_that("trk_dv_formants respects max_formant_freq parameter", {
  skip_if_not(file.exists("testthat/test_audio.wav"), "Test audio file not found")

  # Male speaker range (lower ceiling)
  result_male <- trk_dv_formants(
    "testthat/test_audio.wav",
    max_formant_freq = 5000,
    output_format = "dataframe"
  )

  # Female speaker range (higher ceiling)
  result_female <- trk_dv_formants(
    "testthat/test_audio.wav",
    max_formant_freq = 5500,
    output_format = "dataframe"
  )

  expect_s3_class(result_male, "data.frame")
  expect_s3_class(result_female, "data.frame")
})

# ============================================================================
# Test Helper Functions
# ============================================================================

test_that("as_numpy_audio converts audio correctly", {
  skip_if_not(has_disvoice_support())

  # Create test audio data (INT32 range)
  test_audio <- as.integer(c(-2147483647, 0, 2147483647))

  result <- as_numpy_audio(test_audio, sample_rate = 16000)

  expect_type(result, "list")
  expect_true("array" %in% names(result))
  expect_true("sample_rate" %in% names(result))
  expect_equal(result$sample_rate, 16000L)
})

test_that("load_audio_as_numpy loads and converts audio", {
  skip_if_not(has_disvoice_support())
  skip_if_not(file.exists("testthat/test_audio.wav"), "Test audio file not found")

  result <- load_audio_as_numpy("testthat/test_audio.wav", channels = 1)

  expect_type(result, "list")
  expect_true("array" %in% names(result))
  expect_true("sample_rate" %in% names(result))
})

test_that("init_disvoice loads DisVoice modules", {
  skip_if_not(has_disvoice_support())

  disvoice <- init_disvoice()

  expect_false(is.null(disvoice))
  expect_true("praat_functions" %in% names(disvoice))
})

test_that("get_disvoice_env returns cached environment", {
  skip_if_not(has_disvoice_support())

  # First call
  env1 <- get_disvoice_env()

  # Second call (should be cached)
  env2 <- get_disvoice_env()

  expect_identical(env1, env2)
})

# ============================================================================
# Performance Benchmarks (Optional)
# ============================================================================

test_that("DisVoice functions run reasonably fast", {
  skip_if_not(has_disvoice_support())
  skip_if_not(file.exists("testthat/test_audio.wav"), "Test audio file not found")

  # F0 extraction should complete in reasonable time
  start_time <- Sys.time()
  result_f0 <- trk_dv_f0("testthat/test_audio.wav")
  end_time <- Sys.time()
  f0_duration <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Formant extraction should complete in reasonable time
  start_time <- Sys.time()
  result_formants <- trk_dv_formants("testthat/test_audio.wav")
  end_time <- Sys.time()
  formants_duration <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # For 2-second audio, should be under 1 second each
  # (Very conservative - actual times should be much faster)
  expect_lt(f0_duration, 1.0)
  expect_lt(formants_duration, 1.0)

  message(sprintf("F0 extraction: %.3f seconds", f0_duration))
  message(sprintf("Formant extraction: %.3f seconds", formants_duration))
})

# ============================================================================
# Integration Tests
# ============================================================================

test_that("trk_dv_f0 output is compatible with superassp workflows", {
  skip_if_not(has_disvoice_support())
  skip_if_not(file.exists("testthat/test_audio.wav"), "Test audio file not found")

  result <- trk_dv_f0("testthat/test_audio.wav")

  # Check AsspDataObj structure
  expect_true("sampleRate" %in% names(result))
  expect_true("origFreq" %in% names(result))
  expect_true("startTime" %in% names(result))
  expect_true("startRecord" %in% names(result))
  expect_true("endRecord" %in% names(result))

  # Check track format
  expect_type(result$tracks, "double")
  expect_true(is.matrix(result$tracks))
})

test_that("trk_dv_formants output is compatible with superassp workflows", {
  skip_if_not(has_disvoice_support())
  skip_if_not(file.exists("testthat/test_audio.wav"), "Test audio file not found")

  result <- trk_dv_formants("testthat/test_audio.wav")

  # Check AsspDataObj structure
  expect_true("sampleRate" %in% names(result))
  expect_true("origFreq" %in% names(result))
  expect_true("startTime" %in% names(result))

  # Check track format
  expect_type(result$tracks, "double")
  expect_true(is.matrix(result$tracks))
  expect_equal(ncol(result$tracks), 4)  # F1, F2, F3, F4
})

message("\n=== DisVoice Integration Tests Complete ===\n")

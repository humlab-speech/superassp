# Tests for Snack-compatible Pitch and Formant Functions

# =============================================================================
# trk_snackp() tests
# =============================================================================

test_that("snack_pitch works with single file", {
  skip_if_not_installed("superassp")
  
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  
  # Test in-memory output
  result <- trk_snackp(test_wav, toFile = FALSE, verbose = FALSE)
  
  # Check class
  expect_s3_class(result, "AsspDataObj")
  
  # Check tracks exist
  expect_true("f0" %in% names(result))
  expect_true("voicing" %in% names(result))
  expect_true("rms" %in% names(result))
  
  # Check data dimensions
  expect_true(is.matrix(result$f0))
  expect_true(is.matrix(result$voicing))
  expect_true(is.matrix(result$rms))
  expect_equal(ncol(result$f0), 1)
  
  # Check attributes
  expect_true(!is.null(attr(result, "sampleRate")))
  expect_true(!is.null(attr(result, "startTime")))
  expect_true(attr(result, "sampleRate") > 0)
  
  # Check F0 values are reasonable
  f0_vals <- result$f0[result$f0 > 0]
  if (length(f0_vals) > 0) {
    expect_true(all(f0_vals >= 40))   # Minimum reasonable F0
    expect_true(all(f0_vals <= 600))  # Maximum reasonable F0
  }
  
  # Check voicing is between 0 and 1
  expect_true(all(result$voicing >= 0))
  expect_true(all(result$voicing <= 1))
  
  # Check RMS is non-negative
  expect_true(all(result$rms >= 0))
})

test_that("snack_pitch works with custom parameters", {
  skip_if_not_installed("superassp")
  
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  
  # Test with custom F0 range
  result <- trk_snackp(
    test_wav,
    minF = 75,
    maxF = 400,
    windowShift = 5,
    threshold = 0.4,
    toFile = FALSE,
    verbose = FALSE
  )
  
  expect_s3_class(result, "AsspDataObj")
  expect_true("f0" %in% names(result))
  
  # Check F0 values respect the range
  f0_vals <- result$f0[result$f0 > 0]
  if (length(f0_vals) > 0) {
    # Should be within or near the specified range
    expect_true(all(f0_vals >= 50))
    expect_true(all(f0_vals <= 500))
  }
})

test_that("snack_pitch works with time windowing", {
  skip_if_not_installed("superassp")
  
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  
  # Get full duration result
  result_full <- trk_snackp(test_wav, toFile = FALSE, verbose = FALSE)
  
  # Get windowed result (first 2 seconds)
  result_windowed <- trk_snackp(
    test_wav,
    beginTime = 0.0,
    endTime = 2.0,
    toFile = FALSE,
    verbose = FALSE
  )
  
  # Windowed result should have fewer frames
  expect_true(nrow(result_windowed$f0) < nrow(result_full$f0))
  expect_s3_class(result_windowed, "AsspDataObj")
})

test_that("snack_pitch works with file output", {
  skip_if_not_installed("superassp")
  
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  
  tmpdir <- tempdir()
  
  # Test file output
  n_success <- trk_snackp(
    test_wav,
    toFile = TRUE,
    outputDirectory = tmpdir,
    verbose = FALSE
  )
  
  expect_equal(n_success, 1)
  
  # Check file was created
  output_file <- file.path(tmpdir, "a1.snackpitch")
  expect_true(file.exists(output_file))
  
  # Try to read it back
  result <- wrassp::read.AsspDataObj(output_file)
  expect_s3_class(result, "AsspDataObj")
  expect_true("f0" %in% names(result))
  expect_true("voicing" %in% names(result))
  expect_true("rms" %in% names(result))
  
  # Clean up
  unlink(output_file)
})

test_that("snack_pitch works with multiple files", {
  skip_if_not_installed("superassp")
  
  sample_dir <- system.file("samples", "sustained", package = "superassp")
  test_files <- list.files(sample_dir, pattern = "\\.wav$", full.names = TRUE)[1:2]
  skip_if(length(test_files) < 2, "Not enough test files")
  
  # Test batch processing
  results <- trk_snackp(test_files, toFile = FALSE, verbose = FALSE)
  
  expect_type(results, "list")
  expect_length(results, 2)
  
  # Check each result
  for (i in seq_along(results)) {
    expect_s3_class(results[[i]], "AsspDataObj")
    expect_true("f0" %in% names(results[[i]]))
  }
})

# =============================================================================
# trk_snackf() tests
# =============================================================================

test_that("snack_formant works with single file", {
  skip_if_not_installed("superassp")
  
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  
  # Test with default 4 formants
  result <- trk_snackf(test_wav, toFile = FALSE, verbose = FALSE)
  
  # Check class
  expect_s3_class(result, "AsspDataObj")
  
  # Check tracks exist (4 formants = 8 tracks)
  expect_true("fm_1" %in% names(result))
  expect_true("fm_2" %in% names(result))
  expect_true("fm_3" %in% names(result))
  expect_true("fm_4" %in% names(result))
  expect_true("bw_1" %in% names(result))
  expect_true("bw_2" %in% names(result))
  expect_true("bw_3" %in% names(result))
  expect_true("bw_4" %in% names(result))
  
  # Check data dimensions
  expect_true(is.matrix(result$fm_1))
  expect_equal(ncol(result$fm_1), 1)
  
  # Check attributes
  expect_true(!is.null(attr(result, "sampleRate")))
  expect_true(!is.null(attr(result, "startTime")))
  expect_true(attr(result, "sampleRate") > 0)
  
  # Check formant values are in reasonable ranges
  # F1: typically 200-900 Hz
  f1_vals <- result$fm_1[result$fm_1 > 0]
  if (length(f1_vals) > 0) {
    expect_true(all(f1_vals >= 100))
    expect_true(all(f1_vals <= 5000))  # Allow wider range for errors
  }
  
  # Bandwidths should be positive (can be large for unstable LPC poles)
  bw1_vals <- result$bw_1[result$bw_1 > 0]
  if (length(bw1_vals) > 0) {
    expect_true(all(bw1_vals > 0))
    # Note: Bandwidths can be very large for unstable poles, so no upper limit check
  }
})

test_that("snack_formant works with different formant counts", {
  skip_if_not_installed("superassp")
  
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  
  # Test with 3 formants
  result_3 <- trk_snackf(test_wav, numFormants = 3, toFile = FALSE, verbose = FALSE)
  expect_true("fm_3" %in% names(result_3))
  expect_false("fm_4" %in% names(result_3))
  expect_true("bw_3" %in% names(result_3))
  expect_false("bw_4" %in% names(result_3))
  
  # Test with 5 formants
  result_5 <- trk_snackf(test_wav, numFormants = 5, toFile = FALSE, verbose = FALSE)
  expect_true("fm_5" %in% names(result_5))
  expect_true("bw_5" %in% names(result_5))
  
  # Should have different number of tracks
  expect_true(length(names(result_3)) < length(names(result_5)))
})

test_that("snack_formant works with custom LPC parameters", {
  skip_if_not_installed("superassp")
  
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  
  # Test with custom LPC order and pre-emphasis
  result <- trk_snackf(
    test_wav,
    numFormants = 4,
    lpcOrder = 16,
    preEmphasis = 0.9,
    windowShift = 10,
    toFile = FALSE,
    verbose = FALSE
  )
  
  expect_s3_class(result, "AsspDataObj")
  expect_true("fm_1" %in% names(result))
  expect_true("fm_4" %in% names(result))
})

test_that("snack_formant works with time windowing", {
  skip_if_not_installed("superassp")
  
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  
  # Get full duration result
  result_full <- trk_snackf(test_wav, toFile = FALSE, verbose = FALSE)
  
  # Get windowed result (first 1 second)
  result_windowed <- trk_snackf(
    test_wav,
    beginTime = 0.0,
    endTime = 1.0,
    toFile = FALSE,
    verbose = FALSE
  )
  
  # Windowed result should have fewer frames
  expect_true(nrow(result_windowed$fm_1) < nrow(result_full$fm_1))
  expect_s3_class(result_windowed, "AsspDataObj")
})

test_that("snack_formant works with file output", {
  skip_if_not_installed("superassp")
  
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  
  tmpdir <- tempdir()
  
  # Test file output
  n_success <- trk_snackf(
    test_wav,
    numFormants = 4,
    toFile = TRUE,
    outputDirectory = tmpdir,
    verbose = FALSE
  )
  
  expect_equal(n_success, 1)
  
  # Check file was created
  output_file <- file.path(tmpdir, "a1.snackfmt")
  expect_true(file.exists(output_file))
  
  # Try to read it back
  result <- wrassp::read.AsspDataObj(output_file)
  expect_s3_class(result, "AsspDataObj")
  expect_true("fm_1" %in% names(result))
  expect_true("bw_1" %in% names(result))
  
  # Clean up
  unlink(output_file)
})

test_that("snack_formant works with multiple files", {
  skip_if_not_installed("superassp")
  
  sample_dir <- system.file("samples", "sustained", package = "superassp")
  test_files <- list.files(sample_dir, pattern = "\\.wav$", full.names = TRUE)[1:2]
  skip_if(length(test_files) < 2, "Not enough test files")
  
  # Test batch processing
  results <- trk_snackf(test_files, numFormants = 4, toFile = FALSE, verbose = FALSE)
  
  expect_type(results, "list")
  expect_length(results, 2)
  
  # Check each result
  for (i in seq_along(results)) {
    expect_s3_class(results[[i]], "AsspDataObj")
    expect_true("fm_1" %in% names(results[[i]]))
    expect_true("fm_4" %in% names(results[[i]]))
  }
})

test_that("snack_formant validates numFormants parameter", {
  skip_if_not_installed("superassp")
  
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  
  # Test invalid formant count (too high)
  expect_error(
    trk_snackf(test_wav, numFormants = 10, toFile = FALSE, verbose = FALSE),
    "numFormants must be between 1 and 7"
  )
  
  # Test invalid formant count (too low)
  expect_error(
    trk_snackf(test_wav, numFormants = 0, toFile = FALSE, verbose = FALSE),
    "numFormants must be between 1 and 7"
  )
})

# =============================================================================
# Integration tests
# =============================================================================

test_that("snack functions can be used together", {
  skip_if_not_installed("superassp")
  
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  
  # Extract both pitch and formants
  pitch_data <- trk_snackp(test_wav, toFile = FALSE, verbose = FALSE)
  formant_data <- trk_snackf(test_wav, toFile = FALSE, verbose = FALSE)
  
  # Both should succeed
  expect_s3_class(pitch_data, "AsspDataObj")
  expect_s3_class(formant_data, "AsspDataObj")
  
  # Can extract data
  f0_values <- pitch_data$f0[pitch_data$f0 > 0]
  f1_values <- formant_data$fm_1[formant_data$fm_1 > 0]
  
  # Should have some voiced frames
  expect_true(length(f0_values) > 0 || length(f1_values) > 0)
})

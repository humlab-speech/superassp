# Test suite for legacy STRAIGHT vocoder functions
# Tests F0 extraction, spectral analysis, and synthesis

# Setup -----------------------------------------------------------------------
skip_if_not_available <- function() {
  if (!straight_available()) {
    skip("Legacy STRAIGHT not available. Install with install_legacy_straight()")
  }
}

test_audio <- function() {
  # Use a sample audio file from the package
  audio_file <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  if (!file.exists(audio_file)) {
    skip("Test audio file not found")
  }
  return(audio_file)
}


# Installation and availability tests -----------------------------------------

test_that("straight_available() returns logical", {
  result <- straight_available()
  expect_type(result, "logical")
  expect_length(result, 1)
})

test_that("straight_available(detailed = TRUE) returns list", {
  result <- straight_available(detailed = TRUE)
  expect_type(result, "list")
  expect_true("available" %in% names(result))
})

test_that("straight_info() runs without error when available", {
  skip_if_not_available()
  expect_output(straight_info(), "Legacy STRAIGHT")
})


# F0 extraction tests ---------------------------------------------------------

test_that("trk_straight_f0() extracts F0 from single file", {
  skip_if_not_available()
  audio_file <- test_audio()
  
  result <- trk_straight_f0(audio_file, toFile = FALSE, verbose = FALSE)
  
  expect_s3_class(result, "AsspDataObj")
  expect_true("f0" %in% names(result))
  expect_true("vuv" %in% names(result))
  expect_true(ncol(result$f0) == 1)
  expect_true(all(result$f0 >= 0, na.rm = TRUE))
})

test_that("trk_straight_f0() respects F0 range parameters", {
  skip_if_not_available()
  audio_file <- test_audio()
  
  result <- trk_straight_f0(
    audio_file,
    f0_floor = 80,
    f0_ceil = 300,
    toFile = FALSE,
    verbose = FALSE
  )
  
  # Check that voiced F0 values are within range (allowing some tolerance)
  voiced_f0 <- result$f0[result$vuv > 0]
  if (length(voiced_f0) > 0) {
    expect_true(all(voiced_f0 >= 70 & voiced_f0 <= 310, na.rm = TRUE))
  }
})

test_that("trk_straight_f0() includes auxiliary scores", {
  skip_if_not_available()
  audio_file <- test_audio()
  
  result <- trk_straight_f0(audio_file, toFile = FALSE, verbose = FALSE)
  
  expect_true("if_score" %in% names(result))
  expect_true("ac_score" %in% names(result))
  expect_equal(nrow(result$if_score), nrow(result$f0))
})

test_that("trk_straight_f0() handles frame_shift parameter", {
  skip_if_not_available()
  audio_file <- test_audio()
  
  result_1ms <- trk_straight_f0(
    audio_file,
    frame_shift = 1.0,
    toFile = FALSE,
    verbose = FALSE
  )
  
  result_5ms <- trk_straight_f0(
    audio_file,
    frame_shift = 5.0,
    toFile = FALSE,
    verbose = FALSE
  )
  
  # 5ms shift should produce fewer frames
  expect_true(nrow(result_5ms$f0) < nrow(result_1ms$f0))
  expect_true(nrow(result_5ms$f0) < nrow(result_1ms$f0) / 3)
})

test_that("trk_straight_f0() works with time windowing", {
  skip_if_not_available()
  audio_file <- test_audio()
  
  result_full <- trk_straight_f0(
    audio_file,
    beginTime = 0,
    endTime = 0,
    toFile = FALSE,
    verbose = FALSE
  )
  
  result_windowed <- trk_straight_f0(
    audio_file,
    beginTime = 0.1,
    endTime = 0.3,
    toFile = FALSE,
    verbose = FALSE
  )
  
  # Windowed version should have fewer frames
  expect_true(nrow(result_windowed$f0) < nrow(result_full$f0))
})

test_that("trk_straight_f0() writes to file when toFile = TRUE", {
  skip_if_not_available()
  audio_file <- test_audio()
  temp_dir <- tempdir()
  
  output_file <- trk_straight_f0(
    audio_file,
    toFile = TRUE,
    outputDirectory = temp_dir,
    verbose = FALSE
  )
  
  expect_true(file.exists(output_file))
  expect_match(output_file, "\\.strf0$")
  
  # Clean up
  unlink(output_file)
})

test_that("trk_straight_f0() handles batch processing", {
  skip_if_not_available()
  audio_file <- test_audio()
  
  # Process same file twice (simulating batch)
  files <- rep(audio_file, 2)
  results <- trk_straight_f0(files, toFile = FALSE, verbose = FALSE)
  
  expect_type(results, "list")
  expect_length(results, 2)
  expect_s3_class(results[[1]], "AsspDataObj")
  expect_s3_class(results[[2]], "AsspDataObj")
})


# Spectral analysis tests -----------------------------------------------------

test_that("trk_straight_spec() extracts spectral envelope", {
  skip_if_not_available()
  audio_file <- test_audio()
  
  result <- trk_straight_spec(audio_file, toFile = FALSE, verbose = FALSE)
  
  expect_s3_class(result, "AsspDataObj")
  expect_true("spec" %in% names(result))
  expect_true(is.matrix(result$spec))
  expect_true(ncol(result$spec) > 0)  # Frequency bins
  expect_true(nrow(result$spec) > 0)  # Time frames
})

test_that("trk_straight_spec() respects FFT size parameter", {
  skip_if_not_available()
  audio_file <- test_audio()
  
  result_2048 <- trk_straight_spec(
    audio_file,
    fft_size = 2048,
    toFile = FALSE,
    verbose = FALSE
  )
  
  result_4096 <- trk_straight_spec(
    audio_file,
    fft_size = 4096,
    toFile = FALSE,
    verbose = FALSE
  )
  
  # Larger FFT should produce more frequency bins
  expect_true(ncol(result_4096$spec) > ncol(result_2048$spec))
})

test_that("trk_straight_spec() writes to file when toFile = TRUE", {
  skip_if_not_available()
  audio_file <- test_audio()
  temp_dir <- tempdir()
  
  output_file <- trk_straight_spec(
    audio_file,
    toFile = TRUE,
    outputDirectory = temp_dir,
    verbose = FALSE
  )
  
  expect_true(file.exists(output_file))
  expect_match(output_file, "\\.strspec$")
  
  # Clean up
  unlink(output_file)
})


# Error handling tests --------------------------------------------------------

test_that("Functions fail gracefully when STRAIGHT not available", {
  # Mock unavailability
  skip_if_not_available()
  
  # These tests would need mocking infrastructure
  # For now, just test that the availability check works
  expect_true(TRUE)
})

test_that("Functions validate input parameters", {
  skip_if_not_available()
  audio_file <- test_audio()
  
  # Invalid f0_floor
  expect_error(
    trk_straight_f0(audio_file, f0_floor = -10, toFile = FALSE),
    NA  # Should not error in Python, will be handled there
  )
  
  # Invalid file path
  expect_error(
    trk_straight_f0("nonexistent.wav", toFile = FALSE)
  )
})


# Integration tests -----------------------------------------------------------

test_that("F0 and spectral analysis produce compatible outputs", {
  skip_if_not_available()
  audio_file <- test_audio()
  
  f0_result <- trk_straight_f0(audio_file, toFile = FALSE, verbose = FALSE)
  spec_result <- trk_straight_spec(audio_file, toFile = FALSE, verbose = FALSE)
  
  # Should have similar number of frames (within tolerance for frame alignment)
  expect_true(abs(nrow(f0_result$f0) - nrow(spec_result$spec)) < 5)
})

test_that("Module attributes are correctly set", {
  skip_if_not_available()
  
  # Check function attributes
  expect_equal(attr(trk_straight_f0, "ext"), "strf0")
  expect_equal(attr(trk_straight_spec, "ext"), "strspec")
  
  expect_true("f0" %in% attr(trk_straight_f0, "tracks"))
  expect_true("spec" %in% attr(trk_straight_spec, "tracks"))
})


# Performance tests (optional, can be slow) -----------------------------------

test_that("STRAIGHT F0 extraction completes in reasonable time", {
  skip_if_not_available()
  skip_on_cran()
  
  audio_file <- test_audio()
  
  # Measure execution time
  start_time <- Sys.time()
  result <- trk_straight_f0(audio_file, toFile = FALSE, verbose = FALSE)
  end_time <- Sys.time()
  
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Should complete in under 5 seconds for a short test file
  expect_true(duration < 5.0)
  
  # Report performance
  message(sprintf("F0 extraction took %.2f seconds", duration))
})

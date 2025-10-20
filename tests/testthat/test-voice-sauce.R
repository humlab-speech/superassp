test_that("lst_voice_sauce works with single file (basic)", {
  skip_if_not(voice_sauce_available(), "VoiceSauce not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_voice_sauce(test_wav, verbose = FALSE)

  expect_type(result, "list")

  # Check for key measures
  expect_true("F0" %in% names(result))
  expect_true("F1" %in% names(result))
  expect_true("F2" %in% names(result))
  expect_true("CPP" %in% names(result))
  expect_true("HNR05" %in% names(result))
  expect_true("Energy" %in% names(result))
  expect_true("H1" %in% names(result))
  expect_true("H2" %in% names(result))

  # Values should be numeric vectors
  expect_type(result$F0, "double")
  expect_type(result$CPP, "double")
  expect_true(length(result$F0) > 0)
})

test_that("lst_voice_sauce respects custom F0 range", {
  skip_if_not(voice_sauce_available(), "VoiceSauce not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_voice_sauce(test_wav,
                           f0_min = 80,
                           f0_max = 300,
                           verbose = FALSE)

  expect_type(result, "list")
  expect_true("F0" %in% names(result))

  # F0 values should be within range (excluding NAs)
  valid_f0 <- result$F0[!is.na(result$F0) & result$F0 > 0]
  if (length(valid_f0) > 0) {
    expect_true(all(valid_f0 >= 80 & valid_f0 <= 300))
  }
})

test_that("lst_voice_sauce works with different F0 methods", {
  skip_if_not(voice_sauce_available(), "VoiceSauce not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Test REAPER (default)
  result_reaper <- lst_voice_sauce(test_wav,
                                  f0_method = "reaper",
                                  verbose = FALSE)

  # Test Praat
  result_praat <- lst_voice_sauce(test_wav,
                                 f0_method = "praat",
                                 verbose = FALSE)

  expect_type(result_reaper, "list")
  expect_type(result_praat, "list")

  # Both should have F0
  expect_true("F0" %in% names(result_reaper))
  expect_true("F0" %in% names(result_praat))

  # Results may differ between methods
  expect_true(length(result_reaper$F0) > 0)
  expect_true(length(result_praat$F0) > 0)
})

test_that("lst_voice_sauce supports custom frame shift", {
  skip_if_not(voice_sauce_available(), "VoiceSauce not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result_1ms <- lst_voice_sauce(test_wav,
                               frame_shift = 1.0,
                               verbose = FALSE)

  result_10ms <- lst_voice_sauce(test_wav,
                                frame_shift = 10.0,
                                verbose = FALSE)

  expect_type(result_1ms, "list")
  expect_type(result_10ms, "list")

  # More frames with smaller shift
  expect_true(length(result_1ms$F0) > length(result_10ms$F0))
})

test_that("lst_voice_sauce supports time windowing", {
  skip_if_not(voice_sauce_available(), "VoiceSauce not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result_full <- lst_voice_sauce(test_wav, verbose = FALSE)
  result_window <- lst_voice_sauce(test_wav,
                                   beginTime = 0.5,
                                   endTime = 1.5,
                                   verbose = FALSE)

  expect_type(result_full, "list")
  expect_type(result_window, "list")

  # Both should have same parameters
  expect_true(all(c("F0", "CPP", "Energy") %in% names(result_full)))
  expect_true(all(c("F0", "CPP", "Energy") %in% names(result_window)))

  # Windowed version should have fewer frames
  expect_true(length(result_window$F0) < length(result_full$F0))
})

test_that("lst_voice_sauce handles batch processing", {
  skip_if_not(voice_sauce_available(), "VoiceSauce not installed")

  test_files <- list.files(
    system.file("samples", "sustained", package = "superassp"),
    pattern = "\\.wav$",
    full.names = TRUE
  )

  skip_if(length(test_files) < 2, "Need at least 2 test files")
  test_files <- test_files[1:2]

  results <- lst_voice_sauce(test_files, verbose = FALSE)

  expect_type(results, "list")
  expect_length(results, 2)
  expect_type(results[[1]], "list")
  expect_type(results[[2]], "list")

  # Both should have all parameters
  expect_true("CPP" %in% names(results[[1]]))
  expect_true("CPP" %in% names(results[[2]]))
  expect_true("F0" %in% names(results[[1]]))
  expect_true("F0" %in% names(results[[2]]))
})

test_that("lst_voice_sauce handles missing files gracefully", {
  skip_if_not(voice_sauce_available(), "VoiceSauce not installed")

  expect_warning(
    result <- lst_voice_sauce("nonexistent.wav", verbose = FALSE),
    "File not found"
  )

  expect_null(result)
})

test_that("lst_voice_sauce validates parameters", {
  skip_if_not(voice_sauce_available(), "VoiceSauce not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Invalid frame_shift
  expect_error(
    lst_voice_sauce(test_wav, frame_shift = -1, verbose = FALSE),
    "frame_shift must be a positive number"
  )

  # Invalid F0 range
  expect_error(
    lst_voice_sauce(test_wav, f0_min = 300, f0_max = 100, verbose = FALSE),
    "f0_max must be greater than f0_min"
  )

  # Invalid F0 method
  expect_error(
    lst_voice_sauce(test_wav, f0_method = "invalid", verbose = FALSE),
    "f0_method must be one of"
  )
})

test_that("lst_voice_sauce provides consistent results", {
  skip_if_not(voice_sauce_available(), "VoiceSauce not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result1 <- lst_voice_sauce(test_wav, verbose = FALSE)
  result2 <- lst_voice_sauce(test_wav, verbose = FALSE)

  # Should have same structure
  expect_equal(names(result1), names(result2))

  # Should have same number of frames
  expect_equal(length(result1$F0), length(result2$F0))

  # Values should be identical (deterministic algorithms)
  expect_equal(result1$F0, result2$F0)
  expect_equal(result1$CPP, result2$CPP)
})

test_that("lst_voice_sauce returns expected measures", {
  skip_if_not(voice_sauce_available(), "VoiceSauce not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_voice_sauce(test_wav, verbose = FALSE)

  # Check for all expected measure categories
  expected_measures <- c(
    # F0
    "F0",
    # Formants
    "F1", "F2", "F3", "F4", "F5",
    # Bandwidths
    "B1", "B2", "B3", "B4", "B5",
    # Harmonics
    "H1", "H2", "H4",
    # Formant amplitudes
    "A1", "A2", "A3",
    # Voice quality
    "CPP", "Energy",
    # HNR
    "HNR05", "HNR15", "HNR25", "HNR35",
    # Spectral
    "2K", "5K"
  )

  for (measure in expected_measures) {
    expect_true(measure %in% names(result),
               info = sprintf("Missing measure: %s", measure))
  }

  # All measures should be numeric vectors
  for (measure in expected_measures) {
    expect_type(result[[measure]], "double",
               info = sprintf("Measure %s is not numeric", measure))
  }
})

test_that("lst_voice_sauce returns corrected measures", {
  skip_if_not(voice_sauce_available(), "VoiceSauce not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_voice_sauce(test_wav, verbose = FALSE)

  # Check for Iseli-Alwan corrected measures
  corrected_measures <- c(
    "H1c", "H2c", "H4c",
    "A1c", "A2c", "A3c",
    "H1H2c", "H1A1c", "H1A2c", "H1A3c"
  )

  for (measure in corrected_measures) {
    expect_true(measure %in% names(result),
               info = sprintf("Missing corrected measure: %s", measure))
  }
})

test_that("lst_voice_sauce includes timing information", {
  skip_if_not(voice_sauce_available(), "VoiceSauce not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_voice_sauce(test_wav, verbose = FALSE)

  # Should have times and fs
  expect_true("times" %in% names(result))
  expect_true("fs" %in% names(result))

  # Times should match number of frames
  expect_equal(length(result$times), length(result$F0))

  # Times should be increasing
  expect_true(all(diff(result$times) > 0))
})

test_that("lst_voice_sauce custom window size", {
  skip_if_not(voice_sauce_available(), "VoiceSauce not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result_25ms <- lst_voice_sauce(test_wav,
                                window_size = 25.0,
                                verbose = FALSE)

  result_50ms <- lst_voice_sauce(test_wav,
                                window_size = 50.0,
                                verbose = FALSE)

  expect_type(result_25ms, "list")
  expect_type(result_50ms, "list")

  # Both should have F0
  expect_true("F0" %in% names(result_25ms))
  expect_true("F0" %in% names(result_50ms))
})

test_that("lst_voice_sauce batch with mixed parameters", {
  skip_if_not(voice_sauce_available(), "VoiceSauce not installed")

  test_files <- list.files(
    system.file("samples", "sustained", package = "superassp"),
    pattern = "\\.wav$",
    full.names = TRUE
  )

  skip_if(length(test_files) < 2, "Need at least 2 test files")
  test_files <- test_files[1:2]

  # All files use same parameters
  results <- lst_voice_sauce(test_files,
                            f0_min = 75,
                            f0_max = 400,
                            verbose = FALSE)

  expect_length(results, 2)

  # Both should have valid results
  expect_true(!is.null(results[[1]]))
  expect_true(!is.null(results[[2]]))
})

test_that("lst_voice_sauce summary statistics calculation", {
  skip_if_not(voice_sauce_available(), "VoiceSauce not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_voice_sauce(test_wav, verbose = FALSE)

  # Calculate summary stats
  f0_valid <- result$F0[!is.na(result$F0) & result$F0 > 0]
  cpp_valid <- result$CPP[!is.na(result$CPP)]

  # Should be able to compute summaries
  if (length(f0_valid) > 0) {
    f0_mean <- mean(f0_valid)
    f0_median <- median(f0_valid)
    expect_true(is.numeric(f0_mean))
    expect_true(is.numeric(f0_median))
  }

  if (length(cpp_valid) > 0) {
    cpp_mean <- mean(cpp_valid)
    expect_true(is.numeric(cpp_mean))
  }
})

test_that("lst_voice_sauce harmonic differences are computed", {
  skip_if_not(voice_sauce_available(), "VoiceSauce not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_voice_sauce(test_wav, verbose = FALSE)

  # Check harmonic differences
  expect_true("H1H2" %in% names(result))
  expect_true("H2H4" %in% names(result))

  # Check harmonic-formant differences
  expect_true("H1A1" %in% names(result))
  expect_true("H1A2" %in% names(result))
  expect_true("H1A3" %in% names(result))

  # All should be numeric
  expect_type(result$H1H2, "double")
  expect_type(result$H1A1, "double")
})

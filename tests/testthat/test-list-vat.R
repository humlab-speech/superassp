# Tests for lst_vat - Voice Analysis Toolbox

test_that("voice_analysis_available() works", {
  skip_if_not_installed("superassp")

  # This should return TRUE or FALSE, not error
  result <- voice_analysis_available()
  expect_type(result, "logical")
  expect_length(result, 1)
})

test_that("lst_vat() works with single file", {
  skip_if_not_installed("superassp")
  skip_if(!voice_analysis_available(),
          "voice_analysis module not available. Install with: install_voice_analysis()")

  # Get test file
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Run voice analysis
  result <- lst_vat(test_wav, verbose = FALSE)

  # Check result structure
  expect_type(result, "list")
  expect_named(result, c("measures", "fs", "success", "error", "file"), ignore.order = TRUE)

  # Check success
  expect_true(result$success)
  expect_null(result$error)

  # Check measures
  expect_type(result$measures, "list")
  expect_true(length(result$measures) > 100)  # Should have ~132 measures

  # Check sample rate
  expect_type(result$fs, "integer")
  expect_true(result$fs > 0)

  # Check file path
  expect_equal(result$file, test_wav)

  # Check that measures have reasonable values (not all NA/NULL)
  non_null_measures <- sum(!sapply(result$measures, is.null))
  expect_true(non_null_measures > 50)  # At least half should be computed
})

test_that("lst_vat() computes expected feature categories", {
  skip_if_not_installed("superassp")
  skip_if(!voice_analysis_available(),
          "voice_analysis module not available")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_vat(test_wav, verbose = FALSE)

  measure_names <- names(result$measures)

  # Check for key measures from different categories

  # Jitter measures
  expect_true(any(grepl("jitter", measure_names, ignore.case = TRUE)))

  # Shimmer measures
  expect_true(any(grepl("shimmer", measure_names, ignore.case = TRUE)))

  # HNR/NHR
  expect_true(any(grepl("HNR|NHR", measure_names)))

  # Nonlinear dynamics
  expect_true("DFA" %in% measure_names)
  expect_true("RPDE" %in% measure_names)
  expect_true("PPE" %in% measure_names)

  # GNE
  expect_true(any(grepl("GNE", measure_names)))

  # MFCCs
  expect_true(any(grepl("MFCC", measure_names, ignore.case = TRUE)))

  # Wavelet
  expect_true(any(grepl("wavelet", measure_names, ignore.case = TRUE)))
})

test_that("lst_vat() works with time windowing", {
  skip_if_not_installed("superassp")
  skip_if(!voice_analysis_available(),
          "voice_analysis module not available")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get duration
  info <- av::av_media_info(test_wav)
  duration <- info$duration

  # Analyze with time window
  result <- lst_vat(test_wav,
                   beginTime = 0.1,
                   endTime = min(1.0, duration - 0.1),
                   verbose = FALSE)

  expect_true(result$success)
  expect_type(result$measures, "list")
})

test_that("lst_vat() works with custom F0 range", {
  skip_if_not_installed("superassp")
  skip_if(!voice_analysis_available(),
          "voice_analysis module not available")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Test with female voice range
  result_female <- lst_vat(test_wav,
                          f0_min = 100,
                          f0_max = 500,
                          verbose = FALSE)
  expect_true(result_female$success)

  # Test with male voice range
  result_male <- lst_vat(test_wav,
                        f0_min = 50,
                        f0_max = 250,
                        verbose = FALSE)
  expect_true(result_male$success)

  # Results should differ slightly due to different F0 ranges
  # But both should succeed
  expect_type(result_female$measures, "list")
  expect_type(result_male$measures, "list")
})

test_that("lst_vat() works with different F0 algorithms", {
  skip_if_not_installed("superassp")
  skip_if(!voice_analysis_available(),
          "voice_analysis module not available")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Test SWIPE (default)
  result_swipe <- lst_vat(test_wav,
                         f0_algorithm = "SWIPE",
                         verbose = FALSE)
  expect_true(result_swipe$success)

  # Test PRAAT
  result_praat <- lst_vat(test_wav,
                         f0_algorithm = "PRAAT",
                         verbose = FALSE)
  expect_true(result_praat$success)

  # Both should produce results
  expect_type(result_swipe$measures, "list")
  expect_type(result_praat$measures, "list")
})

test_that("lst_vat() works with return_f0 option", {
  skip_if_not_installed("superassp")
  skip_if(!voice_analysis_available(),
          "voice_analysis module not available")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # With F0 contour
  result_with_f0 <- lst_vat(test_wav,
                           return_f0 = TRUE,
                           verbose = FALSE)

  expect_true("f0" %in% names(result_with_f0))
  expect_type(result_with_f0$f0, "list")
  expect_true(length(result_with_f0$f0) > 0)

  # Without F0 contour (default)
  result_without_f0 <- lst_vat(test_wav,
                              return_f0 = FALSE,
                              verbose = FALSE)

  expect_false("f0" %in% names(result_without_f0))
})

test_that("lst_vat() works with multiple files", {
  skip_if_not_installed("superassp")
  skip_if(!voice_analysis_available(),
          "voice_analysis module not available")

  # Get multiple test files
  test_files <- system.file("samples", "sustained",
                           c("a1.wav", "a32b.wav"),
                           package = "superassp")
  test_files <- test_files[file.exists(test_files)]

  skip_if(length(test_files) < 2, "Not enough test files found")

  # Analyze multiple files
  results <- lst_vat(test_files, verbose = FALSE)

  # Should return list of results
  expect_type(results, "list")
  expect_equal(length(results), length(test_files))

  # Each result should be valid
  for (i in seq_along(results)) {
    expect_type(results[[i]], "list")
    expect_true(results[[i]]$success)
    expect_type(results[[i]]$measures, "list")
    expect_equal(results[[i]]$file, test_files[i])
  }
})

test_that("lst_vat() works with parallel processing", {
  skip_if_not_installed("superassp")
  skip_if(!voice_analysis_available(),
          "voice_analysis module not available")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get system info
  sys_info <- voice_analysis_info()
  n_cores <- min(2, sys_info$recommended_workers)

  # Test with parallel processing
  result <- lst_vat(test_wav,
                   n_cores = n_cores,
                   verbose = FALSE)

  expect_true(result$success)
  expect_type(result$measures, "list")
})

test_that("lst_vat() handles thesis mode", {
  skip_if_not_installed("superassp")
  skip_if(!voice_analysis_available(),
          "voice_analysis module not available")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Test MATLAB mode (default)
  result_matlab <- lst_vat(test_wav,
                          use_thesis_mode = FALSE,
                          verbose = FALSE)
  expect_true(result_matlab$success)

  # Test thesis mode
  result_thesis <- lst_vat(test_wav,
                          use_thesis_mode = TRUE,
                          verbose = FALSE)
  expect_true(result_thesis$success)

  # Both should succeed
  expect_type(result_matlab$measures, "list")
  expect_type(result_thesis$measures, "list")
})

test_that("lst_vat() error handling works", {
  skip_if_not_installed("superassp")
  skip_if(!voice_analysis_available(),
          "voice_analysis module not available")

  # Test with non-existent file
  expect_error(
    lst_vat("/nonexistent/file.wav"),
    "not found"
  )

  # Test with empty file list
  expect_error(
    lst_vat(character(0)),
    "No input files"
  )

  # Test with invalid F0 algorithm
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  expect_error(
    lst_vat(test_wav, f0_algorithm = "INVALID"),
    "should be one of"
  )
})

test_that("lst_vat() handles non-WAV files via av package", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("av")
  skip_if(!voice_analysis_available(),
          "voice_analysis module not available")

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

  # Test that lst_vat can process MP3
  result <- lst_vat(temp_mp3, verbose = FALSE)

  expect_true(result$success)
  expect_type(result$measures, "list")
})

test_that("lst_vat() produces consistent results", {
  skip_if_not_installed("superassp")
  skip_if(!voice_analysis_available(),
          "voice_analysis module not available")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Run twice with same parameters
  result1 <- lst_vat(test_wav,
                    f0_min = 60,
                    f0_max = 400,
                    verbose = FALSE)

  result2 <- lst_vat(test_wav,
                    f0_min = 60,
                    f0_max = 400,
                    verbose = FALSE)

  # Both should succeed
  expect_true(result1$success)
  expect_true(result2$success)

  # Should have same number of measures
  expect_equal(length(result1$measures), length(result2$measures))

  # Key measures should be very similar (allowing for numerical precision)
  if ("DFA" %in% names(result1$measures) && "DFA" %in% names(result2$measures)) {
    expect_equal(result1$measures$DFA, result2$measures$DFA, tolerance = 1e-6)
  }

  if ("RPDE" %in% names(result1$measures) && "RPDE" %in% names(result2$measures)) {
    expect_equal(result1$measures$RPDE, result2$measures$RPDE, tolerance = 1e-6)
  }
})

test_that("voice_analysis_info() works when module is available", {
  skip_if_not_installed("superassp")
  skip_if(!voice_analysis_available(),
          "voice_analysis module not available")

  info <- voice_analysis_info()

  # Check structure
  expect_type(info, "list")

  # Check expected fields
  expect_true("cpu_count" %in% names(info))
  expect_true("cpu_count_physical" %in% names(info))
  expect_true("platform" %in% names(info))
  expect_true("cython_available" %in% names(info))
  expect_true("numba_available" %in% names(info))
  expect_true("recommended_workers" %in% names(info))

  # Check types
  expect_type(info$cpu_count, "integer")
  expect_type(info$recommended_workers, "integer")
  expect_type(info$cython_available, "logical")
  expect_type(info$numba_available, "logical")
})

test_that("lst_vat() measures match expected ranges", {
  skip_if_not_installed("superassp")
  skip_if(!voice_analysis_available(),
          "voice_analysis module not available")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- lst_vat(test_wav, verbose = FALSE)

  measures <- result$measures

  # DFA should be between 0 and 2 (typically 0.5-1.5 for speech)
  if ("DFA" %in% names(measures) && !is.null(measures$DFA)) {
    expect_true(measures$DFA >= 0 && measures$DFA <= 2)
  }

  # RPDE should be positive
  if ("RPDE" %in% names(measures) && !is.null(measures$RPDE)) {
    expect_true(measures$RPDE > 0)
  }

  # PPE should be positive
  if ("PPE" %in% names(measures) && !is.null(measures$PPE)) {
    expect_true(measures$PPE > 0)
  }

  # HNR measures should typically be positive (dB)
  hnr_measures <- names(measures)[grepl("^HNR", names(measures))]
  for (hnr_name in hnr_measures) {
    if (!is.null(measures[[hnr_name]])) {
      # HNR can be negative for very noisy signals, but typically > -10 dB
      expect_true(measures[[hnr_name]] > -20)
    }
  }
})

test_that("lst_vat() handles timeout parameter", {
  skip_if_not_installed("superassp")
  skip_if(!voice_analysis_available(),
          "voice_analysis module not available")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Test with reasonable timeout (should succeed)
  result <- lst_vat(test_wav,
                   timeout = 60,  # 60 seconds
                   verbose = FALSE)

  expect_true(result$success)
  expect_type(result$measures, "list")
})

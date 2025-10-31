# Test script for Phonet integration in superassp
# This tests the R reticulate wrapper for Phonet (both lst_ and trk_ functions)

library(testthat)
library(superassp)

test_that("Phonet availability checks work", {
  skip_if_not(phonet_available(), "Phonet not installed")

  # Test availability function
  expect_true(is.logical(phonet_available()))

  # Test info function
  info <- phonet_info()
  expect_s3_class(info, "phonet_info")
  expect_true("python_version" %in% names(info))
  expect_true("tensorflow_version" %in% names(info))
  expect_true("phonological_classes" %in% names(info))
  expect_equal(info$num_classes, 17)
})

test_that("lst_phonet processes audio files", {
  skip_if_not(phonet_available(), "Phonet not installed")

  # Find test audio file
  test_audio <- system.file("extdata", "test.wav", package = "superassp")
  if (!file.exists(test_audio)) {
    skip("Test audio file not found")
  }

  # Test with specific classes
  result <- lst_phonet(
    test_audio,
    classes = c("vocalic", "consonantal"),
    verbose = FALSE
  )

  # Check result structure
  expect_type(result, "list")
  expect_true("time" %in% names(result))
  expect_true("phoneme" %in% names(result))
  expect_true("vocalic" %in% names(result))
  expect_true("consonantal" %in% names(result))
  expect_true("file" %in% names(result))

  # Check data types and dimensions
  expect_true(is.numeric(result$time))
  expect_true(is.character(result$phoneme))
  expect_true(is.numeric(result$vocalic))
  expect_true(is.numeric(result$consonantal))
  expect_equal(length(result$time), length(result$phoneme))
  expect_equal(length(result$time), length(result$vocalic))

  # Check posterior values are probabilities (0-1)
  expect_true(all(result$vocalic >= 0 & result$vocalic <= 1))
  expect_true(all(result$consonantal >= 0 & result$consonantal <= 1))
})

test_that("lst_phonet handles all phonological classes", {
  skip_if_not(phonet_available(), "Phonet not installed")

  # Find test audio file
  test_audio <- system.file("extdata", "test.wav", package = "superassp")
  if (!file.exists(test_audio)) {
    skip("Test audio file not found")
  }

  # Test with "all" classes
  result <- lst_phonet(
    test_audio,
    classes = "all",
    verbose = FALSE
  )

  # Check that all 17 phonological classes are present
  expected_classes <- c(
    "vocalic", "consonantal", "back", "anterior", "open", "close",
    "nasal", "stop", "continuant", "lateral", "flap", "trill", "voice",
    "strident", "labial", "dental", "velar", "pause"
  )

  for (class_name in expected_classes) {
    expect_true(class_name %in% names(result),
                info = paste("Missing phonological class:", class_name))
    expect_true(is.numeric(result[[class_name]]))
    expect_equal(length(result[[class_name]]), length(result$time))
  }
})

test_that("lst_phonet handles time windowing", {
  skip_if_not(phonet_available(), "Phonet not installed")

  # Find test audio file
  test_audio <- system.file("extdata", "test.wav", package = "superassp")
  if (!file.exists(test_audio)) {
    skip("Test audio file not found")
  }

  # Get full file duration
  audio_info <- av::av_media_info(test_audio)
  duration <- as.numeric(audio_info$duration)

  if (duration < 2.0) {
    skip("Audio file too short for windowing test")
  }

  # Extract middle section
  result_windowed <- lst_phonet(
    test_audio,
    classes = c("nasal"),
    beginTime = 0.5,
    endTime = 1.5,
    verbose = FALSE
  )

  # Extract full file
  result_full <- lst_phonet(
    test_audio,
    classes = c("nasal"),
    verbose = FALSE
  )

  # Windowed result should be shorter
  expect_lt(length(result_windowed$time), length(result_full$time))

  # Windowed result should have roughly 100 frames (1.0 sec * 100 Hz)
  expect_gt(length(result_windowed$time), 80)  # Allow some tolerance
  expect_lt(length(result_windowed$time), 120)
})

test_that("lst_phonet handles multiple files", {
  skip_if_not(phonet_available(), "Phonet not installed")

  # Find test audio files
  test_audio1 <- system.file("extdata", "test.wav", package = "superassp")
  test_audio2 <- system.file("extdata", "test2.wav", package = "superassp")

  if (!file.exists(test_audio1)) {
    skip("Test audio files not found")
  }

  # Create second file if needed (duplicate first)
  if (!file.exists(test_audio2)) {
    test_audio2 <- tempfile(fileext = ".wav")
    file.copy(test_audio1, test_audio2)
    on.exit(unlink(test_audio2), add = TRUE)
  }

  # Process multiple files
  results <- lst_phonet(
    c(test_audio1, test_audio2),
    classes = c("vocalic"),
    verbose = FALSE
  )

  # Check results structure
  expect_type(results, "list")
  expect_equal(length(results), 2)

  # Check each result
  for (result in results) {
    expect_type(result, "list")
    expect_true("time" %in% names(result))
    expect_true("vocalic" %in% names(result))
    expect_true("file" %in% names(result))
  }
})

test_that("lst_phonet validates input parameters", {
  skip_if_not(phonet_available(), "Phonet not installed")

  # Test invalid phonological class
  expect_error(
    lst_phonet("dummy.wav", classes = c("invalid_class")),
    "Invalid phonological class"
  )

  # Test non-existent file
  expect_error(
    lst_phonet("nonexistent_file.wav", classes = c("vocalic")),
    "Unable to find the sound file"
  )
})

test_that("trk_phonet returns SSFF track objects", {
  skip_if_not(phonet_available(), "Phonet not installed")

  # Find test audio file
  test_audio <- system.file("extdata", "test.wav", package = "superassp")
  if (!file.exists(test_audio)) {
    skip("Test audio file not found")
  }

  # Test with toFile=FALSE (returns AsspDataObj)
  result <- trk_phonet(
    test_audio,
    classes = c("vocalic", "nasal"),
    toFile = FALSE,
    verbose = FALSE
  )

  # Check result is AsspDataObj
  expect_s3_class(result, "AsspDataObj")

  # Check tracks are present
  expect_true("vocalic" %in% names(result))
  expect_true("nasal" %in% names(result))

  # Check track data types
  expect_true(is.matrix(result$vocalic))
  expect_true(is.matrix(result$nasal))

  # Check metadata
  expect_true(attr(result, "sampleRate") == 100)  # 10ms frames
  expect_equal(attr(result, "trackFormats"), c("REAL32", "REAL32"))
  expect_true(attr(result, "origFreq") == 16000)  # Resampled to 16kHz

  # Check data dimensions match
  expect_equal(nrow(result$vocalic), nrow(result$nasal))
  expect_equal(nrow(result$vocalic), attr(result, "endRecord"))
})

test_that("trk_phonet writes SSFF files", {
  skip_if_not(phonet_available(), "Phonet not installed")

  # Find test audio file
  test_audio <- system.file("extdata", "test.wav", package = "superassp")
  if (!file.exists(test_audio)) {
    skip("Test audio file not found")
  }

  # Create temporary output directory
  temp_dir <- tempdir()

  # Process with toFile=TRUE
  num_processed <- trk_phonet(
    test_audio,
    classes = c("stop", "vocalic"),
    outputDirectory = temp_dir,
    toFile = TRUE,
    verbose = FALSE
  )

  # Check return value
  expect_equal(num_processed, 1)

  # Check SSFF file was created
  ssff_file <- file.path(temp_dir, sub("\\.wav$", ".phn", basename(test_audio)))
  expect_true(file.exists(ssff_file))

  # Read back SSFF file
  if (file.exists(ssff_file)) {
    result <- read.AsspDataObj(ssff_file)
    expect_s3_class(result, "AsspDataObj")
    expect_true("stop" %in% names(result))
    expect_true("vocalic" %in% names(result))

    # Clean up
    unlink(ssff_file)
  }
})

test_that("trk_phonet and lst_phonet produce consistent results", {
  skip_if_not(phonet_available(), "Phonet not installed")

  # Find test audio file
  test_audio <- system.file("extdata", "test.wav", package = "superassp")
  if (!file.exists(test_audio)) {
    skip("Test audio file not found")
  }

  # Get results from both functions
  lst_result <- lst_phonet(
    test_audio,
    classes = c("nasal"),
    verbose = FALSE
  )

  trk_result <- trk_phonet(
    test_audio,
    classes = c("nasal"),
    toFile = FALSE,
    verbose = FALSE
  )

  # Both should have same number of frames
  expect_equal(length(lst_result$time), nrow(trk_result$nasal))

  # Posteriors should be identical (or very close due to rounding)
  expect_equal(
    lst_result$nasal,
    as.vector(trk_result$nasal),
    tolerance = 0.001
  )
})

test_that("trk_phonet handles all phonological classes", {
  skip_if_not(phonet_available(), "Phonet not installed")

  # Find test audio file
  test_audio <- system.file("extdata", "test.wav", package = "superassp")
  if (!file.exists(test_audio)) {
    skip("Test audio file not found")
  }

  # Track all classes
  result <- trk_phonet(
    test_audio,
    classes = "all",
    toFile = FALSE,
    verbose = FALSE
  )

  # Should have 18 tracks (all phonological classes)
  expect_gte(length(names(result)), 18)

  # Check some expected classes
  expected_classes <- c(
    "vocalic", "consonantal", "nasal", "stop",
    "voice", "strident", "labial", "pause"
  )

  for (cls in expected_classes) {
    expect_true(cls %in% names(result),
                info = paste("Missing class:", cls))
  }
})

# Print summary message
message("\n========================================")
message("Phonet Integration Test Summary")
message("========================================")
if (phonet_available()) {
  info <- phonet_info()
  message("✓ Phonet is available")
  message("  Python: ", info$python_version)
  message("  TensorFlow: ", info$tensorflow_version)
  message("  Phonological classes: ", info$num_classes)
  message("  Functions: lst_phonet(), trk_phonet()")
} else {
  message("✗ Phonet is not installed")
  message("  Install with: install_phonet()")
}
message("========================================\n")

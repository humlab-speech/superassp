test_that("prep_recode works with single WAV file (no conversion)", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- prep_recode(test_wav, format = "wav", verbose = FALSE)

  expect_type(result, "integer")
  expect_true(length(result) > 0)

  # Should have av::read_audio_bin attributes
  expect_true(!is.null(attr(result, "channels")))
  expect_true(!is.null(attr(result, "sample_rate")))

  # Attributes should be numeric/integer
  expect_type(attr(result, "channels"), "integer")
  expect_type(attr(result, "sample_rate"), "integer")
})

test_that("prep_recode validates format argument", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Missing format
  expect_error(
    prep_recode(test_wav, verbose = FALSE),
    "format argument is required"
  )

  # NULL format
  expect_error(
    prep_recode(test_wav, format = NULL, verbose = FALSE),
    "format argument is required"
  )

  # Empty format
  expect_error(
    prep_recode(test_wav, format = "", verbose = FALSE),
    "format argument is required"
  )
})

test_that("prep_recode handles missing files gracefully", {
  expect_warning(
    result <- prep_recode("nonexistent.wav", format = "wav", verbose = FALSE),
    "File not found"
  )

  expect_null(result)
})

test_that("prep_recode supports time windowing", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get full file
  result_full <- prep_recode(test_wav, format = "wav", verbose = FALSE)

  # Get segment
  result_window <- prep_recode(test_wav,
                               format = "wav",
                               start_time = 0.1,
                               end_time = 0.5,
                               verbose = FALSE)

  expect_type(result_window, "integer")

  # Windowed result should be shorter
  expect_true(length(result_window) < length(result_full))

  # Check attributes preserved
  expect_equal(attr(result_window, "channels"), attr(result_full, "channels"))
  expect_equal(attr(result_window, "sample_rate"), attr(result_full, "sample_rate"))
})

test_that("prep_recode supports sample rate conversion", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Convert to 16kHz
  result_16k <- prep_recode(test_wav,
                           format = "wav",
                           sample_rate = 16000,
                           verbose = FALSE)

  expect_type(result_16k, "integer")
  expect_equal(attr(result_16k, "sample_rate"), 16000L)
})

test_that("prep_recode supports channel conversion", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Convert to mono (even if already mono)
  result_mono <- prep_recode(test_wav,
                             format = "wav",
                             channels = 1,
                             verbose = FALSE)

  expect_type(result_mono, "integer")
  expect_equal(attr(result_mono, "channels"), 1L)
})

test_that("prep_recode batch processing works", {
  test_files <- list.files(
    system.file("samples", "sustained", package = "superassp"),
    pattern = "\\.wav$",
    full.names = TRUE
  )

  skip_if(length(test_files) < 2, "Need at least 2 test files")
  test_files <- test_files[1:2]

  results <- prep_recode(test_files, format = "wav", verbose = FALSE)

  expect_type(results, "list")
  expect_length(results, 2)

  # Both should be integer vectors
  expect_type(results[[1]], "integer")
  expect_type(results[[2]], "integer")

  # Both should have attributes
  expect_true(!is.null(attr(results[[1]], "channels")))
  expect_true(!is.null(attr(results[[2]], "channels")))
})

test_that("prep_recode returns same format as av::read_audio_bin", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Direct av::read_audio_bin
  direct <- av::read_audio_bin(test_wav)

  # Via prep_recode (no conversion)
  recoded <- prep_recode(test_wav, format = "wav", verbose = FALSE)

  # Should have same type
  expect_equal(typeof(direct), typeof(recoded))

  # Should have same attributes
  expect_equal(names(attributes(direct)), names(attributes(recoded)))
  expect_equal(attr(direct, "channels"), attr(recoded, "channels"))
  expect_equal(attr(direct, "sample_rate"), attr(recoded, "sample_rate"))

  # Should have same length
  expect_equal(length(direct), length(recoded))

  # Values should be identical (no conversion)
  expect_equal(direct, recoded)
})

test_that("prep_recode with custom bit rate", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Convert with explicit bit rate
  # Note: This triggers re-encoding
  result <- prep_recode(test_wav,
                       format = "wav",
                       bit_rate = 128000,
                       verbose = FALSE)

  # Bit rate parameter triggers re-encoding, so should succeed
  expect_type(result, "integer")
  expect_true(!is.null(attr(result, "channels")))
  expect_true(!is.null(attr(result, "sample_rate")))
})

test_that("prep_recode combined parameters", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Combine sample rate + time window + channels
  result <- prep_recode(test_wav,
                       format = "wav",
                       sample_rate = 16000,
                       start_time = 0.1,
                       end_time = 0.5,
                       channels = 1,
                       verbose = FALSE)

  expect_type(result, "integer")
  expect_equal(attr(result, "sample_rate"), 16000L)
  expect_equal(attr(result, "channels"), 1L)

  # Should be short duration
  duration <- length(result) / attr(result, "channels") / attr(result, "sample_rate")
  expect_true(duration < 0.5)
})

test_that("prep_recode handles files without audio", {
  # Create a text file (no audio)
  temp_txt <- tempfile(fileext = ".txt")
  writeLines("test", temp_txt)
  on.exit(unlink(temp_txt))

  expect_warning(
    result <- prep_recode(temp_txt, format = "wav", verbose = FALSE),
    "Invalid media file"
  )

  expect_null(result)
})

test_that("prep_recode calculates correct duration", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get audio
  result <- prep_recode(test_wav, format = "wav", verbose = FALSE)

  # Calculate duration
  n_samples <- length(result)
  channels <- attr(result, "channels")
  sample_rate <- attr(result, "sample_rate")
  duration <- n_samples / channels / sample_rate

  # Duration should be positive and reasonable
  expect_true(duration > 0)
  expect_true(duration < 60)  # Less than 1 minute for test files
})

test_that("prep_recode with start_time only", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get full file
  result_full <- prep_recode(test_wav, format = "wav", verbose = FALSE)

  # Get from 0.5 seconds onwards
  result_start <- prep_recode(test_wav,
                              format = "wav",
                              start_time = 0.5,
                              verbose = FALSE)

  expect_type(result_start, "integer")

  # Should be shorter than full file
  expect_true(length(result_start) < length(result_full))
})

test_that("prep_recode with end_time only", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get full file
  result_full <- prep_recode(test_wav, format = "wav", verbose = FALSE)

  # Get first 0.5 seconds
  result_end <- prep_recode(test_wav,
                            format = "wav",
                            end_time = 0.5,
                            verbose = FALSE)

  expect_type(result_end, "integer")

  # Should be shorter than full file
  expect_true(length(result_end) < length(result_full))
})

test_that("prep_recode preserves stereo when requested", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Convert to stereo (duplicate mono to stereo)
  result_stereo <- prep_recode(test_wav,
                               format = "wav",
                               channels = 2,
                               verbose = FALSE)

  expect_type(result_stereo, "integer")
  expect_equal(attr(result_stereo, "channels"), 2L)

  # Should have double the samples (interleaved)
  # But we can't compare directly since mono->stereo duplicates
  expect_true(length(result_stereo) > 0)
})

test_that("prep_recode batch with mixed success", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  files <- c(test_wav, "nonexistent.wav", test_wav)

  expect_warning(
    results <- prep_recode(files, format = "wav", verbose = FALSE),
    "File not found"
  )

  expect_type(results, "list")
  expect_length(results, 3)

  # First and third should succeed
  expect_type(results[[1]], "integer")
  expect_null(results[[2]])  # Failed
  expect_type(results[[3]], "integer")
})

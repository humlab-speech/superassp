test_that("yin_cpp works with test audio", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Test basic functionality (toFile = FALSE)
  result <- trk_yin(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("F0" %in% names(result))
  expect_true("prob" %in% names(result))
  expect_true(nrow(result$F0) > 0)
})

test_that("yin_cpp respects F0 range", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_yin(test_wav, minF = 100, maxF = 300, toFile = FALSE, verbose = FALSE)

  # Check that F0 values are within range (excluding zeros/unvoiced)
  f0_values <- result$F0[result$F0 > 0]
  if (length(f0_values) > 0) {
    expect_true(all(f0_values >= 100 | f0_values == 0))
    expect_true(all(f0_values <= 300 | f0_values == 0))
  }
})

test_that("yin_cpp handles time windowing", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get audio duration first
  audio_info <- av::av_media_info(test_wav)
  duration <- as.numeric(audio_info$duration)

  if (duration > 0.5) {
    # Test windowing
    result <- trk_yin(test_wav, beginTime = 0.1, endTime = 0.5,
                      toFile = FALSE, verbose = FALSE)

    expect_s3_class(result, "AsspDataObj")
    expect_true(nrow(result$F0) > 0)

    # Start time should be close to 0.1
    start_time <- attr(result, "startTime")
    expect_true(start_time >= 0.0 && start_time <= 0.2)
  }
})

test_that("yin_cpp handles file output", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Create temp directory
  temp_dir <- tempdir()

  # Process with toFile = TRUE
  result_file <- trk_yin(test_wav, toFile = TRUE, explicitExt = "yip",
                          outputDirectory = temp_dir, verbose = FALSE)

  expect_type(result_file, "character")
  expect_true(file.exists(result_file))

  # Clean up
  unlink(result_file)
})

test_that("yin_cpp handles multiple files", {
  skip_if_not_installed("superassp")

  test_files <- c(
    system.file("samples", "sustained", "a1.wav", package = "superassp"),
    system.file("samples", "sustained", "a32b.wav", package = "superassp")
  )

  skip_if(any(test_files == ""), "Test files not found")
  skip_if(length(test_files[file.exists(test_files)]) < 2, "Need 2 test files")

  # Create temp directory
  temp_dir <- tempdir()

  # Process multiple files
  result_files <- trk_yin(test_files, toFile = TRUE, explicitExt = "yip",
                           outputDirectory = temp_dir, verbose = FALSE)

  expect_type(result_files, "list")
  expect_length(result_files, 2)
  expect_true(all(sapply(result_files, file.exists)))

  # Clean up
  invisible(lapply(result_files, unlink))
})

test_that("yin_cpp handles non-WAV formats via av", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("av")

  # Try to find MP3 test file
  test_mp3 <- system.file("samples", "sustained", "a7.mp3", package = "superassp")

  # Skip if no MP3 available
  skip_if(test_mp3 == "" || !file.exists(test_mp3), "MP3 test file not available")

  result <- trk_yin(test_mp3, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("F0" %in% names(result))
  expect_true("prob" %in% names(result))
})

test_that("yin_cpp validates input parameters", {
  skip_if_not_installed("superassp")

  # Test missing file (S7 dispatch throws method lookup error for NULL)
  expect_error(trk_yin(NULL), "(No input files|Can't find method)")

  # Test non-existent file
  expect_error(trk_yin("nonexistent.wav"), "do not exist")
})

test_that("yin_cpp returns probability track", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_yin(test_wav, toFile = FALSE, verbose = FALSE)

  # Check probability track exists and has valid values [0, 1]
  expect_true("prob" %in% names(result))
  prob_values <- result$prob[,1]
  expect_true(all(prob_values >= 0.0 & prob_values <= 1.0))
})

test_that("yin_cpp custom threshold affects results", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Strict threshold (higher = fewer voiced frames)
  result_strict <- trk_yin(test_wav, threshold = 0.05, toFile = FALSE, verbose = FALSE)

  # Permissive threshold (higher = more voiced frames)
  result_permissive <- trk_yin(test_wav, threshold = 0.2, toFile = FALSE, verbose = FALSE)

  # Count voiced frames (F0 > 0)
  voiced_strict <- sum(result_strict$F0 > 0)
  voiced_permissive <- sum(result_permissive$F0 > 0)

  # Permissive should generally have more or equal voiced frames
  expect_true(voiced_permissive >= voiced_strict * 0.8)  # Allow some variance
})

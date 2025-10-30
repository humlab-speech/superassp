test_that("pyin_cpp works with test audio", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Test basic functionality (toFile = FALSE)
  result <- trk_pyin(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("F0" %in% names(result))
  expect_true("prob" %in% names(result))
  expect_true(nrow(result$F0) > 0)
})

test_that("pyin_cpp respects F0 range", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_pyin(test_wav, minF = 100, maxF = 300, toFile = FALSE, verbose = FALSE)

  # Check that F0 values are within range (excluding zeros/unvoiced)
  f0_values <- result$F0[result$F0 > 0]
  if (length(f0_values) > 0) {
    expect_true(all(f0_values >= 100 | f0_values == 0))
    expect_true(all(f0_values <= 300 | f0_values == 0))
  }
})

test_that("pyin_cpp handles file output", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Create temp directory
  temp_dir <- tempdir()

  # Process with toFile = TRUE
  result_file <- trk_pyin(test_wav, toFile = TRUE, explicitExt = "pyp",
                           outputDirectory = temp_dir, verbose = FALSE)

  expect_type(result_file, "character")
  expect_true(file.exists(result_file))

  # Clean up
  unlink(result_file)
})

test_that("pyin_cpp produces similar results to yin_cpp", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get results from both algorithms with same parameters
  result_yin <- trk_yin(test_wav, minF = 70, maxF = 200, windowShift = 5,
                         windowSize = 30, threshold = 0.1,
                         toFile = FALSE, verbose = FALSE)

  result_pyin <- trk_pyin(test_wav, minF = 70, maxF = 200, windowShift = 5,
                           windowSize = 30, threshold = 0.1,
                           toFile = FALSE, verbose = FALSE)

  # Results should be similar (simplified pYIN uses same algorithm as YIN)
  # Check that at least 90% of frames agree
  f0_yin <- result_yin$F0[,1]
  f0_pyin <- result_pyin$F0[,1]

  # Frames should have same length
  expect_equal(length(f0_yin), length(f0_pyin))

  # Voiced/unvoiced decisions should be very similar
  voiced_yin <- f0_yin > 0
  voiced_pyin <- f0_pyin > 0
  agreement <- sum(voiced_yin == voiced_pyin) / length(voiced_yin)
  expect_true(agreement >= 0.90)
})

test_that("pyin_cpp handles time windowing", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get audio duration first
  audio_info <- av::av_media_info(test_wav)
  duration <- as.numeric(audio_info$duration)

  if (duration > 0.5) {
    # Test windowing
    result <- trk_pyin(test_wav, beginTime = 0.1, endTime = 0.5,
                        toFile = FALSE, verbose = FALSE)

    expect_s3_class(result, "AsspDataObj")
    expect_true(nrow(result$F0) > 0)

    # Start time should be close to 0.1
    start_time <- attr(result, "startTime")
    expect_true(start_time >= 0.0 && start_time <= 0.2)
  }
})

test_that("pyin_cpp validates input parameters", {
  skip_if_not_installed("superassp")

  # Test missing file
  expect_error(trk_pyin(NULL), "No input files")

  # Test non-existent file
  expect_error(trk_pyin("nonexistent.wav"), "do not exist")
})

test_that("pyin_cpp returns probability track", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_pyin(test_wav, toFile = FALSE, verbose = FALSE)

  # Check probability track exists and has valid values [0, 1]
  expect_true("prob" %in% names(result))
  prob_values <- result$prob[,1]
  expect_true(all(prob_values >= 0.0 & prob_values <= 1.0))
})

test_that("pyin_cpp handles multiple files", {
  skip_if_not_installed("superassp")

  test_files <- c(
    system.file("samples", "sustained", "a1.wav", package = "superassp"),
    system.file("samples", "sustained", "a2.wav", package = "superassp")
  )

  skip_if(any(test_files == ""), "Test files not found")
  skip_if(length(test_files[file.exists(test_files)]) < 2, "Need 2 test files")

  # Create temp directory
  temp_dir <- tempdir()

  # Process multiple files
  result_files <- trk_pyin(test_files, toFile = TRUE, explicitExt = "pyp",
                            outputDirectory = temp_dir, verbose = FALSE)

  expect_type(result_files, "list")
  expect_length(result_files, 2)
  expect_true(all(sapply(result_files, file.exists)))

  # Clean up
  invisible(lapply(result_files, unlink))
})

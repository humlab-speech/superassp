test_that("trk_tandem works with single WAV file", {
  skip_if_not_installed("superassp")
  
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  
  result <- trk_tandem(test_wav, toFile = FALSE, verbose = FALSE)
  
  expect_s3_class(result, "AsspDataObj")
  expect_true("pitch" %in% names(result))
  expect_true("voicing_prob" %in% names(result))
  expect_true(length(result$pitch) > 0)
  expect_true(all(result$voicing_prob >= 0 & result$voicing_prob <= 1, na.rm = TRUE))
})

test_that("trk_tandem handles custom F0 range", {
  skip_if_not_installed("superassp")
  
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  
  result <- trk_tandem(test_wav, minF = 80, maxF = 400, 
                       toFile = FALSE, verbose = FALSE)
  
  expect_s3_class(result, "AsspDataObj")
  pitch_values <- result$pitch[!is.na(result$pitch)]
  
  # Check that pitched frames are within range (allowing some tolerance)
  if (length(pitch_values) > 0) {
    expect_true(min(pitch_values) >= 70)  # Allow 10 Hz tolerance
    expect_true(max(pitch_values) <= 410)
  }
})

test_that("trk_tandem handles non-WAV formats via av", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("av")
  
  # Try to find MP3 test file
  test_mp3 <- system.file("samples", "test.mp3", package = "superassp")
  if (test_mp3 == "") {
    skip("MP3 test file not available")
  }
  
  result <- trk_tandem(test_mp3, toFile = FALSE, verbose = FALSE)
  
  expect_s3_class(result, "AsspDataObj")
  expect_true("pitch" %in% names(result))
})

test_that("trk_tandem handles resampling", {
  skip_if_not_installed("superassp")
  skip_if_not_installed("av")
  
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  
  # TANDEM requires 20 kHz - function should handle resampling
  result <- trk_tandem(test_wav, target_sample_rate = 20000, 
                       toFile = FALSE, verbose = FALSE)
  
  expect_s3_class(result, "AsspDataObj")
  expect_equal(attr(result, "sampleRate"), 100)  # 100 Hz frame rate
})

test_that("trk_tandem batch processing works", {
  skip_if_not_installed("superassp")
  
  test_files <- c(
    system.file("samples", "sustained", "a1.wav", package = "superassp"),
    system.file("samples", "sustained", "a2.wav", package = "superassp")
  )
  
  # Filter out missing files
  test_files <- test_files[file.exists(test_files)]
  skip_if(length(test_files) < 2, "Need at least 2 test files")
  
  results <- trk_tandem(test_files, toFile = FALSE, verbose = FALSE)
  
  expect_type(results, "list")
  expect_equal(length(results), length(test_files))
  expect_true(all(sapply(results, function(x) "AsspDataObj" %in% class(x))))
})

test_that("trk_tandem file output works", {
  skip_if_not_installed("superassp")
  
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  
  temp_dir <- tempdir()
  output_path <- trk_tandem(test_wav, toFile = TRUE, 
                            outputDirectory = temp_dir, 
                            verbose = FALSE)
  
  expect_type(output_path, "character")
  expect_true(file.exists(output_path))
  expect_true(grepl("\\.tnd$", output_path))
  
  # Cleanup
  unlink(output_path)
})

test_that("trk_tandem attributes are correct", {
  # Get the character method which has the attributes
  char_method <- tryCatch({
    S7::method(trk_tandem, S7::class_character)
  }, error = function(e) {
    trk_tandem
  })
  
  expect_equal(attr(char_method, "ext"), "tnd")
  expect_equal(attr(char_method, "tracks"), c("pitch", "voicing_prob"))
  expect_equal(attr(char_method, "outputType"), "SSFF")
  expect_equal(attr(char_method, "nativeFiletypes"), c("wav"))
})

test_that("trk_tandem error handling works", {
  expect_error(trk_tandem("nonexistent_file.wav"), "not found")
})

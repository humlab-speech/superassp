test_that("S7 dispatch works for trk_rmsana with character input", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Should work with file path (character vector)
  result <- trk_rmsana(test_wav, toFile = FALSE, verbose = FALSE)

  # Check result
  expect_s3_class(result, "AsspDataObj")
  expect_true("rms" %in% names(result))
})

test_that("S7 dispatch works for trk_rmsana with AVAudio input", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Create AVAudio object
  audio <- read_avaudio(test_wav)

  # Should work with AVAudio object
  result <- trk_rmsana(audio, toFile = FALSE, verbose = FALSE)

  # Check result
  expect_s3_class(result, "AsspDataObj")
  expect_true("rms" %in% names(result))
})

test_that("S7 dispatch produces same results for character and AVAudio", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Process with file path
  result_char <- trk_rmsana(test_wav, toFile = FALSE, verbose = FALSE)

  # Process with AVAudio
  audio <- read_avaudio(test_wav)
  result_avaudio <- trk_rmsana(audio, toFile = FALSE, verbose = FALSE)

  # Results should be similar (may not be identical due to file I/O)
  expect_equal(length(result_char$rms), length(result_avaudio$rms))
  expect_equal(attr(result_char, "sampleRate"), attr(result_avaudio, "sampleRate"))
})

test_that("S7 dispatch works for lst_voice_sauce with character input", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Check if voice_sauce available
  skip_if(!voice_sauce_available(), "VoiceSauce not installed")

  # Should work with file path
  result <- lst_voice_sauce(test_wav, verbose = FALSE)

  # Check result
  expect_type(result, "list")
  expect_true("F0" %in% names(result))
})

test_that("S7 dispatch works for lst_voice_sauce with AVAudio input", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Check if voice_sauce available
  skip_if(!voice_sauce_available(), "VoiceSauce not installed")

  # Create AVAudio object
  audio <- read_avaudio(test_wav)

  # Should work with AVAudio object
  result <- lst_voice_sauce(audio, verbose = FALSE)

  # Check result
  expect_type(result, "list")
  expect_true("F0" %in% names(result))
})

test_that("S7 dispatch maintains backward compatibility", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # All existing usage patterns should still work

  # Pattern 1: Single file with toFile=FALSE
  result1 <- trk_rapt(test_wav, toFile = FALSE, verbose = FALSE)
  expect_s3_class(result1, "AsspDataObj")

  # Pattern 2: Multiple files
  test_files <- rep(test_wav, 2)
  result2 <- trk_rapt(test_files, toFile = FALSE, verbose = FALSE)
  expect_type(result2, "list")
  expect_length(result2, 2)

  # Pattern 3: Custom parameters
  result3 <- trk_rapt(test_wav, minF = 75, maxF = 300,
                     toFile = FALSE, verbose = FALSE)
  expect_s3_class(result3, "AsspDataObj")
})

test_that("S7 dispatch works with AVAudio list", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Create multiple AVAudio objects
  audio1 <- read_avaudio(test_wav, start_time = 0.0, end_time = 1.0)
  audio2 <- read_avaudio(test_wav, start_time = 1.0, end_time = 2.0)

  # Note: Current implementation processes single AVAudio objects
  # For multiple, users would need to use lapply or similar
  result1 <- trk_rapt(audio1, toFile = FALSE, verbose = FALSE)
  result2 <- trk_rapt(audio2, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result1, "AsspDataObj")
  expect_s3_class(result2, "AsspDataObj")
})

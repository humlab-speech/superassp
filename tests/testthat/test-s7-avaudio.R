test_that("AVAudio class can be created", {
  # Skip if test file not found
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Create AVAudio from file
  audio <- read_avaudio(test_wav)

  # Check class
  expect_true(is_avaudio(audio))
  # S7 objects have package-qualified class names
  expect_true(inherits(audio, "superassp::AVAudio"))

  # Check properties
  expect_true(audio@sample_rate > 0)
  expect_true(audio@channels > 0)
  expect_true(length(audio@samples) > 0)
})

test_that("AVAudio can be created from prep_recode output", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get audio data from prep_recode
  audio_data <- prep_recode(test_wav, format = "wav", verbose = FALSE)

  # Convert to AVAudio
  audio <- as_avaudio(audio_data, file_path = test_wav)

  # Check it worked
  expect_true(is_avaudio(audio))
  expect_equal(audio@file_path, test_wav)
})

test_that("AVAudio can be converted back to av format", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Create AVAudio
  audio <- read_avaudio(test_wav)

  # Convert back
  audio_vec <- avaudio_to_av(audio)

  # Check attributes
  expect_type(audio_vec, "integer")
  expect_equal(attr(audio_vec, "sample_rate"), audio@sample_rate)
  expect_equal(attr(audio_vec, "channels"), audio@channels)
  expect_equal(length(audio_vec), length(audio@samples))
})

test_that("AVAudio can be converted to temp file", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Create AVAudio
  audio <- read_avaudio(test_wav)

  # Convert to temp file
  temp_file <- avaudio_to_tempfile(audio, verbose = FALSE)

  # Check file exists
  expect_true(file.exists(temp_file))

  # Check it's a WAV file
  expect_match(temp_file, "\\.wav$")

  # Clean up
  unlink(temp_file)
})

test_that("AVAudio print method works", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio <- read_avaudio(test_wav)

  # Should not error
  expect_output(print(audio), "<AVAudio>")
  expect_output(print(audio), "Sample rate:")
})

test_that("AVAudio summary method works", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio <- read_avaudio(test_wav)

  # Should not error
  expect_output(summary(audio), "<AVAudio Summary>")
  expect_output(summary(audio), "statistics")
})

test_that("AVAudio with time windowing works", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get full file
  audio_full <- read_avaudio(test_wav)

  # Get windowed version
  audio_window <- read_avaudio(test_wav, start_time = 0.1, end_time = 0.5)

  # Windowed should be shorter
  expect_true(length(audio_window@samples) < length(audio_full@samples))
})

test_that("AVAudio with resampling works", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Resample to 16kHz
  audio_16k <- read_avaudio(test_wav, sample_rate = 16000)

  # Check sample rate
  expect_equal(audio_16k@sample_rate, 16000L)
})

test_that("AVAudio validation works", {
  # Should error with invalid data
  expect_error(
    AVAudio(samples = integer(0), sample_rate = 0L, channels = 1L),
    "sample_rate must be positive"
  )

  expect_error(
    AVAudio(samples = integer(10), sample_rate = 16000L, channels = 0L),
    "channels must be positive"
  )

  expect_error(
    AVAudio(samples = integer(11), sample_rate = 16000L, channels = 2L),
    "multiple of channels"
  )
})

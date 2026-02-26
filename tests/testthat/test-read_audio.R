test_that("read_audio reads native WAV via ASSP C path", {
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wav == "", "test wav not found")

  result <- read_audio(wav)
  expect_s3_class(result, "AsspDataObj")
  expect_true("audio" %in% names(result))
  expect_true(attr(result, "sampleRate") > 0)
})

test_that("read_audio falls back to av for MP3", {
  skip_if_not_installed("av")
  mp3 <- system.file("samples", "sustained", "a7.mp3", package = "superassp")
  skip_if(mp3 == "" || !file.exists(mp3), "test mp3 not found")

  result <- read_audio(mp3)
  expect_s3_class(result, "AsspDataObj")
  expect_true(attr(result, "sampleRate") > 0)
})

test_that("read_audio respects begin/end (seconds)", {
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wav == "", "test wav not found")

  full <- read_audio(wav)
  part <- read_audio(wav, begin = 0.1, end = 0.3)
  expect_lt(numRecs.AsspDataObj(part), numRecs.AsspDataObj(full))
})

test_that("read_audio samples=TRUE approximates correctly for WAV", {
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wav == "", "test wav not found")

  full <- read_audio(wav)
  sr   <- attr(full, "sampleRate")
  part <- read_audio(wav, begin = round(0.1 * sr), end = round(0.3 * sr), samples = TRUE)
  expect_lt(numRecs.AsspDataObj(part), numRecs.AsspDataObj(full))
})

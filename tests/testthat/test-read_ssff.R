test_that("read_ssff reads WAV file like read.AsspDataObj", {
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wav == "", "test wav not found")

  result <- read_ssff(wav)
  expect_s3_class(result, "AsspDataObj")
  expect_true("audio" %in% names(result))
  expect_true(attr(result, "sampleRate") > 0)
})

test_that("read_ssff respects begin/end windowing", {
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wav == "", "test wav not found")

  full  <- read_ssff(wav)
  part  <- read_ssff(wav, begin = 0.1, end = 0.3)
  expect_lt(numRecs.AsspDataObj(part), numRecs.AsspDataObj(full))
})

test_that("read_ssff samples=TRUE uses sample-based indexing", {
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wav == "", "test wav not found")

  full <- read_ssff(wav)
  sr   <- attr(full, "sampleRate")
  part <- read_ssff(wav, begin = round(0.1 * sr), end = round(0.3 * sr), samples = TRUE)
  expect_lt(numRecs.AsspDataObj(part), numRecs.AsspDataObj(full))
})

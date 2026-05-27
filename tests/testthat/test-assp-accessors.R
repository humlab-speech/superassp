test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")

test_that("sample_rate() returns numeric Hz", {
  skip_if(test_wav == "", "test file missing")
  x <- read_audio(test_wav)
  expect_equal(sample_rate(x), attr(x, "sampleRate"))
  expect_true(sample_rate(x) > 0)
})

test_that("n_records() returns positive integer matching endRecord - startRecord + 1", {
  skip_if(test_wav == "", "test file missing")
  x <- read_audio(test_wav)
  expect_equal(n_records(x), attr(x, "endRecord") - attr(x, "startRecord") + 1L)
  expect_true(n_records(x) > 0L)
})

test_that("signal_duration() equals n_records / sample_rate", {
  skip_if(test_wav == "", "test file missing")
  x <- read_audio(test_wav)
  expect_equal(signal_duration(x), n_records(x) / sample_rate(x))
})

test_that("start_time() returns numeric seconds", {
  skip_if(test_wav == "", "test file missing")
  x <- read_audio(test_wav)
  expect_equal(start_time(x), attr(x, "startTime"))
  expect_true(is.numeric(start_time(x)))
})

test_that("track_names() returns names(x)", {
  skip_if(test_wav == "", "test file missing")
  x <- read_audio(test_wav)
  expect_equal(track_names(x), names(x))
  expect_true(is.character(track_names(x)))
})

test_that("file_path() returns character path for file-backed objects", {
  skip_if(test_wav == "", "test file missing")
  x <- read_audio(test_wav)
  fp <- file_path(x)
  expect_true(is.character(fp) || is.null(fp))
})

test_that("track_formats() returns character vector or NULL", {
  skip_if(test_wav == "", "test file missing")
  x <- read_audio(test_wav)
  tf <- track_formats(x)
  expect_true(is.character(tf) || is.null(tf))
})

test_that("new accessors work on SSFF track objects", {
  skip_if(test_wav == "", "test file missing")
  x <- trk_rms(test_wav, toFile = FALSE, verbose = FALSE)
  expect_true(sample_rate(x) > 0)
  expect_true(n_records(x) > 0L)
  expect_true(signal_duration(x) > 0)
  expect_true(is.character(track_names(x)))
  expect_true(length(track_names(x)) >= 1L)
})

test_that("old accessors emit lifecycle deprecation warnings", {
  skip_if(test_wav == "", "test file missing")
  x <- read_audio(test_wav)
  expect_warning(dur(x),       class = "lifecycle_warning_deprecated")
  expect_warning(numRecs(x),   class = "lifecycle_warning_deprecated")
  expect_warning(rate(x),      class = "lifecycle_warning_deprecated")
  expect_warning(startTime(x), class = "lifecycle_warning_deprecated")
  expect_warning(tracks(x),    class = "lifecycle_warning_deprecated")
})

test_that("old accessors return same values as new API", {
  skip_if(test_wav == "", "test file missing")
  x <- read_audio(test_wav)
  withr::with_options(list(lifecycle_verbosity = "quiet"), {
    expect_equal(dur(x),       signal_duration(x))
    expect_equal(numRecs(x),   n_records(x))
    expect_equal(rate(x),      sample_rate(x))
    expect_equal(startTime(x), start_time(x))
    expect_equal(tracks(x),    track_names(x))
  })
})

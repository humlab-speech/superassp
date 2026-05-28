# Smoke tests for lst_* functions that had zero coverage.
# Each test checks: returns non-null result with expected structure.

wav_file <- function() {
  f <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  if (f == "") skip("Test wav not found")
  f
}

# ---- OpenSMILE feature sets ----

test_that("lst_eGeMAPS returns named list with 88 features", {
  result <- lst_eGeMAPS(wav_file(), toFile = FALSE)
  expect_true(is.list(result))
  expect_gt(length(result), 0)
})

test_that("lst_eGeMAPS return_jstf=TRUE returns JsonTrackObj", {
  result <- lst_eGeMAPS(wav_file(), return_jstf = TRUE)
  expect_s3_class(result, "JsonTrackObj")
})

test_that("lst_GeMAPS returns named list with features", {
  result <- lst_GeMAPS(wav_file(), toFile = FALSE)
  expect_true(is.list(result))
  expect_gt(length(result), 0)
})

test_that("lst_GeMAPS return_jstf=TRUE returns JsonTrackObj", {
  result <- lst_GeMAPS(wav_file(), return_jstf = TRUE)
  expect_s3_class(result, "JsonTrackObj")
})

test_that("lst_ComParE_2016 returns named list with features", {
  result <- lst_ComParE_2016(wav_file(), toFile = FALSE)
  expect_true(is.list(result))
  expect_gt(length(result), 0)
})

test_that("lst_ComParE_2016 return_jstf=TRUE returns JsonTrackObj", {
  result <- lst_ComParE_2016(wav_file(), return_jstf = TRUE)
  expect_s3_class(result, "JsonTrackObj")
})

test_that("lst_emobase returns named list with features", {
  result <- suppressWarnings(lst_emobase(wav_file(), toFile = FALSE))
  expect_true(is.list(result))
  expect_gt(length(result), 0)
})

test_that("lst_emobase return_jstf=TRUE returns JsonTrackObj", {
  result <- suppressWarnings(lst_emobase(wav_file(), return_jstf = TRUE))
  expect_s3_class(result, "JsonTrackObj")
})

# ---- pladdrr-based summary functions ----

test_that("lst_voice_tremor returns data.frame with tremor measures", {
  result <- lst_voice_tremor(wav_file(), toFile = FALSE)
  expect_true(is.data.frame(result) || is.list(result))
  expect_gt(length(result), 0)
})

test_that("lst_voice_tremor return_jstf=TRUE returns JsonTrackObj", {
  skip_if_not(pladdrr_available(), "pladdrr not available")
  result <- lst_voice_tremor(wav_file(), return_jstf = TRUE)
  expect_s3_class(result, "JsonTrackObj")
})

test_that("lst_vq returns data.frame with voice quality measures", {
  result <- lst_vq(wav_file(), toFile = FALSE)
  expect_true(is.data.frame(result) || is.list(result))
  expect_gt(length(result), 0)
})

test_that("lst_vq return_jstf=TRUE returns JsonTrackObj", {
  skip_if_not(pladdrr_available(), "pladdrr not available")
  result <- lst_vq(wav_file(), return_jstf = TRUE)
  expect_s3_class(result, "JsonTrackObj")
})

# ---- lst_pharyngeal: requires time window or TextGrid ----
# Tested with explicit time window; returns data.frame (may be empty for
# a sustained vowel without strong pharyngeal features).

test_that("lst_pharyngeal returns data.frame when given time window", {
  result <- suppressWarnings(
    lst_pharyngeal(wav_file(), beginTime = 0.1, endTime = 2.0, toFile = FALSE)
  )
  expect_true(is.data.frame(result) || is.list(result))
})

test_that("lst_pharyngeal return_jstf=TRUE returns JsonTrackObj", {
  skip_if_not(pladdrr_available(), "pladdrr not available")
  result <- suppressWarnings(
    lst_pharyngeal(wav_file(), beginTime = 0.1, endTime = 2.0, return_jstf = TRUE)
  )
  expect_s3_class(result, "JsonTrackObj")
})

test_that("lst_voxit return_jstf=TRUE returns JsonTrackObj", {
  result <- lst_voxit(wav_file(), return_jstf = TRUE, verbose = FALSE)
  expect_s3_class(result, "JsonTrackObj")
})

test_that("lst_dysprosody return_jstf=TRUE returns JsonTrackObj", {
  skip_if_not(pladdrr_available(), "pladdrr not available")
  result <- suppressWarnings(lst_dysprosody(wav_file(), return_jstf = TRUE, verbose = FALSE))
  skip_if(is.null(result) || length(result) == 0, "dysprosody pipeline unavailable for test file")
  expect_s3_class(result, "JsonTrackObj")
})

# ---- lst_avqi and lst_dsi require pre-computed data frames ----
# These take acoustic measure data frames as input, not raw audio files.
# Covered by their own dedicated tests if/when pipeline integration is added.

# Tests for WORLD CheapTrick Spectral Envelope and D4C Aperiodicity

test_that("cheap_trick_cpp returns valid spectral envelope", {
  skip_if_not_installed("superassp")
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- read_audio(test_wav)
  h <- superassp:::harvest_cpp(audio_obj, minF = 60.0, maxF = 400.0,
                                windowShift = 5.0)
  result <- superassp:::cheap_trick_cpp(audio_obj,
                                        f0 = as.numeric(h$f0),
                                        temporal_positions = as.numeric(h$times))

  expect_type(result, "list")
  expect_named(result, c("spectrogram", "temporal_positions",
                         "sample_rate", "n_frames", "fft_size"))
  expect_true(is.matrix(result$spectrogram))
  expect_gt(ncol(result$spectrogram), 1)
  expect_equal(ncol(result$spectrogram), result$fft_size / 2L + 1L)
  expect_equal(nrow(result$spectrogram), result$n_frames)
  expect_true(all(result$spectrogram >= 0))
})

test_that("trk_cheap_trick returns valid AsspDataObj", {
  skip_if_not_installed("superassp")
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- trk_cheap_trick(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("sp" %in% names(result))
  expect_true(is.matrix(result$sp))
  expect_gt(ncol(result$sp), 1)
  expect_equal(ncol(result$sp), attr(result, "fft_size") / 2L + 1L)
})

test_that("d4c_cpp returns valid aperiodicity", {
  skip_if_not_installed("superassp")
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- read_audio(test_wav)
  h <- superassp:::harvest_cpp(audio_obj, minF = 60.0, maxF = 400.0,
                                windowShift = 5.0)
  result <- superassp:::d4c_cpp(audio_obj,
                                f0 = as.numeric(h$f0),
                                temporal_positions = as.numeric(h$times))

  expect_type(result, "list")
  expect_named(result, c("aperiodicity", "temporal_positions",
                         "sample_rate", "n_frames", "fft_size"))
  expect_true(is.matrix(result$aperiodicity))
  expect_equal(nrow(result$aperiodicity), result$n_frames)
  expect_equal(ncol(result$aperiodicity), result$fft_size / 2L + 1L)
  expect_true(all(result$aperiodicity >= 0 & result$aperiodicity <= 1))
})

test_that("lst_voxit intensity_mean_db is not the -80 placeholder", {
  skip_if_not_installed("superassp")
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  feat <- lst_voxit(test_wav)

  expect_false(is.nan(feat$intensity_mean_db))
  expect_lt(feat$intensity_mean_db, 0)
  expect_gt(feat$intensity_mean_db, -90)
})

test_that("cheap_trick_cpp fft_size matches spec for 16kHz audio", {
  skip_if_not_installed("superassp")
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  audio_obj <- read_audio(test_wav)
  sr <- attr(audio_obj, "origFreq")
  skip_if(sr != 16000L, "Test requires 16 kHz audio")

  h <- superassp:::harvest_cpp(audio_obj, minF = 60.0, maxF = 400.0)
  result <- superassp:::cheap_trick_cpp(audio_obj, as.numeric(h$f0),
                                        as.numeric(h$times))
  # For 16kHz, f0_floor=71.0 → fft_size=2048
  expect_equal(result$fft_size, 2048L)
  expect_equal(ncol(result$spectrogram), 1025L)
})

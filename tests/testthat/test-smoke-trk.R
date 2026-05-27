# Smoke tests for trk_* functions that had zero coverage.
# Each test checks: returns AsspDataObj, has expected track, has > 0 rows.

wav_file <- function() {
  f <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  if (f == "") skip("Test wav not found")
  f
}

# ---- C-ASSP time-domain tracks ----

test_that("trk_acf returns AsspDataObj with ACF track", {
  result <- trk_acf(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_true("ACF" %in% names(result))
  expect_gt(nrow(result$ACF), 0)
})

test_that("trk_rms returns AsspDataObj with rms track", {
  result <- trk_rms(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_gt(length(names(result)), 0)
  expect_gt(nrow(result[[1]]), 0)
})

test_that("trk_zcr returns AsspDataObj with zcr track", {
  result <- trk_zcr(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_gt(length(names(result)), 0)
  expect_gt(nrow(result[[1]]), 0)
})

test_that("trk_lpc returns AsspDataObj with lpc track", {
  result <- trk_lpc(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_true("LPC" %in% names(result))
  expect_gt(nrow(result$LPC), 0)
})

test_that("trk_rfc returns AsspDataObj with rfc track", {
  result <- trk_rfc(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_true("RFC" %in% names(result))
  expect_gt(nrow(result$RFC), 0)
})

test_that("trk_arf returns AsspDataObj with arf track", {
  result <- trk_arf(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_true("ARF" %in% names(result) || "arf" %in% names(result))
  expect_gt(nrow(result[[1]]), 0)
})

test_that("trk_lar returns AsspDataObj with lar track", {
  result <- trk_lar(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_true(length(names(result)) > 0)
  expect_gt(nrow(result[[1]]), 0)
})

# ---- Spectrum tracks ----

test_that("trk_dft_spectrum returns AsspDataObj with DFT track", {
  result <- trk_dft_spectrum(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_true(any(grepl("DFT", names(result), ignore.case = TRUE)))
  expect_gt(nrow(result[[1]]), 0)
})

test_that("trk_cepstrum returns AsspDataObj with cepstrum track", {
  result <- trk_cepstrum(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_gt(nrow(result[[1]]), 0)
})

# ---- SPTK pitch trackers ----

test_that("trk_pitch_dio returns AsspDataObj with f0 track", {
  result <- trk_pitch_dio(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_true("f0" %in% names(result))
  expect_gt(nrow(result$f0), 0)
})

test_that("trk_pitch_harvest returns AsspDataObj with f0 track", {
  result <- trk_pitch_harvest(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_true("f0" %in% names(result))
  expect_gt(nrow(result$f0), 0)
})

test_that("trk_pitch_yin returns AsspDataObj with f0 track", {
  result <- trk_pitch_yin(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_true("F0" %in% names(result))
  expect_gt(nrow(result$F0), 0)
})

# ---- MFCC ----

test_that("trk_mfcc returns AsspDataObj with MFCC track", {
  result <- trk_mfcc(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_gt(nrow(result[[1]]), 0)
})

# ---- WORLD vocoder ----

test_that("trk_d4c returns AsspDataObj with aperiodicity track", {
  result <- trk_d4c(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_gt(nrow(result[[1]]), 0)
})

test_that("trk_cheap_trick returns AsspDataObj with spectrogram track", {
  result <- trk_cheap_trick(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_gt(nrow(result[[1]]), 0)
})

# ---- Formant trackers ----

test_that("trk_formant_burg returns AsspDataObj with formant tracks", {
  result <- trk_formant_burg(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_gt(nrow(result[[1]]), 0)
})

test_that("trk_formant_forest returns AsspDataObj with formant tracks", {
  result <- trk_formant_forest(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_gt(nrow(result[[1]]), 0)
})

# ---- VUV and intensity ----

test_that("trk_vuv returns a non-null result", {
  result <- trk_vuv(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_true(!is.null(result))
})

test_that("trk_intensity returns AsspDataObj with intensity track", {
  result <- trk_intensity(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_gt(nrow(result[[1]]), 0)
})

test_that("trk_cpps returns AsspDataObj or numeric result", {
  result <- trk_cpps(wav_file(), toFile = FALSE, verbose = FALSE)
  expect_true(inherits(result, "AsspDataObj") || is.numeric(result) || is.list(result))
})

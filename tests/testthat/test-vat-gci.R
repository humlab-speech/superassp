test_that("trk_gci_vat returns GCIs at the audio sample rate", {
  skip_if_not_installed("voiceanalysis")
  wav <- system.file("samples/sustained/a1.wav", package = "superassp")
  skip_if(wav == "", "a1.wav not bundled")
  obj <- trk_gci_vat(wav, toFile = FALSE)
  expect_s3_class(obj, "AsspDataObj")
  expect_named(obj, c("gci_sample", "residual"))
  expect_gt(nrow(obj$gci_sample), 10)
  fs <- attr(obj, "sampleRate")
  # GCIs must be inside the residual length
  expect_true(all(obj$gci_sample[, 1] >= 1 & obj$gci_sample[, 1] <= nrow(obj$residual)))
  # Median inter-GCI interval in [fs/500, fs/50] (50–500 Hz)
  igi <- diff(obj$gci_sample[, 1])
  expect_true(median(igi) > fs / 500 && median(igi) < fs / 50)
})

test_that("trk_gci_vat var_f0 = TRUE also runs", {
  skip_if_not_installed("voiceanalysis")
  wav <- system.file("samples/sustained/a1.wav", package = "superassp")
  skip_if(wav == "", "a1.wav not bundled")
  obj <- trk_gci_vat(wav, var_f0 = TRUE, toFile = FALSE)
  expect_s3_class(obj, "AsspDataObj")
  expect_gt(nrow(obj$gci_sample), 10)
})

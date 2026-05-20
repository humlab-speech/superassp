test_that("trk_pitch_vat returns AsspDataObj with expected tracks", {
  skip_if_not_installed("voiceanalysis")
  wav <- system.file("samples/sustained/a1.wav", package = "superassp")
  skip_if(wav == "", "a1.wav not bundled")
  obj <- trk_pitch_vat(wav, toFile = FALSE)
  expect_s3_class(obj, "AsspDataObj")
  expect_named(obj, c("f0", "vad", "srh_val"))
  expect_equal(attr(obj, "sampleRate"), 100)
  voiced <- obj$f0[obj$vad == 1]
  if (length(voiced) > 5) {
    expect_true(all(voiced > 50 & voiced < 500))
  }
})

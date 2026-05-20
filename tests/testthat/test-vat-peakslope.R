test_that("trk_peakslope_vat returns a single REAL32 track at 100 Hz", {
  skip_if_not_installed("voiceanalysis")
  wav <- system.file("samples/sustained/a1.wav", package = "superassp")
  skip_if(wav == "", "a1.wav not bundled")
  obj <- trk_peakslope_vat(wav, toFile = FALSE)
  expect_s3_class(obj, "AsspDataObj")
  expect_named(obj, "peak_slope")
  expect_equal(attr(obj, "sampleRate"), 100)
  expect_true(all(is.finite(obj$peak_slope)))
  # peak slope typically negative for voiced speech, in [-1, 0.5]
  expect_lt(median(obj$peak_slope), 0.5)
  expect_gt(median(obj$peak_slope), -2.0)
})

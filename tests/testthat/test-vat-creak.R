test_that("trk_creak_vat returns bounded posterior + binary decision", {
  skip_if_not_installed("voiceanalysis")
  wav <- system.file("samples/sustained/a1.wav", package = "superassp")
  skip_if(wav == "", "a1.wav not bundled")
  obj <- trk_creak_vat(wav, toFile = FALSE)
  expect_s3_class(obj, "AsspDataObj")
  expect_named(obj, c("creak_pp", "creak_bin"))
  expect_equal(attr(obj, "sampleRate"), 100)
  expect_true(all(obj$creak_pp >= 0 & obj$creak_pp <= 1))
  expect_true(all(obj$creak_bin %in% c(0, 1)))
  # decision must respect threshold
  expect_true(all((obj$creak_pp > 0.3) == (obj$creak_bin == 1)))
})

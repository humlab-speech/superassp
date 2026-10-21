test_that("trk_affilter defaults to wrassp's 4 kHz high-pass", {
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wav == "", "Test file not found")
  def <- trk_affilter(wav, verbose = FALSE)
  expect_s3_class(def, "AsspDataObj")
  expect_identical(def, trk_affilter(wav, highPass = 4000, lowPass = 0, verbose = FALSE))
  # NULL and 0 both disable a band edge
  expect_identical(trk_affilter(wav, highPass = NULL, lowPass = 1000, verbose = FALSE),
                   trk_affilter(wav, highPass = 0, lowPass = 1000, verbose = FALSE))
  expect_error(trk_affilter(wav, highPass = 0, verbose = FALSE), "must be positive")
})

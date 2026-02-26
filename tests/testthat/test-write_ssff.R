test_that("write_ssff roundtrips an AsspDataObj", {
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wav == "", "test wav not found")

  obj  <- read_ssff(wav)
  tmp  <- tempfile(fileext = ".wav")
  on.exit(unlink(tmp))

  write_ssff(obj, tmp)
  expect_true(file.exists(tmp))

  obj2 <- read_ssff(tmp)
  expect_s3_class(obj2, "AsspDataObj")
  expect_equal(attr(obj, "sampleRate"), attr(obj2, "sampleRate"))
})

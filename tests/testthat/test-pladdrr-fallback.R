test_that("av_load_for_pladdrr falls back to av for non-native formats", {
  skip_if(!superassp:::pladdrr_available(), "pladdrr not installed")
  skip_if_not_installed("av")

  # OGG is not natively supported by pladdrr but av can handle it
  ogg <- system.file("samples", "sustained", "a8.ogg", package = "superassp")
  skip_if(ogg == "" || !file.exists(ogg), "test ogg not found")

  result <- superassp:::av_load_for_pladdrr(ogg)
  expect_true(!is.null(result))
})

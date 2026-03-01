# Tests for Deprecation Warnings
# Tests that deprecated functions properly warn users about replacements

test_that("trk_yaapt triggers deprecation warning", {
  skip_if_not_installed("superassp")
  skip_if_not(reticulate::py_module_available("pyyaapt"), "pyyaapt not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  expect_warning(
    superassp::trk_yaapt(test_wav, toFile = FALSE),
    "deprecated"
  )
})

test_that("trk_sacc triggers deprecation warning", {
  skip_if_not_installed("superassp")
  skip_if_not(reticulate::py_module_available("sacc"), "sacc not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  expect_warning(
    superassp::trk_sacc(test_wav, toFile = FALSE),
    "deprecated"
  )
})

test_that("trk_straight_f0 triggers deprecation warning", {
  skip_if_not_installed("superassp")
  skip_if_not(superassp::straight_available(), "STRAIGHT not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  expect_warning(
    superassp::trk_straight_f0(test_wav, toFile = FALSE),
    "deprecated"
  )
})

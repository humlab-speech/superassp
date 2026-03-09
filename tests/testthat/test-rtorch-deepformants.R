## Tests for R torch DeepFormants implementation
## (feature extraction + model inference without Python)

test_that("deepformants_available() returns a single logical", {
  result <- deepformants_available()
  expect_true(is.logical(result))
  expect_length(result, 1L)
})

test_that("df_specps_features() returns 50 finite numeric values", {
  set.seed(42L)
  frame <- as.integer(rnorm(480L) * 1000L)
  result <- df_specps_features(frame)
  expect_length(result, 50L)
  expect_true(is.numeric(result))
  expect_false(any(is.nan(result)))
  expect_false(any(is.infinite(result)))
})

test_that("df_arspec_features() returns 30 finite values for each LPC order 8-17", {
  set.seed(42L)
  frame <- as.double(rnorm(480L) * 1000L)
  for (order in 8:17) {
    result <- df_arspec_features(frame, order)
    expect_length(result, 30L)
    expect_true(is.numeric(result))
    expect_false(any(is.nan(result)), info = paste("NaN for order", order))
    expect_false(any(is.infinite(result)), info = paste("Inf for order", order))
  }
})

test_that("extract_deepformants_features() tracking mode: (n_frames, 350) matrix", {
  skip_if_not(requireNamespace("av", quietly = TRUE))
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wav == "", "test audio not found")

  mat <- extract_deepformants_features(wav)
  expect_true(is.matrix(mat))
  expect_equal(ncol(mat), 350L)
  expect_gt(nrow(mat), 0L)
  expect_false(any(is.nan(mat)))
})

test_that("extract_deepformants_features() estimation mode: (1, 350) matrix", {
  skip_if_not(requireNamespace("av", quietly = TRUE))
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wav == "", "test audio not found")

  mat <- extract_deepformants_features(wav, begin = 0.0, end = 0.1)
  expect_true(is.matrix(mat))
  expect_equal(ncol(mat), 350L)
  expect_equal(nrow(mat), 1L)
  expect_false(any(is.nan(mat)))
})

## Tasks 5-6 tests (model loading + inference) — run once other agent's work lands
## test_that("load_deepformants_estimator() loads model", { ... })
## test_that("load_deepformants_tracker() loads model", { ... })
## test_that("trk_deepformants() returns AsspDataObj with R torch", { ... })
## test_that("lst_deepformants() returns named list with R torch", { ... })

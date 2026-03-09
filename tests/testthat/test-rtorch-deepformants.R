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

## Task 5: model loading + inference

test_that("load_deepformants_estimator() loads MLP and forward pass is (1,4)", {
  skip_if_not(requireNamespace("torch", quietly = TRUE))
  model <- load_deepformants_estimator()
  expect_true(inherits(model, "nn_module"))
  x <- torch::torch_randn(c(1L, 350L))
  torch::with_no_grad(out <- model(x))
  expect_equal(out$shape, c(1L, 4L))
})

test_that("load_deepformants_tracker() loads LSTM and forward pass is (1,n,4)", {
  skip_if_not(requireNamespace("torch", quietly = TRUE))
  model <- load_deepformants_tracker()
  expect_true(inherits(model, "nn_module"))
  x <- torch::torch_randn(c(1L, 10L, 350L))
  torch::with_no_grad(out <- model(x))
  expect_equal(out$shape, c(1L, 10L, 4L))
})

test_that("run_deepformants_tracker() returns n_frames x 4 Hz matrix", {
  skip_if_not(requireNamespace("torch", quietly = TRUE))
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wav == "", "test audio not found")
  feat_mat <- extract_deepformants_features(wav)
  result   <- run_deepformants_tracker(feat_mat)
  expect_true(is.matrix(result))
  expect_equal(ncol(result), 4L)
  expect_equal(nrow(result), nrow(feat_mat))
  expect_true(all(result > 0, na.rm = TRUE))
})

test_that("run_deepformants_estimator() returns length-4 positive Hz vector", {
  skip_if_not(requireNamespace("torch", quietly = TRUE))
  wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(wav == "", "test audio not found")
  feat_mat <- extract_deepformants_features(wav, begin = 0.0, end = 0.1)
  result   <- run_deepformants_estimator(feat_mat)
  expect_length(result, 4L)
  expect_true(all(result > 0))
})

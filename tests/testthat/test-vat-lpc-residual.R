test_that("vat_lpc_residual_cpp produces a normalized, finite residual of matching length", {
  skip_if_not_installed("superassp")
  set.seed(11)
  wave <- rnorm(4000)
  res <- superassp:::vat_lpc_residual_cpp(wave, 400L, 100L, 12L)

  expect_length(res, length(wave))
  expect_true(all(is.finite(res)))
  expect_lte(max(abs(res)), 1 + 1e-9)
})

test_that("vat_lpc_residual_cpp is deterministic across repeated calls (golden-master baseline)", {
  skip_if_not_installed("superassp")
  set.seed(11)
  wave <- rnorm(4000)
  res1 <- superassp:::vat_lpc_residual_cpp(wave, 400L, 100L, 12L)
  res2 <- superassp:::vat_lpc_residual_cpp(wave, 400L, 100L, 12L)
  expect_identical(res1, res2)
})

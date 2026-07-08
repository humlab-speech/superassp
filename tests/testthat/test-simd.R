# Faithfulness tests for the double-precision SIMD primitives in
# src/simd_utils.hpp. The vectorized path must match an independent scalar
# reference (base R) within double-precision tolerance. Bindings are internal.

test_that("simd_dot matches crossprod", {
  set.seed(1)
  for (n in c(1L, 3L, 4L, 7L, 8L, 100L, 1023L)) {
    a <- runif(n, -1, 1)
    b <- runif(n, -1, 1)
    expect_equal(
      superassp:::simd_dot_cpp(a, b),
      as.numeric(crossprod(a, b)),
      info = paste("n =", n)
    )
  }
})

test_that("simd_energy matches sum of squares", {
  set.seed(2)
  for (n in c(1L, 5L, 8L, 256L, 999L)) {
    x <- rnorm(n)
    expect_equal(
      superassp:::simd_energy_cpp(x),
      sum(x^2),
      info = paste("n =", n)
    )
  }
})

test_that("simd_fir matches causal filter(b, 1, x)", {
  set.seed(3)
  x <- rnorm(64)
  b <- c(0.2, 0.5, 0.3, -0.1)
  # stats::filter(sides = 1) yields NA during ramp-up; compare only full-support
  ref <- stats::filter(x, b, method = "convolution", sides = 1)
  got <- superassp:::simd_fir_cpp(x, b)
  full <- length(b):length(x)
  expect_equal(got[full], as.numeric(ref[full]))

  # Single-tap FIR is a pure scale — exact.
  expect_equal(superassp:::simd_fir_cpp(x, 2), 2 * x)
})

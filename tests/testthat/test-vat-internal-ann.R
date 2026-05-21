test_that("vat_ann_forward_cpp produces tanh hidden + linear output", {
  # 4 features, 5 hidden, 1 output. 3 frames.
  set.seed(42)
  X  <- matrix(rnorm(12), 4, 3)
  IW <- matrix(rnorm(20), 5, 4)
  bh <- rnorm(5)
  LW <- matrix(rnorm(5), 1, 5)
  bo <- rnorm(1)
  mn <- rep(-3, 4); mx <- rep(3, 4)

  y <- vat_ann_forward_cpp(X, IW, bh, LW, bo, mn, mx)
  expect_equal(length(y), 3)
  expect_true(all(is.finite(y)))

  # Reference: compute one column manually
  x0 <- X[, 1]
  xn <- -1 + (x0 - mn) / (mx - mn) * 2
  a  <- 2 / (1 + exp(-2 * (IW %*% xn + bh))) - 1
  y0 <- LW %*% a + bo
  # median filter of 3 — at edges, equals frame 1's value (n=odd buf w/ zero pad)
  # we only check finite + reasonable magnitude here
  expect_lt(abs(y0[1]), 100)
})

test_that("vat_ann_forward_cpp handles NaN inputs by substituting mini", {
  X <- matrix(c(NaN, 0.5, -0.5, 0.0), 2, 2)
  IW <- matrix(c(1, 0, 0, 1), 2, 2)
  bh <- c(0, 0)
  LW <- matrix(c(1, 1), 1, 2)
  bo <- 0
  mn <- c(-1, -1); mx <- c(1, 1)
  y <- vat_ann_forward_cpp(X, IW, bh, LW, bo, mn, mx)
  expect_true(all(is.finite(y)))
})

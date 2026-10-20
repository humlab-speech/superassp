# The Windows batch path reuses one PSOCK cluster across calls. PSOCK also
# works on Unix, so the cache logic is exercised everywhere.

test_that("PSOCK cluster is reused, resized and replaced when dead", {
  skip_on_cran()
  withr::defer(superassp:::.stop_superassp_cluster())

  cl1 <- superassp:::.superassp_cluster(2L)
  expect_identical(superassp:::.superassp_cluster(2L), cl1)
  expect_equal(unlist(parallel::clusterEvalQ(cl1, 1 + 1)), c(2, 2))

  # A dead cluster is detected and replaced by a working one
  parallel::stopCluster(cl1)
  cl2 <- superassp:::.superassp_cluster(2L)
  expect_false(identical(cl2, cl1))
  expect_equal(unlist(parallel::clusterEvalQ(cl2, 1 + 1)), c(2, 2))

  # A different size starts a new cluster
  cl3 <- superassp:::.superassp_cluster(1L)
  expect_length(cl3, 1L)
  expect_equal(unlist(parallel::clusterEvalQ(cl3, 1 + 1)), 2)

  superassp:::.stop_superassp_cluster()
  expect_null(superassp:::.cluster_cache$cl)
})

test_that("run_parallel_files sequential path preserves order and values", {
  process_fn <- function(i) i * 10L
  results <- superassp:::run_parallel_files(
    n_files = 5, process_single_file = process_fn,
    parallel = FALSE, verbose = FALSE
  )
  expect_equal(unlist(results), c(10L, 20L, 30L, 40L, 50L))
})

test_that("run_parallel_files parallel path matches sequential path", {
  skip_on_cran()
  skip_on_os("windows")

  process_fn <- function(i) i^2
  seq_results <- superassp:::run_parallel_files(
    n_files = 6, process_single_file = process_fn,
    parallel = FALSE, verbose = FALSE
  )
  par_results <- superassp:::run_parallel_files(
    n_files = 6, process_single_file = process_fn,
    parallel = TRUE, n_cores = 2, verbose = FALSE
  )
  expect_equal(unlist(par_results), unlist(seq_results))
})

test_that("run_parallel_files handles n_files = 1 without a progress bar", {
  process_fn <- function(i) "ok"
  results <- superassp:::run_parallel_files(
    n_files = 1, process_single_file = process_fn, verbose = TRUE
  )
  expect_equal(results[[1]], "ok")
})

test_that("run_parallel_files propagates per-file errors caught inside the closure", {
  process_fn <- function(i) {
    tryCatch({
      if (i == 2) stop("boom")
      i
    }, error = function(e) NA_integer_)
  }
  results <- superassp:::run_parallel_files(
    n_files = 3, process_single_file = process_fn,
    parallel = FALSE, verbose = FALSE
  )
  expect_equal(unlist(results), c(1L, NA_integer_, 3L))
})

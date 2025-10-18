context("Parallel Processing and Thread Safety")

library(testthat)
library(superassp)

# Get test files
get_test_files <- function(n = 3) {
  files <- list.files(
    system.file('samples', 'sustained', package='superassp'),
    pattern = "\\.wav$",
    full.names = TRUE
  )
  if (length(files) == 0) {
    skip("No test files available")
  }
  rep(files[1:min(length(files), n)], ceiling(10/length(files)))[1:10]
}

test_that("Parallel processing produces identical results to sequential", {
  skip_on_cran()

  test_files <- get_test_files()

  # Process sequentially (single file at a time ensures sequential)
  results_seq <- lapply(test_files, function(f) {
    trk_rmsana(f, toFile = FALSE, verbose = FALSE)
  })

  # Process in parallel (batch triggers parallel)
  results_par <- trk_rmsana(test_files, toFile = FALSE, verbose = FALSE)

  # Check we got results
  expect_true(is.list(results_par))
  expect_equal(length(results_par), length(test_files))

  # Compare results - should be identical
  for (i in seq_along(test_files)) {
    # Check both are AsspDataObj
    expect_true(is.AsspDataObj(results_seq[[i]]))
    expect_true(is.AsspDataObj(results_par[[i]]))

    # Check same number of records
    expect_equal(
      numRecs.AsspDataObj(results_seq[[i]]),
      numRecs.AsspDataObj(results_par[[i]])
    )

    # Check same sample rate
    expect_equal(
      rate.AsspDataObj(results_seq[[i]]),
      rate.AsspDataObj(results_par[[i]])
    )

    # Check data values are identical (within floating point tolerance)
    expect_equal(
      results_seq[[i]][["RMS[dB]"]],
      results_par[[i]][["RMS[dB]"]],
      tolerance = 1e-10
    )
  }
})

test_that("All DSP functions work with parallel processing", {
  skip_on_cran()

  test_files <- get_test_files(2)

  # Test a representative set of DSP functions
  # Note: Using simple functions that don't require track names
  functions_to_test <- list(
    rmsana = trk_rmsana,
    zcrana = trk_zcrana
  )

  for (fname in names(functions_to_test)) {
    f <- functions_to_test[[fname]]

    # Should not error
    expect_silent({
      results <- f(test_files, toFile = FALSE, verbose = FALSE)
    })

    # Should return list of results
    expect_true(is.list(results))
    expect_equal(length(results), length(test_files))

    # Each should be AsspDataObj
    for (r in results) {
      expect_true(is.AsspDataObj(r),
                  info = paste(fname, "should return AsspDataObj"))
    }
  }
})

test_that("Parallel processing handles errors gracefully", {
  skip_on_cran()

  # Mix of valid and invalid files
  test_files <- c(
    get_test_files(1),
    "/nonexistent/file.wav"
  )

  # Should error appropriately, not crash
  expect_error({
    trk_rmsana(test_files, toFile = FALSE, verbose = FALSE)
  })
})

test_that("Single file processing remains sequential", {
  skip_on_cran()

  test_file <- get_test_files(1)[1]

  # Process single file - should work fine
  result <- trk_rmsana(test_file, toFile = FALSE, verbose = FALSE)

  expect_true(is.AsspDataObj(result))
  expect_gt(numRecs.AsspDataObj(result), 0)
})

test_that("Thread safety: No race conditions in parallel processing", {
  skip_on_cran()
  skip_on_os("windows")  # Fork-based test

  test_files <- get_test_files(5)

  # Run same batch multiple times - results should be identical
  results1 <- trk_rmsana(test_files, toFile = FALSE, verbose = FALSE)
  results2 <- trk_rmsana(test_files, toFile = FALSE, verbose = FALSE)
  results3 <- trk_rmsana(test_files, toFile = FALSE, verbose = FALSE)

  # All runs should produce identical results
  for (i in seq_along(test_files)) {
    expect_equal(
      results1[[i]][["RMS[dB]"]],
      results2[[i]][["RMS[dB]"]],
      tolerance = 1e-10
    )
    expect_equal(
      results2[[i]][["RMS[dB]"]],
      results3[[i]][["RMS[dB]"]],
      tolerance = 1e-10
    )
  }
})

test_that("Memory is properly managed in parallel processing", {
  skip_on_cran()

  test_files <- get_test_files(5)

  # Get baseline memory
  gc()
  mem_before <- as.numeric(gc()[2, 2])  # max used

  # Process batch
  results <- trk_rmsana(test_files, toFile = FALSE, verbose = FALSE)

  # Clear results and force garbage collection
  rm(results)
  gc()
  gc()  # Call twice to ensure full collection
  mem_after <- as.numeric(gc()[2, 2])

  # Memory growth should be reasonable (less than 100MB for small test files)
  mem_growth_mb <- mem_after - mem_before
  expect_lt(mem_growth_mb, 100,
            label = paste("Memory growth:", mem_growth_mb, "MB"))
})

test_that("Parallel processing respects file order", {
  skip_on_cran()

  test_files <- get_test_files(3)

  # Add identifiable time windows
  results <- lapply(seq_along(test_files), function(i) {
    trk_rmsana(test_files[i], toFile = FALSE, verbose = FALSE)
  })

  results_par <- trk_rmsana(test_files, toFile = FALSE, verbose = FALSE)

  # Results should be in same order
  for (i in seq_along(test_files)) {
    expect_equal(
      results[[i]][["RMS[dB]"]],
      results_par[[i]][["RMS[dB]"]],
      tolerance = 1e-10
    )
  }
})

test_that("Process function is truly thread-safe (stress test)", {
  skip_on_cran()
  skip_if(parallel::detectCores() < 4, "Need at least 4 cores for stress test")

  test_files <- get_test_files(3)

  # Create larger batch by repeating files
  large_batch <- rep(test_files, 10)  # 30 files

  # Process large batch - should not crash or corrupt data
  expect_silent({
    results <- trk_rmsana(large_batch, toFile = FALSE, verbose = FALSE)
  })

  # All results should be valid
  expect_equal(length(results), length(large_batch))
  for (r in results) {
    expect_true(is.AsspDataObj(r))
    expect_gt(numRecs.AsspDataObj(r), 0)
  }

  # Verify results are consistent (same file should give same result)
  # Check first occurrence vs last occurrence of first file
  first_idx <- 1
  last_idx <- length(large_batch) - length(test_files) + 1

  expect_equal(
    results[[first_idx]][["RMS[dB]"]],
    results[[last_idx]][["RMS[dB]"]],
    tolerance = 1e-10
  )
})

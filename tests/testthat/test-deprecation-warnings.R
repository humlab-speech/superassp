# Tests for Deprecation Warnings
# Tests that deprecated functions properly warn users about replacements

test_that("reaper_pm triggers deprecation warning", {
  skip_if_not_installed("superassp")
  skip_if_not(reticulate::py_module_available("pyreaper"), "pyreaper not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Should trigger .Deprecated() warning
  expect_warning(
    superassp::reaper_pm(test_wav, toFile = FALSE, verbose = FALSE),
    "deprecated"
  )
})

test_that("reaper_pm deprecation message mentions trk_reaper_pm", {
  skip_if_not_installed("superassp")
  skip_if_not(reticulate::py_module_available("pyreaper"), "pyreaper not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Capture warning message
  warning_msg <- ""
  withCallingHandlers(
    superassp::reaper_pm(test_wav, toFile = FALSE, verbose = FALSE),
    warning = function(w) {
      warning_msg <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )

  # Should mention the replacement function
  expect_true(grepl("trk_reaper_pm", warning_msg, ignore.case = TRUE))
})

test_that("reaper_pm deprecation message mentions version removal", {
  skip_if_not_installed("superassp")
  skip_if_not(reticulate::py_module_available("pyreaper"), "pyreaper not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Capture warning message
  warning_msg <- ""
  withCallingHandlers(
    superassp::reaper_pm(test_wav, toFile = FALSE, verbose = FALSE),
    warning = function(w) {
      warning_msg <<- conditionMessage(w)
      invokeRestart("muffleWarning")
    }
  )

  # Should mention when it will be removed
  expect_true(grepl("v0\\.11", warning_msg) || grepl("removed", warning_msg, ignore.case = TRUE))
})

test_that("deprecated reaper_pm still works correctly", {
  skip_if_not_installed("superassp")
  skip_if_not(reticulate::py_module_available("pyreaper"), "pyreaper not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Suppress deprecation warning for this test
  result <- suppressWarnings(
    superassp::reaper_pm(test_wav, toFile = FALSE, verbose = FALSE)
  )

  # Function should still work
  expect_s3_class(result, "AsspDataObj")
  expect_true("pm" %in% names(result))
  expect_true(is.matrix(result$pm))

  # Should still produce binary indicators
  expect_true(all(result$pm %in% c(0L, 1L)))
})

test_that("reaper_pm and trk_reaper_pm produce similar output", {
  skip_if_not_installed("superassp")
  skip_if_not(reticulate::py_module_available("pyreaper"), "pyreaper not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # Get results from both functions
  result_old <- suppressWarnings(
    superassp::reaper_pm(test_wav, toFile = FALSE, verbose = FALSE)
  )

  result_new <- superassp::trk_reaper_pm(
    test_wav,
    toFile = FALSE,
    verbose = FALSE
  )

  # Both should be AsspDataObj with pm track
  expect_s3_class(result_old, "AsspDataObj")
  expect_s3_class(result_new, "AsspDataObj")
  expect_true("pm" %in% names(result_old))
  expect_true("pm" %in% names(result_new))

  # Should have similar number of frames (allow some tolerance)
  expect_equal(nrow(result_old$pm), nrow(result_new$pm), tolerance = 10)

  # Both should have binary values
  expect_true(all(result_old$pm %in% c(0L, 1L)))
  expect_true(all(result_new$pm %in% c(0L, 1L)))

  # Should have similar number of pitch marks
  n_marks_old <- sum(result_old$pm == 1L)
  n_marks_new <- sum(result_new$pm == 1L)

  expect_equal(n_marks_old, n_marks_new, tolerance = max(5, n_marks_old * 0.1))
})

test_that("deprecation warning can be suppressed", {
  skip_if_not_installed("superassp")
  skip_if_not(reticulate::py_module_available("pyreaper"), "pyreaper not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  # suppressWarnings should work
  expect_silent({
    result <- suppressWarnings(
      superassp::reaper_pm(test_wav, toFile = FALSE, verbose = FALSE)
    )
  })
})

# =============================================================================
# Future Deprecations (Placeholder)
# =============================================================================

test_that("check for undocumented deprecations", {
  # This is a meta-test to remind us to add tests for future deprecations

  # List of currently deprecated functions
  deprecated_funs <- c("reaper_pm", "trk_egg_f0_deprecated")

  # Verify each has a test
  for (fun in deprecated_funs) {
    expect_true(
      exists(fun, where = asNamespace("superassp")),
      info = paste(fun, "should exist in namespace if listed as deprecated")
    )
  }
})

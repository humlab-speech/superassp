# Helper function to normalize track names for comparison
# Removes bracket notation like [dB], [Hz], etc. and converts to lowercase
normalize_track_name <- function(name) {
  # Remove anything in brackets and convert to lowercase
  tolower(gsub("\\[.*?\\]", "", name))
}

# NOTE: The "Parselmouth optimized functions" test block was removed — the
# praat_*_opt functions were deleted in the Python/reticulate purge.

test_that("SuperASP and wrassp functions produce compatible outputs", {
  skip_on_cran()
  # superassp is self-contained; the wrassp equivalence comparison is retired
  # (loading wrassp masks superassp's S3 methods for the rest of the suite).
  skip("wrassp comparison retired: superassp is self-contained")

  # Get a test file
  test_file <- list.files(
    path = testthat::test_path("..", "signalfiles"),
    pattern = "\\.wav$",
    recursive = TRUE,
    full.names = TRUE
  )[1]

  skip_if(is.na(test_file) || !file.exists(test_file),
          "Test file not found")

  # Test rmsana
  result_superassp <- superassp::trk_rms(test_file, toFile = FALSE)
  result_wrassp <- wrassp::trk_rms(test_file, toFile = FALSE)

  expect_s3_class(result_superassp, "AsspDataObj")
  expect_s3_class(result_wrassp, "AsspDataObj")
  # Compare track names (normalized to remove bracket notation and case)
  expect_equal(
    normalize_track_name(names(result_superassp)),
    normalize_track_name(names(result_wrassp))
  )
  expect_equal(dim(result_superassp[[1]]), dim(result_wrassp[[1]]))

  # Test acfana
  result_superassp <- superassp::trk_acf(test_file, toFile = FALSE)
  result_wrassp <- wrassp::trk_acf(test_file, toFile = FALSE)

  # Compare track names (normalized to remove bracket notation and case)
  expect_equal(
    normalize_track_name(names(result_superassp)),
    normalize_track_name(names(result_wrassp))
  )
  expect_equal(nrow(result_superassp[[1]]), nrow(result_wrassp[[1]]))
})

# Helper function to normalize track names for comparison
# Removes bracket notation like [dB], [Hz], etc. and converts to lowercase
normalize_track_name <- function(name) {
  # Remove anything in brackets and convert to lowercase
  tolower(gsub("\\[.*?\\]", "", name))
}

test_that("Parselmouth optimized functions return valid AsspDataObj", {
  skip_if_not_installed("reticulate")

  # Setup Python
  Sys.unsetenv("RETICULATE_PYTHON")
  reticulate::use_python("/usr/bin/python3", required = FALSE)

  # Check if parselmouth is available
  skip_if_not(reticulate::py_module_available("parselmouth"),
              "Parselmouth not available")

  # Get a test file
  test_file <- list.files(
    path = file.path("..", "signalfiles"),
    pattern = "sv1\\.wav$",
    recursive = TRUE,
    full.names = TRUE
  )[1]

  skip_if(is.na(test_file) || !file.exists(test_file),
          "Test file not found")

  # Test praat_formant_burg_opt
  result <- praat_formant_burg_opt(test_file, toFile = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_true(length(names(result)) > 0)
  expect_true(nrow(result[[1]]) > 0)

  # Test praat_pitch_opt
  result <- praat_pitch_opt(test_file, toFile = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_true("cc" %in% names(result) || "ac" %in% names(result))

  # Test praat_intensity_opt
  result <- praat_intensity_opt(test_file, toFile = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_true("intensity" %in% names(result))
  expect_true(nrow(result$intensity) > 0)
})

test_that("SuperASP and wrassp functions produce compatible outputs", {
  skip_on_cran()
  skip_if_not_installed("wrassp")

  # Get a test file
  test_file <- list.files(
    path = file.path("..", "signalfiles"),
    pattern = "\\.wav$",
    recursive = TRUE,
    full.names = TRUE
  )[1]

  skip_if(is.na(test_file) || !file.exists(test_file),
          "Test file not found")

  # Test rmsana
  result_superassp <- superassp::rmsana(test_file, toFile = FALSE)
  result_wrassp <- wrassp::rmsana(test_file, toFile = FALSE)

  expect_s3_class(result_superassp, "AsspDataObj")
  expect_s3_class(result_wrassp, "AsspDataObj")
  # Compare track names (normalized to remove bracket notation and case)
  expect_equal(
    normalize_track_name(names(result_superassp)),
    normalize_track_name(names(result_wrassp))
  )
  expect_equal(dim(result_superassp[[1]]), dim(result_wrassp[[1]]))

  # Test acfana
  result_superassp <- superassp::acfana(test_file, toFile = FALSE)
  result_wrassp <- wrassp::acfana(test_file, toFile = FALSE)

  # Compare track names (normalized to remove bracket notation and case)
  expect_equal(
    normalize_track_name(names(result_superassp)),
    normalize_track_name(names(result_wrassp))
  )
  expect_equal(nrow(result_superassp[[1]]), nrow(result_wrassp[[1]]))
})

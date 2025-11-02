# Tests for dysprosody prosodic assessment functions
#
# DEPRECATED: Tests commented out because lst_dysprosody has been deprecated
# due to tuneR dependency removal.
#
# To re-enable these tests:
# 1. Re-implement lst_dysprosody without tuneR dependency
# 2. Uncomment the tests below
# 3. Run devtools::test()
#
# ============================================================================

# # Tests: lst_dysprosody(), install_dysprosody(), dysprosody_available(), dysprosody_info()
#
# test_that("dysprosody_available() returns logical", {
#   result <- dysprosody_available()
#   expect_type(result, "logical")
#   expect_length(result, 1)
# })
#
# test_that("dysprosody_info() returns NULL or list", {
#   if (dysprosody_available()) {
#     info <- dysprosody_info()
#     expect_type(info, "list")
#     expect_true("version" %in% names(info))
#     expect_true("optimized" %in% names(info))
#     expect_true("dependencies" %in% names(info))
#   } else {
#     info <- dysprosody_info()
#     expect_null(info)
#   }
# })
#
# test_that("lst_dysprosody works with single file", {
#   skip_if_not(dysprosody_available(), "dysprosody not available")
#
#   test_wav <- system.file("samples", "sustained", "a32b.wav", package = "superassp")
#   skip_if(test_wav == "", "Test file not found")
#
#   result <- lst_dysprosody(test_wav, verbose = FALSE)
#
#   # Check basic structure
#   expect_type(result, "list")
#   expect_true(length(result) > 0)
#
#   # Check required prosodic metadata features
#   expect_true("Duration" %in% names(result))
#   expect_true("PitchMean" %in% names(result))
#   expect_true("PitchRange" %in% names(result))
#   expect_true("PitchKey" %in% names(result))
#   expect_true("IntsIntLabels" %in% names(result))
#   expect_true("UniqueIntsInt" %in% names(result))
#
#   # Check feature count (should be 193)
#   expect_equal(length(result), 193)
#
#   # Check feature types
#   expect_type(result$Duration, "double")
#   expect_type(result$PitchMean, "double")
#   expect_type(result$IntsIntLabels, "double")
#
#   # Check Duration is reasonable
#   expect_gt(result$Duration, 0)
#   expect_lt(result$Duration, 10)  # Test file should be < 10 seconds
# })
#
# test_that("lst_dysprosody handles time windowing", {
#   skip_if_not(dysprosody_available(), "dysprosody not available")
#
#   test_wav <- system.file("samples", "sustained", "a32b.wav", package = "superassp")
#   skip_if(test_wav == "", "Test file not found")
#
#   # Get duration of full file
#   result_full <- lst_dysprosody(test_wav, verbose = FALSE)
#   skip_if(is.null(result_full), "Full file too short")
#
#   full_duration <- result_full$Duration
#
#   # Test time windowing
#   result_windowed <- lst_dysprosody(
#     test_wav,
#     beginTime = 0.5,
#     endTime = min(2.0, full_duration - 0.1),
#     verbose = FALSE
#   )
#
#   if (!is.null(result_windowed)) {
#     expect_type(result_windowed, "list")
#     expect_true("Duration" %in% names(result_windowed))
#     expect_lt(result_windowed$Duration, full_duration)
#     expect_gt(result_windowed$Duration, 0)
#   }
# })
#
# test_that("lst_dysprosody skips very short files", {
#   skip_if_not(dysprosody_available(), "dysprosody not available")
#
#   test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
#   skip_if(test_wav == "", "Test file not found")
#
#   # Try to process only 0.5 seconds (should be skipped as < 1 second)
#   result <- lst_dysprosody(
#     test_wav,
#     beginTime = 0,
#     endTime = 0.5,
#     verbose = FALSE
#   )
#
#   # Should return NULL for files < 1 second
#   expect_null(result)
# })
#
# test_that("lst_dysprosody handles batch processing (sequential)", {
#   skip_if_not(dysprosody_available(), "dysprosody not available")
#
#   test_files <- list.files(
#     system.file("samples", "sustained", package = "superassp"),
#     pattern = "\\.wav$",
#     full.names = TRUE
#   )
#
#   skip_if(length(test_files) < 2, "Not enough test files")
#
#   # Take first 3 files
#   files <- test_files[1:min(3, length(test_files))]
#
#   results <- lst_dysprosody(
#     files,
#     verbose = FALSE,
#     parallel = FALSE  # Sequential for predictability
#   )
#
#   expect_type(results, "list")
#   expect_true(length(results) > 0)
#
#   # Check each result
#   for (result in results) {
#     if (!is.null(result)) {
#       expect_type(result, "list")
#       expect_true("Duration" %in% names(result))
#       expect_true("PitchMean" %in% names(result))
#     }
#   }
# })
#
# test_that("lst_dysprosody handles missing files", {
#   skip_if_not(dysprosody_available(), "dysprosody not available")
#
#   expect_error(
#     lst_dysprosody("nonexistent_file.wav", verbose = FALSE),
#     regexp = "not found|does not exist",
#     ignore.case = TRUE
#   )
# })
#
# test_that("lst_dysprosody handles empty input", {
#   skip_if_not(dysprosody_available(), "dysprosody not available")
#
#   expect_error(
#     lst_dysprosody(character(0), verbose = FALSE),
#     regexp = "No files|empty"
#   )
# })
#
# test_that("lst_dysprosody validates time parameters", {
#   skip_if_not(dysprosody_available(), "dysprosody not available")
#
#   test_wav <- system.file("samples", "sustained", "a32b.wav", package = "superassp")
#   skip_if(test_wav == "", "Test file not found")
#
#   # Multiple files with mismatched time parameters should error
#   files <- rep(test_wav, 3)
#   expect_error(
#     lst_dysprosody(files, beginTime = c(0, 0.5), verbose = FALSE),
#     regexp = "same length|must match"
#   )
# })
#
# test_that("lst_dysprosody handles custom F0 range", {
#   skip_if_not(dysprosody_available(), "dysprosody not available")
#
#   test_wav <- system.file("samples", "sustained", "a32b.wav", package = "superassp")
#   skip_if(test_wav == "", "Test file not found")
#
#   # Test with custom F0 range
#   result <- lst_dysprosody(
#     test_wav,
#     minF = 100,
#     maxF = 500,
#     verbose = FALSE
#   )
#
#   if (!is.null(result)) {
#     expect_type(result, "list")
#     expect_true("PitchMean" %in% names(result))
#
#     # Pitch mean should be reasonable for the given range
#     if (!is.na(result$PitchMean) && result$PitchMean > 0) {
#       expect_gte(result$PitchMean, 50)   # Allow some margin
#       expect_lte(result$PitchMean, 600)
#     }
#   }
# })
#
# test_that("lst_dysprosody returns consistent feature names", {
#   skip_if_not(dysprosody_available(), "dysprosody not available")
#
#   test_files <- list.files(
#     system.file("samples", "sustained", package = "superassp"),
#     pattern = "\\.wav$",
#     full.names = TRUE
#   )[1:2]
#
#   skip_if(length(test_files) < 2, "Not enough test files")
#
#   results <- lst_dysprosody(test_files, verbose = FALSE, parallel = FALSE)
#
#   # Get non-NULL results
#   valid_results <- Filter(Negate(is.null), results)
#   skip_if(length(valid_results) < 2, "Not enough valid results")
#
#   # Check feature names are consistent across files
#   names1 <- names(valid_results[[1]])
#   names2 <- names(valid_results[[2]])
#
#   expect_equal(length(names1), length(names2))
#   expect_equal(sort(names1), sort(names2))
# })
#
# test_that("lst_dysprosody spectral features are present", {
#   skip_if_not(dysprosody_available(), "dysprosody not available")
#
#   test_wav <- system.file("samples", "sustained", "a32b.wav", package = "superassp")
#   skip_if(test_wav == "", "Test file not found")
#
#   result <- lst_dysprosody(test_wav, verbose = FALSE)
#   skip_if(is.null(result), "File too short")
#
#   # Check spectral tilt features
#   expect_true("L2L1_mean" %in% names(result) ||
#               "L2L1" %in% names(result))
#   expect_true("SLF_mean" %in% names(result) ||
#               "SLF" %in% names(result))
#   expect_true("C1_mean" %in% names(result) ||
#               "C1" %in% names(result))
# })
#
# test_that("lst_dysprosody handles non-WAV formats", {
#   skip_if_not(dysprosody_available(), "dysprosody not available")
#
#   # This test would require MP3/MP4 test files
#   # Skip if not available
#   skip("Non-WAV format test requires MP3/MP4 test files")
# })
#
# test_that("lst_dysprosody parallel processing works", {
#   skip_if_not(dysprosody_available(), "dysprosody not available")
#   skip_on_cran()  # Parallel processing may be unstable on CRAN
#
#   test_files <- list.files(
#     system.file("samples", "sustained", package = "superassp"),
#     pattern = "\\.wav$",
#     full.names = TRUE
#   )[1:3]
#
#   skip_if(length(test_files) < 3, "Not enough test files")
#
#   # Test parallel processing
#   results_parallel <- lst_dysprosody(
#     test_files,
#     verbose = FALSE,
#     parallel = TRUE,
#     n_cores = 2
#   )
#
#   # Test sequential processing
#   results_sequential <- lst_dysprosody(
#     test_files,
#     verbose = FALSE,
#     parallel = FALSE
#   )
#
#   # Results should be similar (allowing for minor numerical differences)
#   expect_equal(length(results_parallel), length(results_sequential))
#
#   # Check that at least some files were processed
#   n_valid_parallel <- sum(sapply(results_parallel, function(x) !is.null(x)))
#   n_valid_sequential <- sum(sapply(results_sequential, function(x) !is.null(x)))
#
#   expect_equal(n_valid_parallel, n_valid_sequential)
# })

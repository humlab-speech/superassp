test_that("lst_voice_report works with single file", {
  skip_if_not(pladdrr_available(), "pladdrr not available")
  
  test_file <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_file == "", "Test file not found")
  
  result <- lst_voice_report(test_file, verbose = FALSE)
  
  # Check structure
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1)
  expect_equal(ncol(result), 31)  # file + 30 measures
  
  # Check column names
  expected_cols <- c("file", "start_time", "end_time", "selection_start", "selection_end",
                     "median_pitch", "mean_pitch", "sd_pitch", "min_pitch", "max_pitch",
                     "num_pulses", "num_periods", "mean_period", "sd_period",
                     "fraction_unvoiced", "num_voice_breaks", "degree_voice_breaks",
                     "jitter_local_percent", "jitter_local_abs", "jitter_rap_percent",
                     "jitter_ppq5_percent", "jitter_ddp_percent",
                     "shimmer_local_percent", "shimmer_local_db", "shimmer_apq3_percent",
                     "shimmer_apq5_percent", "shimmer_apq11_percent", "shimmer_dda_percent",
                     "mean_autocorrelation", "mean_nhr", "mean_hnr")
  expect_equal(names(result), expected_cols)
  
  # Check reasonable values for sustained /a/
  expect_gt(result$mean_pitch, 50)
  expect_lt(result$mean_pitch, 500)
  expect_gt(result$mean_hnr, 10)  # Sustained vowel should have good HNR
  expect_lt(result$jitter_local_percent, 5)  # Low jitter expected
  expect_lt(result$shimmer_local_percent, 10)  # Low shimmer expected
})

test_that("lst_voice_report works with multiple files", {
  skip_if_not(pladdrr_available(), "pladdrr not available")
  
  sustained_dir <- system.file("samples", "sustained", package = "superassp")
  files <- list.files(sustained_dir, pattern = ".wav", full.names = TRUE)
  skip_if(length(files) < 2, "Need at least 2 test files")
  
  results <- lst_voice_report(files[1:2], verbose = FALSE)
  
  expect_s3_class(results, "data.frame")
  expect_equal(nrow(results), 2)
  expect_equal(ncol(results), 31)
  expect_true(all(!is.na(results$mean_hnr)))
})

test_that("lst_voice_report supports time windowing", {
  skip_if_not(pladdrr_available(), "pladdrr not available")
  
  test_file <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_file == "", "Test file not found")
  
  # Full file
  result_full <- lst_voice_report(test_file, verbose = FALSE)
  
  # Windowed (1.0-2.0s)
  result_window <- lst_voice_report(test_file, beginTime = 1.0, endTime = 2.0, verbose = FALSE)
  
  expect_equal(result_window$start_time, 1.0)
  expect_equal(result_window$end_time, 2.0)
  
  # Pitch should be similar but not identical
  expect_lt(abs(result_window$mean_pitch - result_full$mean_pitch), 20)
})

test_that("lst_voice_report supports selection offset/length", {
  skip_if_not(pladdrr_available(), "pladdrr not available")
  
  test_file <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_file == "", "Test file not found")
  
  # Selection: 0.5s offset, 1.0s length
  result <- lst_voice_report(test_file, 
                               selectionOffset = 0.5, 
                               selectionLength = 1.0, 
                               verbose = FALSE)
  
  expect_equal(result$selection_start, 0.5)
  expect_equal(result$selection_end, 1.5)
})

test_that("lst_voice_report writes JSTF files", {
  skip_if_not(pladdrr_available(), "pladdrr not available")
  
  test_file <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_file == "", "Test file not found")
  
  temp_dir <- tempdir()
  
  output <- lst_voice_report(test_file, toFile = TRUE, 
                               outputDirectory = temp_dir, 
                               verbose = FALSE)
  
  expect_true(file.exists(output))
  expect_match(output, "\\.pvr$")
  
  # Read back
  track <- read_jstf(output)
  expect_s3_class(track, "JsonTrackObj")
  
  # Convert to data.frame
  df <- as.data.frame(track)
  expect_s3_class(df, "data.frame")
  expect_equal(nrow(df), 1)
  expect_true("median_pitch" %in% names(df))
  expect_true("mean_hnr" %in% names(df))
  
  # Check values match in-memory version
  result_mem <- lst_voice_report(test_file, verbose = FALSE)
  expect_equal(df$median_pitch, result_mem$median_pitch, tolerance = 0.01)
  expect_equal(df$mean_hnr, result_mem$mean_hnr, tolerance = 0.01)
  
  # Cleanup
  unlink(output)
})

test_that("lst_voice_report handles pitch parameters", {
  skip_if_not(pladdrr_available(), "pladdrr not available")
  
  test_file <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_file == "", "Test file not found")
  
  # Test with different F0 range
  result_low <- lst_voice_report(test_file, minF = 50, maxF = 300, verbose = FALSE)
  result_high <- lst_voice_report(test_file, minF = 100, maxF = 600, verbose = FALSE)
  
  # Both should work
  expect_s3_class(result_low, "data.frame")
  expect_s3_class(result_high, "data.frame")
  
  # Values might differ slightly
  expect_true(is.numeric(result_low$mean_pitch))
  expect_true(is.numeric(result_high$mean_pitch))
})

test_that("lst_voice_report validates inputs", {
  skip_if_not(pladdrr_available(), "pladdrr not available")
  
  # Invalid window shape
  test_file <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_file == "", "Test file not found")
  
  expect_error(
    lst_voice_report(test_file, windowShape = "invalid"),
    "Invalid windowShape"
  )
  
  # Missing file
  expect_error(
    lst_voice_report("nonexistent.wav"),
    "Files not found"
  )
})

# ---- Task 4: return_jstf parameter ----

test_that("lst_voice_report return_jstf=TRUE returns JsonTrackObj", {
  skip_if_not(pladdrr_available(), "pladdrr not available")
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "test wav not found")
  result <- lst_voice_report(test_wav, return_jstf = TRUE, verbose = FALSE)
  expect_s3_class(result, "JsonTrackObj")
  expect_equal(result$function_name, "lst_voice_report")
  expect_true("median_pitch" %in% track_names(result))
})

test_that("lst_voice_report return_jstf=TRUE multi-file returns list of JsonTrackObj", {
  skip_if_not(pladdrr_available(), "pladdrr not available")
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "test wav not found")
  results <- lst_voice_report(c(test_wav, test_wav), return_jstf = TRUE, verbose = FALSE)
  expect_type(results, "list")
  expect_s3_class(results[[1]], "JsonTrackObj")
  expect_s3_class(results[[2]], "JsonTrackObj")
})

test_that("lst_voice_report toFile=TRUE + return_jstf=TRUE writes file and returns obj", {
  skip_if_not(pladdrr_available(), "pladdrr not available")
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "test wav not found")
  tmp <- tempdir()
  result <- lst_voice_report(test_wav, toFile = TRUE, return_jstf = TRUE,
                              outputDirectory = tmp, verbose = FALSE)
  expect_s3_class(result, "JsonTrackObj")
  expect_true(file.exists(file.path(tmp, "a1.pvr")))
})

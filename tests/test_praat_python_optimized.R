#!/usr/bin/env Rscript
# Test script for Python/Parselmouth-based Praat implementations

library(superassp)

# Source the optimized functions
source("../R/praat_python_optimized.R")

# Check if reticulate and parselmouth are available
if (!requireNamespace("reticulate", quietly = TRUE)) {
  cat("reticulate package not available - skipping Python tests\n")
  quit(status = 0)
}

# Use system Python3
Sys.unsetenv("RETICULATE_PYTHON")
reticulate::use_python("/usr/bin/python3", required = FALSE)

# Try to source a Python script to check if environment works
tryCatch({
  reticulate::py_run_string("import parselmouth")
  cat("✓ Parselmouth is available\n\n")
}, error = function(e) {
  cat("Python/Parselmouth not available or not properly configured\n")
  cat("Error:", conditionMessage(e), "\n")
  cat("Please install: pip3 install praat-parselmouth\n")
  cat("Skipping Python/Parselmouth tests\n")
  quit(status = 0)
})

test_file <- "signalfiles/AVQI/input/sv1.wav"

cat("Testing Python/Parselmouth-based Praat implementations\n")
cat("=======================================================\n\n")

# Test 1: praat_formant_burg_opt
cat("Test 1: praat_formant_burg_opt - Formant analysis\n")
cat("-------------------------------------------------\n")
tryCatch({
  result1 <- praat_formant_burg_opt(
    test_file,
    beginTime = 0.0,
    endTime = 0.0,
    timeStep = 0.005,  # Use timeStep not time_step
    number_of_formants = 5,
    maxHzFormant = 5500.0,
    toFile = FALSE
  )
  cat("✓ Function executed successfully\n")
  cat("  Result class:", class(result1), "\n")
  cat("  Tracks:", names(result1), "\n")
  if (length(result1) > 0) {
    cat("  Track 1 dimensions:", dim(result1[[1]]), "\n")
  }
}, error = function(e) {
  cat("✗ Error:", conditionMessage(e), "\n")
})
cat("\n")

# Test 2: praat_pitch_opt
cat("Test 2: praat_pitch_opt - Pitch tracking\n")
cat("-----------------------------------------\n")
tryCatch({
  result2 <- praat_pitch_opt(
    test_file,
    beginTime = 0.0,
    endTime = 0.0,
    timeStep = 0.005,
    minF = 75.0,
    maxF = 600.0,
    toFile = FALSE
  )
  cat("✓ Function executed successfully\n")
  cat("  Result class:", class(result2), "\n")
  cat("  Tracks:", names(result2), "\n")
  if (length(result2) > 0) {
    cat("  Track 1 dimensions:", dim(result2[[1]]), "\n")
  }
}, error = function(e) {
  cat("✗ Error:", conditionMessage(e), "\n")
})
cat("\n")

# Test 3: praat_intensity_opt
cat("Test 3: praat_intensity_opt - Intensity analysis\n")
cat("-------------------------------------------------\n")
tryCatch({
  result3 <- praat_intensity_opt(
    test_file,
    beginTime = 0.0,
    endTime = 0.0,
    minF = 100.0,
    toFile = FALSE
  )
  cat("✓ Function executed successfully\n")
  cat("  Result class:", class(result3), "\n")
  cat("  Tracks:", names(result3), "\n")
  if (length(result3) > 0) {
    cat("  Track 1 dimensions:", dim(result3[[1]]), "\n")
  }
}, error = function(e) {
  cat("✗ Error:", conditionMessage(e), "\n")
})
cat("\n")

# Test 4: praat_spectral_moments_opt
cat("Test 4: praat_spectral_moments_opt - Spectral shape\n")
cat("----------------------------------------------------\n")
tryCatch({
  result4 <- praat_spectral_moments_opt(
    test_file,
    beginTime = 0.0,
    endTime = 0.0,
    frameLength = 0.030,
    toFile = FALSE
  )
  cat("✓ Function executed successfully\n")
  cat("  Result class:", class(result4), "\n")
  cat("  Tracks:", names(result4), "\n")
  if (length(result4) > 0) {
    cat("  Track 1 dimensions:", dim(result4[[1]]), "\n")
  }
}, error = function(e) {
  cat("✗ Error:", conditionMessage(e), "\n")
})
cat("\n")

# Test 5: praat_formantpath_burg_opt
cat("Test 5: praat_formantpath_burg_opt - FormantPath\n")
cat("-------------------------------------------------\n")
tryCatch({
  result5 <- praat_formantpath_burg_opt(
    test_file,
    beginTime = 0.0,
    endTime = 0.0,
    timeStep = 0.005,
    maxNumFormants = 5,
    toFile = FALSE
  )
  cat("✓ Function executed successfully\n")
  cat("  Result class:", class(result5), "\n")
  cat("  Tracks:", names(result5), "\n")
  if (length(result5) > 0) {
    cat("  Track 1 dimensions:", dim(result5[[1]]), "\n")
  }
}, error = function(e) {
  cat("✗ Error:", conditionMessage(e), "\n")
})
cat("\n")

# Test 6: Time windowing
cat("Test 6: Time windowing with praat_pitch_opt\n")
cat("--------------------------------------------\n")
tryCatch({
  result6 <- praat_pitch_opt(
    test_file,
    beginTime = 0.5,
    endTime = 1.5,
    timeStep = 0.005,
    toFile = FALSE
  )
  cat("✓ Time windowing works\n")
  cat("  Track 1 dimensions:", dim(result6[[1]]), "\n")
  cat("  Expected ~200 frames for 1 second at 200Hz\n")
}, error = function(e) {
  cat("✗ Error:", conditionMessage(e), "\n")
})
cat("\n")

# Test 7: File output
cat("Test 7: Writing output to file\n")
cat("-------------------------------\n")
temp_dir <- tempdir()
tryCatch({
  result7 <- praat_intensity_opt(
    test_file,
    toFile = TRUE,
    outputDirectory = temp_dir
  )
  output_file <- file.path(temp_dir, paste0(tools::file_path_sans_ext(basename(test_file)), ".int"))
  if (file.exists(output_file)) {
    cat("✓ File written successfully\n")
    cat("  Output file:", output_file, "\n")
    cat("  File size:", file.size(output_file), "bytes\n")
    unlink(output_file)
  } else {
    cat("✗ Output file not found\n")
  }
}, error = function(e) {
  cat("✗ Error:", conditionMessage(e), "\n")
})
cat("\n")

cat("All tests completed!\n")
cat("====================\n")

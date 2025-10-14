#!/usr/bin/env Rscript
# Test script to verify Python/Parselmouth implementations work correctly
# and produce valid output

library(superassp)

# Check if reticulate and parselmouth are available
if (!requireNamespace("reticulate", quietly = TRUE)) {
  cat("reticulate package not available - skipping Python tests\n")
  quit(status = 0)
}

# Use system Python3
Sys.unsetenv("RETICULATE_PYTHON")
reticulate::use_python("/usr/bin/python3", required = FALSE)

# Try to load parselmouth
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

if (!file.exists(test_file)) {
  cat("Test file not found:", test_file, "\n")
  quit(status = 1)
}

cat("Testing Python/Parselmouth-based Praat implementations\n")
cat("=======================================================\n\n")

# Test 1: praat_formant_burg_opt
cat("Test 1: praat_formant_burg_opt - Basic formant analysis\n")
cat("--------------------------------------------------------\n")
tryCatch({
  result1 <- praat_formant_burg_opt(
    test_file,
    beginTime = 0.0,
    endTime = 0.0,
    timeStep = 0.005,
    number_of_formants = 5,
    maxHzFormant = 5500.0,
    toFile = FALSE
  )
  cat("✓ Function executed successfully\n")
  cat("  Result class:", class(result1), "\n")
  cat("  Tracks:", names(result1), "\n")
  if (length(result1) > 0) {
    cat("  Track 'fm1' dimensions:", dim(result1$fm1), "\n")
    cat("  Sample rate:", attr(result1, "sampleRate"), "\n")
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
    time_step = 0.005,
    window_length = 0.040,
    minimum_f0 = 75.0,
    maximum_f0 = 600.0,
    toFile = FALSE
  )
  cat("✓ Function executed successfully\n")
  cat("  Result class:", class(result2), "\n")
  cat("  Tracks:", names(result2), "\n")
  if (length(result2) > 0) {
    cat("  Track 'cc' dimensions:", dim(result2$cc), "\n")
    cat("  Sample rate:", attr(result2, "sampleRate"), "\n")
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
    time_step = 0.0,  # Automatic
    minimal_f0_frequency = 50.0,
    toFile = FALSE
  )
  cat("✓ Function executed successfully\n")
  cat("  Result class:", class(result3), "\n")
  cat("  Tracks:", names(result3), "\n")
  if (length(result3) > 0) {
    cat("  Track 'intensity' dimensions:", dim(result3$intensity), "\n")
    cat("  Sample rate:", attr(result3, "sampleRate"), "\n")
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
    windowLength = 0.005,
    time_step = 0.005,
    toFile = FALSE
  )
  cat("✓ Function executed successfully\n")
  cat("  Result class:", class(result4), "\n")
  cat("  Tracks:", names(result4), "\n")
  if (length(result4) > 0) {
    cat("  Track 'cog' dimensions:", dim(result4$cog), "\n")
    cat("  Sample rate:", attr(result4, "sampleRate"), "\n")
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
    time_step = 0.005,
    number_of_formants = 5,
    maxHzFormant = 5500.0,
    toFile = FALSE
  )
  cat("✓ Function executed successfully\n")
  cat("  Result class:", class(result5), "\n")
  cat("  Tracks:", names(result5), "\n")
  if (length(result5) > 0) {
    cat("  Track 'fm1' dimensions:", dim(result5$fm1), "\n")
    cat("  Sample rate:", attr(result5, "sampleRate"), "\n")
  }
}, error = function(e) {
  cat("✗ Error:", conditionMessage(e), "\n")
})
cat("\n")

# Test 6: Time windowing with praat_pitch_opt
cat("Test 6: Time windowing with praat_pitch_opt\n")
cat("--------------------------------------------\n")
tryCatch({
  result6 <- praat_pitch_opt(
    test_file,
    beginTime = 0.5,
    endTime = 1.5,
    time_step = 0.005,
    minimum_f0 = 75.0,
    maximum_f0 = 600.0,
    toFile = FALSE
  )
  cat("✓ Time windowing works\n")
  cat("  Track 'cc' dimensions:", dim(result6$cc), "\n")
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
  num_files <- praat_intensity_opt(
    test_file,
    toFile = TRUE,
    outputDirectory = temp_dir
  )
  output_file <- file.path(temp_dir, "sv1.int")
  if (file.exists(output_file)) {
    cat("✓ File written successfully\n")
    cat("  Output file:", output_file, "\n")
    cat("  File size:", file.size(output_file), "bytes\n")
    cat("  Files processed:", num_files, "\n")

    # Read it back to verify
    read_obj <- wrassp::read.AsspDataObj(output_file)
    cat("  Read back - tracks:", names(read_obj), "\n")

    unlink(output_file)
  } else {
    cat("✗ Output file not found\n")
  }
}, error = function(e) {
  cat("✗ Error:", conditionMessage(e), "\n")
})
cat("\n")

cat("All basic tests completed!\n")
cat("==========================\n")

#!/usr/bin/env Rscript
# Test script for av package integration and memory-based processing

library(superassp)

# Source the updated av_helpers functions
source("../R/av_helpers.R")

test_file <- "signalfiles/AVQI/input/sv1.wav"

cat("Testing av integration and memory-based processing\n")
cat("===================================================\n\n")

# Test 1: Traditional file-based processing
cat("Test 1: File-based RMS analysis (baseline)\n")
cat("-------------------------------------------\n")
result_file <- rmsana(test_file, toFile = FALSE)
cat("Result class:", class(result_file), "\n")
cat("Result structure:\n")
str(result_file, max.level = 1)
cat("\n")

# Test 2: Check if av package is available
cat("Test 2: Check av package availability\n")
cat("--------------------------------------\n")
av_available <- requireNamespace("av", quietly = TRUE)
if (av_available) {
  cat("av package is available\n")

  # Test 3: Convert audio using av
  cat("\nTest 3: Convert WAV to AsspDataObj using av\n")
  cat("--------------------------------------------\n")
  tryCatch({
    audio_obj <- av_to_asspDataObj(test_file)
    cat("Successfully created AsspDataObj from av\n")
    cat("Sample rate:", attr(audio_obj, "sampleRate"), "Hz\n")
    cat("Number of frames:", attr(audio_obj, "endRecord") - attr(audio_obj, "startRecord") + 1, "\n")
    cat("Track format:", attr(audio_obj, "trackFormats"), "\n")

    # Test 4: Process using memory-based analysis
    cat("\nTest 4: Memory-based RMS analysis\n")
    cat("-----------------------------------\n")
    result_memory <- rmsana_memory(audio_obj, windowShift = 5)
    cat("Successfully performed memory-based analysis\n")
    cat("Result class:", class(result_memory), "\n")
    str(result_memory, max.level = 1)

    # Test 5: Compare results
    cat("\nTest 5: Compare file-based vs memory-based results\n")
    cat("---------------------------------------------------\n")
    cat("File-based dims:", dim(result_file[[1]]), "\n")
    cat("Memory-based dims:", dim(result_memory[[1]]), "\n")

    # Since we're using the same parameters, results should be identical
    if (identical(dim(result_file[[1]]), dim(result_memory[[1]]))) {
      cat("✓ Dimensions match\n")

      # Compare values (allowing for small numerical differences)
      max_diff <- max(abs(result_file[[1]] - result_memory[[1]]), na.rm = TRUE)
      cat("Maximum difference:", max_diff, "\n")

      if (max_diff < 1e-6) {
        cat("✓ Results are numerically identical\n")
      } else if (max_diff < 0.01) {
        cat("✓ Results are very similar (within tolerance)\n")
      } else {
        cat("✗ Results differ significantly\n")
      }
    } else {
      cat("✗ Dimensions do not match\n")
    }

    # Test 6: Process with different time window
    cat("\nTest 6: Process specific time segment\n")
    cat("--------------------------------------\n")
    audio_segment <- av_to_asspDataObj(test_file, start_time = 0.1, end_time = 0.5)
    result_segment <- rmsana_memory(audio_segment, windowShift = 5)
    cat("Processed segment from 0.1s to 0.5s\n")
    cat("Number of frames in segment:", nrow(result_segment[[1]]), "\n")

    # Test 7: Test convenience function
    cat("\nTest 7: Test process_media_file convenience function\n")
    cat("-----------------------------------------------------\n")
    result_convenience <- process_media_file(test_file, "rmsana", windowShift = 5)
    cat("Successfully used convenience function\n")
    cat("Result dimensions:", dim(result_convenience[[1]]), "\n")

  }, error = function(e) {
    cat("✗ Error in av-based processing:", conditionMessage(e), "\n")
  })

} else {
  cat("av package is NOT available\n")
  cat("Install with: install.packages('av')\n")
  cat("Skipping av-based tests\n")
}

cat("\n")
cat("Test suite completed\n")
cat("====================\n")

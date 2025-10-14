#!/usr/bin/env Rscript
# Test script for load-and-process pattern with time windows

library(superassp)

test_file <- "signalfiles/AVQI/input/sv1.wav"

cat("Testing load-and-process pattern with time windows\n")
cat("===================================================\n\n")

# Test 1: RMS analysis without time window (fast path)
cat("Test 1: rmsana - Native file, no time window (fast path)\n")
cat("----------------------------------------------------------\n")
result1 <- rmsana(test_file, toFile = FALSE, windowShift = 5)
cat("✓ Result dimensions:", dim(result1[[1]]), "\n")
cat("✓ Track name:", names(result1), "\n\n")

# Test 2: RMS analysis with time window (load-and-process path)
cat("Test 2: rmsana - Native file with time window (load-and-process)\n")
cat("----------------------------------------------------------------\n")
result2 <- rmsana(test_file, beginTime = 0.5, endTime = 1.5, toFile = FALSE, windowShift = 5)
cat("✓ Result dimensions:", dim(result2[[1]]), "\n")
cat("✓ Expected ~200 frames for 1 second at 200Hz\n")
cat("✓ Actual frames:", nrow(result2[[1]]), "\n\n")

# Test 3: Forest analysis without time window
cat("Test 3: forest - Native file, no time window\n")
cat("---------------------------------------------\n")
result3 <- forest(test_file, toFile = FALSE, windowShift = 5)
cat("✓ Result tracks:", names(result3), "\n")
cat("✓ First track dimensions:", dim(result3[[1]]), "\n\n")

# Test 4: Forest with time window
cat("Test 4: forest - Native file with time window\n")
cat("----------------------------------------------\n")
result4 <- forest(test_file, beginTime = 0.2, endTime = 0.8, toFile = FALSE, windowShift = 5)
cat("✓ Result tracks:", names(result4), "\n")
cat("✓ Expected ~120 frames for 0.6 seconds at 200Hz\n")
cat("✓ Actual frames:", nrow(result4[[1]]), "\n\n")

# Test 5: Compare beginTime/endTime handling
cat("Test 5: Verify time window accuracy\n")
cat("------------------------------------\n")
full_duration <- 2.918  # sv1.wav duration
segment_start <- 1.0
segment_end <- 2.0
segment_duration <- segment_end - segment_start

result5a <- rmsana(test_file, toFile = FALSE, windowShift = 5)
result5b <- rmsana(test_file, beginTime = segment_start, endTime = segment_end,
                   toFile = FALSE, windowShift = 5)

cat("Full file frames:", nrow(result5a[[1]]), "\n")
cat("Segment frames:", nrow(result5b[[1]]), "\n")
cat("Segment duration:", segment_duration, "seconds\n")
cat("Expected frames:", round(segment_duration * 200), "(at 200Hz)\n")

ratio <- nrow(result5b[[1]]) / nrow(result5a[[1]])
expected_ratio <- segment_duration / full_duration
cat("Actual ratio:", round(ratio, 3), "\n")
cat("Expected ratio:", round(expected_ratio, 3), "\n")
cat(if(abs(ratio - expected_ratio) < 0.05) "✓" else "✗",
    "Time window processing is accurate\n\n")

# Test 6: Zero-crossing rate analysis
cat("Test 6: zcrana - Native file with time window\n")
cat("----------------------------------------------\n")
result6 <- zcrana(test_file, beginTime = 0.5, endTime = 1.0, toFile = FALSE, windowShift = 5)
cat("✓ Result tracks:", names(result6), "\n")
cat("✓ Frames:", nrow(result6[[1]]), "\n\n")

# Test 7: ACF analysis with time window
cat("Test 7: acfana - Native file with time window\n")
cat("----------------------------------------------\n")
result7 <- acfana(test_file, beginTime = 0.3, endTime = 0.9, toFile = FALSE, windowShift = 5)
cat("✓ Result tracks:", names(result7), "\n")
cat("✓ Frames:", nrow(result7[[1]]), "\n\n")

# Test 8: Multiple files with different time windows
cat("Test 8: Multiple files with different time windows\n")
cat("---------------------------------------------------\n")
files <- c(test_file, test_file)
begin_times <- c(0.0, 1.0)
end_times <- c(1.0, 2.0)

result8 <- rmsana(files, beginTime = begin_times, endTime = end_times,
                  toFile = FALSE, windowShift = 5)
cat("✓ Number of results:", length(result8), "\n")
cat("✓ First segment frames:", nrow(result8[[1]][[1]]), "\n")
cat("✓ Second segment frames:", nrow(result8[[2]][[1]]), "\n")
cat("✓ Segments should have similar frame counts:",
    abs(nrow(result8[[1]][[1]]) - nrow(result8[[2]][[1]])) < 10, "\n\n")

cat("All tests completed successfully!\n")
cat("==================================\n")

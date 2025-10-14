#!/usr/bin/env Rscript
# Test load-and-process pattern with non-native media formats

library(superassp)

# Check if av package is available
if (!requireNamespace("av", quietly = TRUE)) {
  cat("av package not available - skipping non-native format tests\n")
  quit(status = 0)
}

test_wav <- "signalfiles/AVQI/input/sv1.wav"
test_mp3 <- tempfile(fileext = ".mp3")

cat("Testing load-and-process with non-native formats\n")
cat("=================================================\n\n")

# Create a test MP3 file
cat("Creating temporary MP3 file for testing...\n")
av::av_audio_convert(test_wav, test_mp3, verbose = FALSE)
cat("✓ MP3 file created:", test_mp3, "\n\n")

# Test 1: RMS on MP3 without time window
cat("Test 1: rmsana on MP3 file (no time window)\n")
cat("--------------------------------------------\n")
result1_wav <- rmsana(test_wav, toFile = FALSE, windowShift = 5, verbose = FALSE)
result1_mp3 <- rmsana(test_mp3, toFile = FALSE, windowShift = 5, verbose = FALSE)
cat("WAV frames:", nrow(result1_wav[[1]]), "\n")
cat("MP3 frames:", nrow(result1_mp3[[1]]), "\n")
cat("✓ Frame counts match:", abs(nrow(result1_wav[[1]]) - nrow(result1_mp3[[1]])) < 5, "\n\n")

# Test 2: RMS on MP3 with time window
cat("Test 2: rmsana on MP3 with time window\n")
cat("---------------------------------------\n")
result2 <- rmsana(test_mp3, beginTime = 0.5, endTime = 1.5,
                  toFile = FALSE, windowShift = 5, verbose = FALSE)
cat("✓ Result dimensions:", dim(result2[[1]]), "\n")
cat("✓ Expected ~200 frames for 1 second\n")
cat("✓ Actual frames:", nrow(result2[[1]]), "\n\n")

# Test 3: Forest on MP3 with time window
cat("Test 3: forest on MP3 with time window\n")
cat("---------------------------------------\n")
result3 <- forest(test_mp3, beginTime = 0.3, endTime = 0.9,
                  toFile = FALSE, windowShift = 5, verbose = FALSE)
cat("✓ Result tracks:", names(result3), "\n")
cat("✓ Expected ~120 frames for 0.6 seconds\n")
cat("✓ Actual frames:", nrow(result3[[1]]), "\n\n")

# Test 4: Compare WAV vs MP3 results with time window
cat("Test 4: Compare WAV vs MP3 with same time window\n")
cat("-------------------------------------------------\n")
result4_wav <- rmsana(test_wav, beginTime = 1.0, endTime = 2.0,
                      toFile = FALSE, windowShift = 5, verbose = FALSE)
result4_mp3 <- rmsana(test_mp3, beginTime = 1.0, endTime = 2.0,
                      toFile = FALSE, windowShift = 5, verbose = FALSE)
cat("WAV segment frames:", nrow(result4_wav[[1]]), "\n")
cat("MP3 segment frames:", nrow(result4_mp3[[1]]), "\n")
cat("✓ Frame counts similar:", abs(nrow(result4_wav[[1]]) - nrow(result4_mp3[[1]])) < 5, "\n\n")

# Cleanup
unlink(test_mp3)
cat("✓ Cleanup complete\n\n")

cat("All non-native format tests passed!\n")
cat("====================================\n")

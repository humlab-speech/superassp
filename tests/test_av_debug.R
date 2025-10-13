#!/usr/bin/env Rscript
# Detailed debug test for av integration

library(superassp)

test_file <- "signalfiles/AVQI/input/sv1.wav"

cat("Step 1: Check av package\n")
if (!requireNamespace("av", quietly = TRUE)) {
  cat("av package not available - stopping\n")
  quit(status = 1)
}
cat("✓ av package available\n\n")

cat("Step 2: Read audio info using av\n")
info <- av::av_media_info(test_file)
cat("Audio info:\n")
print(info$audio)
cat("\n")

cat("Step 3: Read audio data using av\n")
audio_data <- av::read_audio_bin(test_file, start_time = 0, end_time = NULL,
                                  sample_rate = info$audio$sample_rate)
cat("Audio data length:", length(audio_data), "bytes\n")
cat("Audio data class:", class(audio_data), "\n\n")

cat("Step 4: Scale from 32-bit to 16-bit\n")
# av::read_audio_bin returns 32-bit signed integers
# We need to scale to 16-bit range for AsspDataObj
samples_int16 <- as.integer(audio_data / 65536)
cat("Number of samples:", length(samples_int16), "\n")
cat("Sample range (32-bit):", range(audio_data), "\n")
cat("Sample range (16-bit):", range(samples_int16), "\n\n")

cat("Step 5: Create sample matrix\n")
channels <- info$audio$channels
if (channels > 1) {
  n_frames <- length(samples_int16) / channels
  sample_matrix <- matrix(samples_int16, nrow = n_frames, ncol = channels, byrow = TRUE)
} else {
  sample_matrix <- matrix(samples_int16, ncol = 1)
}
cat("Matrix dimensions:", dim(sample_matrix), "\n\n")

cat("Step 6: Create AsspDataObj\n")
result <- list()
result$audio <- sample_matrix

n_frames <- nrow(sample_matrix)
attr(result, "sampleRate") <- as.numeric(info$audio$sample_rate)
attr(result, "trackFormats") <- "INT16"
attr(result, "startTime") <- 0.0
attr(result, "startRecord") <- 1L
attr(result, "endRecord") <- as.integer(n_frames)
attr(result, "origFreq") <- 0.0  # numeric 0 for converted audio
attr(result, "filePath") <- test_file
attr(result, "fileInfo") <- c(21L, 2L)  # FF_WAVE=21, FDF_BIN=2

class(result) <- "AsspDataObj"

cat("AsspDataObj created\n")
cat("Class:", class(result), "\n")
cat("Names:", names(result), "\n")
cat("Sample rate:", attr(result, "sampleRate"), "\n")
cat("Track formats:", attr(result, "trackFormats"), "\n")
cat("Start/end record:", attr(result, "startRecord"), "/", attr(result, "endRecord"), "\n")
cat("File info:", attr(result, "fileInfo"), "\n\n")

cat("Step 7: Compare with file-based read\n")
obj_from_file <- read.AsspDataObj(test_file)
cat("File-based object:\n")
cat("  Class:", class(obj_from_file), "\n")
cat("  Names:", names(obj_from_file), "\n")
cat("  Sample rate:", attr(obj_from_file, "sampleRate"), "\n")
cat("  Track formats:", attr(obj_from_file, "trackFormats"), "\n")
cat("  Start/end record:", attr(obj_from_file, "startRecord"), "/", attr(obj_from_file, "endRecord"), "\n")
cat("  Attributes:\n")
print(names(attributes(obj_from_file)))
cat("\n")

cat("Step 8: Write av-based AsspDataObj to temp file\n")
temp_wav <- tempfile(fileext = ".wav")
cat("Temp file:", temp_wav, "\n")
tryCatch({
  write.AsspDataObj(result, temp_wav)
  cat("✓ Successfully wrote to temp file\n")
  cat("  File size:", file.size(temp_wav), "bytes\n")
}, error = function(e) {
  cat("✗ Error writing to file:", conditionMessage(e), "\n")
  cat("Error details:\n")
  print(e)
})

cat("\nStep 9: Read back the temp file\n")
if (file.exists(temp_wav)) {
  tryCatch({
    readback <- read.AsspDataObj(temp_wav)
    cat("✓ Successfully read back temp file\n")
    cat("  Dimensions match original:", identical(dim(readback$samples), dim(result$samples)), "\n")
  }, error = function(e) {
    cat("✗ Error reading back:", conditionMessage(e), "\n")
  })
  unlink(temp_wav)
}

cat("\nDebug test completed\n")

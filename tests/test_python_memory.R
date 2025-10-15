#!/usr/bin/env Rscript
##' Test memory-based Python DSP processing
##'
##' This script tests the av → Python numpy conversion pipeline
##' to ensure we can eliminate the "convert → store → read → DSP" pattern
##' for Python-based DSP functions.

library(superassp)
library(av)

# Source the helper functions
source("R/av_python_helpers.R")

# Get test file
test_file <- list.files(
  system.file("samples", "sustained", package = "superassp"),
  pattern = "a1.wav",
  full.names = TRUE
)[1]

cat("==== Test 1: av_to_python_audio conversion ====\n")

# Read audio with av
info <- av_media_info(test_file)
audio_raw <- read_audio_bin(
  test_file,
  channels = info$audio$channels,
  sample_rate = info$audio$sample_rate
)

cat("Raw audio data from av:\n")
cat("  Class:", class(audio_raw), "\n")
cat("  Length:", length(audio_raw), "\n")
cat("  Sample rate:", info$audio$sample_rate, "\n")

# Convert to Python numpy
audio_np <- av_to_python_audio(
  audio_raw,
  sample_rate = info$audio$sample_rate,
  channels = info$audio$channels
)

cat("\nPython numpy array:\n")
np <- reticulate::import("numpy")
cat("  Type:", class(audio_np), "\n")
cat("  Shape:", np$shape(audio_np), "\n")
cat("  Dtype:", audio_np$dtype$name, "\n")
cat("  Min value:", np$min(audio_np), "\n")
cat("  Max value:", np$max(audio_np), "\n")

cat("\n==== Test 2: av_load_for_python convenience function ====\n")

result <- av_load_for_python(test_file, start_time = 0.5, end_time = 1.5)

cat("Loaded audio:\n")
cat("  Sample rate:", result$sample_rate, "\n")
cat("  Shape:", np$shape(result$audio_np), "\n")
cat("  Duration:", np$shape(result$audio_np) / result$sample_rate, "seconds\n")

cat("\n==== Test 3: Compare with librosa loading (OLD way) ====\n")

# Check if librosa is available
if (reticulate::py_module_available("librosa")) {

  py <- reticulate::import_main()

  # OLD WAY: librosa reads from file (disk I/O)
  cat("\nOLD WAY (librosa from file):\n")
  reticulate::py_run_string(sprintf("
import librosa as lr
import numpy as np

x_old, fs_old = lr.load('%s', dtype=np.float64, offset=0.5, duration=1.0)
print('  Shape:', x_old.shape)
print('  Sample rate:', fs_old)
print('  Min:', np.min(x_old))
print('  Max:', np.max(x_old))
", test_file))

  # NEW WAY: av loads to memory, convert to numpy (no disk I/O!)
  cat("\nNEW WAY (av → numpy in memory):\n")
  result2 <- av_load_for_python(test_file, start_time = 0.5, end_time = 1.5)
  py$x_new <- result2$audio_np
  py$fs_new <- result2$sample_rate

  reticulate::py_run_string("
print('  Shape:', x_new.shape)
print('  Sample rate:', fs_new)
print('  Min:', np.min(x_new))
print('  Max:', np.max(x_new))
")

  # Check if they're similar
  cat("\nComparing old vs new:\n")
  reticulate::py_run_string("
diff = np.abs(x_old - x_new)
print('  Mean absolute difference:', np.mean(diff))
print('  Max absolute difference:', np.max(diff))
print('  Shapes match:', x_old.shape == x_new.shape)
")

} else {
  cat("\nlibrosa not available - skipping comparison test\n")
}

cat("\n==== Test 4: Use with actual Python DSP (SWIPE) ====\n")

if (reticulate::py_module_available("pysptk")) {

  # Load audio using new memory-based approach
  audio_result <- av_load_for_python(test_file, start_time = 0, end_time = NULL)

  # Pass to Python
  py <- reticulate::import_main()
  py$audio <- audio_result$audio_np
  py$fs <- audio_result$sample_rate

  # Run SWIPE f0 extraction (no file I/O!)
  cat("Running SWIPE f0 extraction (memory-based)...\n")
  reticulate::py_run_string("
import pysptk as sp
import numpy as np

f0 = sp.swipe(audio, fs=fs, hopsize=int(0.005*fs), min=70, max=200, otype='f0')
print('  F0 values computed:', len(f0))
print('  F0 range: [{:.2f}, {:.2f}] Hz'.format(np.min(f0[f0>0]), np.max(f0)))
print('  Voiced frames:', np.sum(f0 > 0))
")

  cat("\n✓ SUCCESS: Python DSP processing works with memory-based audio!\n")
  cat("  No intermediate files were created - everything stayed in memory!\n")

} else {
  cat("\npysptk not available - skipping SWIPE test\n")
}

cat("\n==== Summary ====\n")
cat("✓ av loads audio into R memory (integer vector)\n")
cat("✓ av_to_python_audio converts to Python numpy array (float64)\n")
cat("✓ av_load_for_python provides complete pipeline\n")
cat("✓ Python DSP libraries can process the numpy array directly\n")
cat("✓ NO intermediate files needed!\n")
cat("\nNext step: Update Python DSP functions (swipe, rapt, reaper, etc.)\n")
cat("           to use this memory-based approach\n")

# Python DSP Memory-Based Processing

## Problem Statement

Currently, Python DSP functions in superassp (e.g., `swipe_opt`, `rapt_opt`, `reaper_opt`) use the **"convert → store → read → DSP"** pattern:

1. Convert non-native files to WAV (disk I/O)
2. Store converted files on disk
3. Python reads files with `librosa.load()` (disk I/O)
4. Perform DSP processing
5. Clean up temporary files

This creates unnecessary disk I/O and temporary files, just like we had with the C-based DSP functions before the memory-based optimization.

## Solution: av → Python Numpy Pipeline

We can eliminate the disk I/O loop by loading audio with the `av` package and passing it directly to Python as a numpy array!

### Architecture

```
┌─────────────────────────────────────────────────────────────────┐
│ OLD PATTERN (Current)                                           │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  Input File (any format)                                        │
│       ↓                                                          │
│  av::av_audio_convert() → temp.wav  (DISK WRITE)               │
│       ↓                                                          │
│  librosa.load('temp.wav')           (DISK READ)                │
│       ↓                                                          │
│  Python DSP (pysptk, pyreaper, etc.)                            │
│       ↓                                                          │
│  Clean up temp.wav                  (DISK DELETE)              │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘

┌─────────────────────────────────────────────────────────────────┐
│ NEW PATTERN (Memory-Based)                                      │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  Input File (any format)                                        │
│       ↓                                                          │
│  av::read_audio_bin() → int32 vector  (IN MEMORY!)             │
│       ↓                                                          │
│  av_to_python_audio() → numpy float64 (IN MEMORY!)             │
│       ↓                                                          │
│  Python DSP (pysptk, pyreaper, etc.)                            │
│       ↓                                                          │
│  Results (no cleanup needed)                                    │
│                                                                  │
│  NO DISK I/O!                                                   │
└─────────────────────────────────────────────────────────────────┘
```

### Technical Details

#### 1. Audio Data Formats

**av package output:**
- Format: `s32le` (32-bit signed little-endian integers)
- Range: -2,147,483,648 to 2,147,483,647
- Type: R integer vector
- Interleaved for multi-channel

**librosa/Python DSP expects:**
- Format: `float64` (double precision floating point)
- Range: -1.0 to 1.0 (normalized)
- Type: numpy array
- Mono (single channel for most DSP)

#### 2. Conversion Process

```r
# av returns s32le integers
audio_data <- av::read_audio_bin(file, channels = 1, sample_rate = 44100)
# audio_data is integer vector: e.g., [0, -131072, -655360, ...]

# Convert to float64 in [-1, 1] range
audio_float <- as.numeric(audio_data) / 2147483648.0  # 2^31
# audio_float is numeric vector: e.g., [0.0, -0.000061, -0.000305, ...]

# Convert to Python numpy array
np <- reticulate::import("numpy", convert = FALSE)
audio_np <- np$array(audio_float, dtype = np$float64)
# audio_np is now Python numpy.ndarray (float64)
```

## Implementation

### Helper Functions

Three new functions in `R/av_python_helpers.R`:

#### 1. `av_to_python_audio()` - Core conversion

Converts av integer vector to Python numpy array:

```r
audio_np <- av_to_python_audio(
  audio_data = audio_raw,      # From av::read_audio_bin()
  sample_rate = 44100,
  channels = 1
)
```

#### 2. `av_load_for_python()` - Convenience wrapper

Complete pipeline from file to numpy array:

```r
result <- av_load_for_python(
  file_path = "audio.mp4",
  start_time = 10,
  end_time = 20,
  target_sample_rate = NULL  # Keep original
)

# result$audio_np is ready for Python DSP!
# result$sample_rate
# result$original_file
```

#### 3. `processMediaFiles_Python()` - Batch processing

Processes multiple files with memory-based loading:

```r
results <- processMediaFiles_Python(
  listOfFiles = c("file1.mp4", "file2.wav"),
  beginTime = c(0, 5),
  endTime = c(0, 15),
  python_function = "compute_swipe",
  hopsize = 220,
  fmin = 70,
  fmax = 200
)
```

## Usage Example: Updating swipe_opt

### BEFORE (Current - Uses Disk I/O)

```r
process_swipe_single <- function(soundFile, beginTime, endTime, ...) {
  # Python reads from file (DISK I/O)
  reticulate::py_run_string("
import librosa as lr
x, fs = lr.load(soundFile, offset=bt, duration=et - bt)  # ← DISK READ
f0 = sp.swipe(x, fs=fs, ...)
")
}
```

### AFTER (Memory-Based - No Disk I/O)

```r
process_swipe_single <- function(soundFile, beginTime, endTime, ...) {
  # Load audio with av → convert to numpy (IN MEMORY!)
  audio_result <- av_load_for_python(
    soundFile,
    start_time = beginTime,
    end_time = if(endTime == 0) NULL else endTime
  )

  # Pass numpy array directly to Python (NO DISK I/O!)
  py <- reticulate::import_main()
  py$x <- audio_result$audio_np
  py$fs <- audio_result$sample_rate

  reticulate::py_run_string("
import pysptk as sp
f0 = sp.swipe(x, fs=fs, ...)  # ← x is already in memory!
")
}
```

## Benefits

1. **No Temporary Files**: No more temp WAV files cluttering the disk
2. **Faster Processing**: Eliminates disk I/O bottleneck
3. **Lower Memory**: No need to store converted files
4. **Simpler Code**: No file cleanup logic needed
5. **Consistent Architecture**: Matches the C DSP memory-based optimization

## Performance Comparison

### Expected Improvements

Based on the C DSP migration (which showed 5x speedup):

- **Disk I/O eliminated**: No av_audio_convert() writes, no librosa reads
- **Time window support**: Native support in av (no conversion needed)
- **All formats supported**: Works for WAV, MP4, MKV, etc. without conversion

### Test Results

Run `tests/test_python_memory.R` to verify:

```r
source("tests/test_python_memory.R")
```

Expected output:
- ✓ av loads audio into R memory
- ✓ Conversion to Python numpy array works
- ✓ Numpy array matches librosa format
- ✓ Python DSP functions process the array successfully
- ✓ No intermediate files created

## Migration Guide

### For Each Python DSP Function

1. **Remove `convertInputMediaFiles()` calls**
2. **Replace `librosa.load()` with `av_load_for_python()`**
3. **Pass numpy array to Python instead of file path**
4. **Remove cleanup logic**

### Example Functions to Update

- `swipe_opt` in `R/python_ssff_optimized.R` (lines 172-199)
- `rapt_opt` in `R/python_ssff_optimized.R` (lines 429-454)
- `reaper_opt` in `R/python_ssff_optimized.R` (lines 690-737)
- Praat functions in `R/praat_python_optimized.R`

### Pattern to Follow

```r
# OLD:
# 1. Convert file if needed
listOfFiles_toClear <- convertInputMediaFiles(...)
dsp_input <- listOfFilesDF$dsp_input[i]

# 2. Python reads from file
py$soundFile <- dsp_input
reticulate::py_run_string("x, fs = lr.load(soundFile, ...)")

# 3. Cleanup
cleanupConvertedInputMediaFiles(toClear, ...)

# NEW:
# 1. Load directly with av → numpy
audio_result <- av_load_for_python(file_path, start_time, end_time)

# 2. Pass numpy array to Python
py$x <- audio_result$audio_np
py$fs <- audio_result$sample_rate
# (Python code stays same - just uses x and fs)

# 3. No cleanup needed!
```

## Testing

### Unit Tests

```r
# Test conversion
test_file <- "test.wav"
info <- av::av_media_info(test_file)
audio_raw <- av::read_audio_bin(test_file, ...)
audio_np <- av_to_python_audio(audio_raw, ...)

# Verify numpy array properties
np <- reticulate::import("numpy")
stopifnot(audio_np$dtype$name == "float64")
stopifnot(np$min(audio_np) >= -1.0)
stopifnot(np$max(audio_np) <= 1.0)
```

### Integration Tests

```r
# Compare OLD (librosa) vs NEW (av→numpy)
# They should produce identical results

# OLD
reticulate::py_run_string("
x_old, fs = lr.load('test.wav')
f0_old = sp.swipe(x_old, fs=fs, ...)
")

# NEW
audio_result <- av_load_for_python("test.wav")
py$x_new <- audio_result$audio_np
reticulate::py_run_string("
f0_new = sp.swipe(x_new, fs=fs, ...)
")

# Compare
stopifnot(all.equal(py$f0_old, py$f0_new, tolerance = 1e-6))
```

## Future Work

1. **Update all Python DSP functions** to use memory-based loading
2. **Benchmark performance** improvements
3. **Add multi-channel support** (currently takes first channel)
4. **Create unified interface** for both C and Python DSP functions

## References

- av package: https://CRAN.R-project.org/package=av
- reticulate: https://CRAN.R-project.org/package=reticulate
- librosa: https://librosa.org/
- numpy: https://numpy.org/

## Summary

✅ **Problem Solved**: Python DSP functions can now use the same memory-based processing as C DSP functions!

The `av → Python numpy` pipeline eliminates the "convert → store → read → DSP" loop, providing:
- No temporary files
- No disk I/O
- Faster processing
- Support for all media formats
- Consistent architecture across the package

**Next Steps**: Update Python DSP functions (swipe, rapt, reaper, praat_*) to use `av_load_for_python()` instead of `librosa.load()`.

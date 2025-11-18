# Migration Guide: File-Based to In-Memory Processing

**Date**: 2025-10-30
**Target Version**: v0.9.0+
**Priority**: High

This guide provides step-by-step instructions for migrating file-based DSP functions to in-memory processing using the av package.

---

## Overview

**Goal**: Migrate 12 functions from file-based processing to in-memory processing

**Benefits**:
- ✅ 20-40% faster (no disk I/O overhead)
- ✅ Universal media format support (WAV, MP3, MP4, video, etc.)
- ✅ Better memory management
- ✅ Simpler code (fewer temp file operations)
- ✅ Thread-safe (no file locking issues)

---

## MIGRATION PATTERNS

### Pattern 1: Python Functions Using librosa.load()

**Current Code Pattern**:
```r
# R/ssff_python_example.R (OLD)
trk_example <- function(listOfFiles, ...) {
  # Python code loads file directly
  py_run_string("
import librosa
import numpy as np

def process_audio(file_path, sr=16000):
    # Load audio file from disk
    y, sr = librosa.load(file_path, sr=sr)
    # Process...
    return result
  ")

  result <- py$process_audio(file_path)
  return(result)
}
```

**New Code Pattern**:
```r
# R/ssff_cpp_example.R or R/ssff_python_example.R (NEW)
trk_example <- function(listOfFiles,
                        beginTime = 0.0,
                        endTime = 0.0,
                        ...,
                        toFile = FALSE,
                        verbose = TRUE) {

  # Load audio with av package (in-memory)
  audio_data <- av::read_audio_bin(
    audio = file_path,
    start_time = if (beginTime > 0) beginTime else NULL,
    end_time = if (endTime > 0) endTime else NULL,
    channels = 1
  )

  sample_rate <- attr(audio_data, "sample_rate")

  # Convert to float for Python
  audio_float <- as.numeric(audio_data) / 2147483647.0  # INT32_MAX

  # Pass to Python as numpy array
  np <- reticulate::import("numpy")
  audio_np <- np$array(audio_float, dtype = "float32")

  # Python processes in-memory numpy array
  result <- py$process_audio_array(audio_np, sample_rate)

  # Convert to AsspDataObj
  return(result)
}
```

### Pattern 2: Parselmouth Functions with Temp Files

**Current Code Pattern**:
```r
# R/ssff_python_pm_example.R (OLD)
trk_example_p <- function(listOfFiles, ...) {
  # Create temp file
  temp_file <- tempfile(fileext = ".wav")

  # Convert/copy to temp
  file.copy(listOfFiles, temp_file)

  # Python loads temp file
  py_run_string("
import parselmouth as pm

def process_file(file_path):
    sound = pm.Sound(file_path)  # Loads from disk
    # Process...
    return result
  ")

  result <- py$process_file(temp_file)
  unlink(temp_file)
  return(result)
}
```

**New Code Pattern (Option A - R creates Sound)**:
```r
# R/ssff_python_pm_example.R (NEW - Recommended)
trk_example_p <- function(listOfFiles,
                          beginTime = 0.0,
                          endTime = 0.0,
                          ...) {

  # Load audio and create Parselmouth Sound in R
  sound <- av_load_for_parselmouth(
    file_path = file_path,
    start_time = if (beginTime > 0) beginTime else NULL,
    end_time = if (endTime > 0) endTime else NULL,
    channels = 1
  )

  # Python processes Sound object directly (no file I/O)
  result <- py$process_sound(sound, ...)

  return(result)
}
```

**New Code Pattern (Option B - Python creates Sound from numpy)**:
```r
# R/ssff_python_pm_example.R (NEW - Alternative)
trk_example_p <- function(listOfFiles, ...) {

  # Load audio with av
  audio_data <- av_load_for_python(
    file_path = file_path,
    start_time = bt,
    end_time = if (et > 0) et else NULL
  )

  # Python creates Sound from numpy array (in-memory)
  py_run_string("
import parselmouth as pm
import numpy as np

def process_from_array(audio_np, sample_rate):
    # Create Sound from numpy array (no file)
    sound = pm.Sound(audio_np, sampling_frequency=sample_rate)
    # Process...
    return result
  ")

  result <- py$process_from_array(audio_data$audio_np, audio_data$sample_rate)
  return(result)
}
```

---

## STEP-BY-STEP MIGRATION: trk_yaapt Example

### Step 1: Analyze Current Implementation

```r
# Current: R/ssff_python_yaapt.R
trk_yaapt <- function(listOfFiles, ...) {
  # Uses librosa.load() internally in Python
  py$yaapt_pitch(file_path, ...)
}
```

**Issues**:
- ❌ Python loads file from disk
- ❌ No universal media format support
- ❌ Temp file operations in Python

### Step 2: Load Audio with av Package

```r
trk_yaapt <- function(listOfFiles,
                      beginTime = 0.0,
                      endTime = 0.0,
                      minF = 60.0,
                      maxF = 400.0,
                      windowShift = 10.0,
                      toFile = FALSE,
                      verbose = TRUE) {

  # Validate inputs
  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No input files specified")
  }

  n_files <- length(listOfFiles)
  results <- vector("list", n_files)

  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[i]

    tryCatch({
      # NEW: Load with av package (in-memory)
      audio_data <- av::read_audio_bin(
        audio = file_path,
        start_time = if (bt > 0) bt else NULL,
        end_time = if (et > 0) et else NULL,
        channels = 1
      )

      sample_rate <- attr(audio_data, "sample_rate")

      # Convert to float for Python
      audio_float <- as.numeric(audio_data) / 2147483647.0

      # Pass to Python as numpy array
      np <- reticulate::import("numpy")
      audio_np <- np$array(audio_float, dtype = "float32")

      # MODIFIED: Python processes numpy array (not file)
      result <- py$yaapt_from_array(audio_np, sample_rate, minF, maxF, windowShift)

      # Convert to AsspDataObj
      out_obj <- create_yaapt_asspobj(result, windowShift)

      # Handle output
      if (toFile) {
        # Write to file
        results[[i]] <- write_output(out_obj, file_path)
      } else {
        results[[i]] <- out_obj
      }

    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(file_path)}}: {e$message}")
      results[[i]] <- NULL
    })
  }

  return(if (n_files == 1) results[[1]] else results)
}
```

### Step 3: Update Python Code

```python
# inst/python/yaapt_wrapper.py (OLD)
def yaapt_pitch(file_path, minF, maxF, frame_shift):
    import librosa
    y, sr = librosa.load(file_path, sr=16000)  # Loads from disk
    # Process...
    return result

# inst/python/yaapt_wrapper.py (NEW)
def yaapt_from_array(audio_np, sample_rate, minF, maxF, frame_shift):
    """
    Process audio from numpy array (no file I/O).

    Parameters:
    -----------
    audio_np : np.ndarray
        Audio samples as float32 array
    sample_rate : int
        Sampling rate in Hz
    """
    # Process in-memory numpy array
    # No file loading needed!
    result = yaapt.process(audio_np, sr=sample_rate, ...)
    return result
```

### Step 4: Add Tests

```r
# tests/testthat/test-yaapt-memory.R
test_that("trk_yaapt works with in-memory processing", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")

  # Test in-memory (toFile = FALSE)
  result <- trk_yaapt(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("F0" %in% names(result))
  expect_true(nrow(result$F0) > 0)
})

test_that("trk_yaapt handles MP3 files", {
  test_mp3 <- system.file("samples", "sustained", "a7.mp3", package = "superassp")
  skip_if(test_mp3 == "", "MP3 not available")

  # Should work with any av-supported format
  result <- trk_yaapt(test_mp3, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
})
```

### Step 5: Verify Performance

```r
# Benchmark old vs new
library(microbenchmark)

old_time <- system.time(trk_yaapt_old("audio.wav"))
new_time <- system.time(trk_yaapt("audio.wav"))

cat("Old:", old_time[3], "s\n")
cat("New:", new_time[3], "s\n")
cat("Speedup:", old_time[3] / new_time[3], "x\n")

# Expected: 1.2-1.4x speedup (20-40% faster)
```

---

## MIGRATION CHECKLIST

For each function being migrated:

### Planning Phase
- [ ] Identify current file I/O operations
- [ ] Check Python dependencies (librosa, parselmouth, etc.)
- [ ] Review existing tests
- [ ] Estimate complexity (simple/medium/complex)

### Implementation Phase
- [ ] Add `av::read_audio_bin()` or `av_load_for_parselmouth()` call
- [ ] Update Python code to accept numpy arrays (not file paths)
- [ ] Remove temp file operations
- [ ] Add time windowing support (beginTime/endTime)
- [ ] Update function documentation

### Testing Phase
- [ ] Run existing tests (ensure backward compatibility)
- [ ] Add new tests for media formats (MP3, MP4, FLAC, OGG)
- [ ] Test time windowing functionality
- [ ] Benchmark performance (should be 20-40% faster)
- [ ] Test memory usage (should be <2x increase)

### Documentation Phase
- [ ] Update roxygen2 documentation
- [ ] Add examples with different media formats
- [ ] Document performance improvements
- [ ] Update NEWS.md

### Validation Phase
- [ ] Run R CMD check (should pass)
- [ ] Run full test suite (should pass)
- [ ] Verify no temp file creation (check /tmp during execution)
- [ ] Test with large files (ensure memory efficiency)

---

## PRIORITY LIST (12 Functions to Migrate)

### High Priority (5 functions)
1. **trk_yaapt** - Popular pitch tracker
2. **trk_snackp** - Snack pitch tracking
3. **trk_snackf** - Snack formant tracking
4. **trk_formantp** - Praat formant tracking (high usage)
5. **trk_formantpathp** - Praat formant path tracking

### Medium Priority (5 functions)
6. **trk_intensityp** - Praat intensity
7. **trk_praat_sauce** - Voice quality (specialized)
8. **trk_spectral_momentsp** - Spectral moments
9. **trk_excite** - Source-filter decomposition
10. **trk_seenc** - Spectral envelope

### Low Priority (2 functions)
11. **trk_aperiodicities** - Aperiodicity (consider deprecation)
12. **trk_straight_synth** - STRAIGHT synthesis (consider deprecation)

---

## HELPER FUNCTIONS AVAILABLE

The package provides several helper functions for migration:

### av_to_asspDataObj()
```r
# Load any media format to AsspDataObj
audio_obj <- av_to_asspDataObj(
  file_path = "audio.mp3",
  start_time = 1.0,
  end_time = 3.0
)
```

### av_load_for_python()
```r
# Load audio for Python processing
audio_data <- av_load_for_python(
  file_path = "audio.mp4",
  start_time = 0,
  end_time = NULL,
  channels = 1,
  target_sample_rate = 16000
)

# Returns list with:
# - audio_np: numpy array
# - sample_rate: int
```

### av_load_for_parselmouth()
```r
# Load audio and create Parselmouth Sound object
sound <- av_load_for_parselmouth(
  file_path = "audio.wav",
  start_time = 1.0,
  end_time = 3.0,
  channels = 1
)

# Returns parselmouth.Sound object (no file I/O)
```

---

## COMMON PITFALLS & SOLUTIONS

### Pitfall 1: Forgetting INT32 to Float Conversion

```r
# WRONG (causes clipping)
audio_np <- np$array(as.numeric(audio_data))

# CORRECT
audio_np <- np$array(as.numeric(audio_data) / 2147483647.0)
```

### Pitfall 2: Not Handling NULL endTime

```r
# WRONG (av throws error)
audio_data <- av::read_audio_bin(audio, end_time = 0.0)

# CORRECT
audio_data <- av::read_audio_bin(
  audio,
  end_time = if (endTime > 0) endTime else NULL
)
```

### Pitfall 3: Hardcoding Sample Rate

```r
# WRONG (doesn't respect original sample rate)
audio_data <- av::read_audio_bin(audio, sample_rate = 16000)

# CORRECT (use original rate or make it a parameter)
audio_data <- av::read_audio_bin(audio, sample_rate = target_sr)
```

### Pitfall 4: Not Cleaning Up Temp Files in Tests

```r
# WRONG (leaves temp files)
test_that("function works", {
  temp_file <- tempfile()
  result <- my_function(test_wav, toFile = TRUE, outputDirectory = dirname(temp_file))
  # Test assertions...
  # Temp file not deleted!
})

# CORRECT
test_that("function works", {
  temp_dir <- tempdir()
  result <- my_function(test_wav, toFile = TRUE, outputDirectory = temp_dir)
  # Test assertions...
  unlink(result)  # Clean up
})
```

---

## PERFORMANCE EXPECTATIONS

Based on migrations completed:

| Function | Old (File) | New (Memory) | Speedup |
|----------|-----------|--------------|---------|
| trk_yin | ~110ms | ~35ms | 3.1x |
| trk_dysprosody | ~440ms | ~270ms | 1.6x |
| Average | N/A | N/A | **1.2-1.4x** |

**Typical improvements**:
- Simple functions: 20-30% faster
- I/O-heavy functions: 30-50% faster
- Multi-pass functions: 40-60% faster

---

## TESTING TEMPLATE

```r
# tests/testthat/test-function-memory.R

test_that("function works with in-memory processing", {
  skip_if_not_installed("superassp")
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- my_function(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true(length(names(result)) > 0)
})

test_that("function handles multiple media formats", {
  formats <- c("a1.wav", "a7.mp3", "a2.flac", "a8.ogg")

  for (fmt in formats) {
    test_file <- system.file("samples", "sustained", fmt, package = "superassp")
    if (test_file != "" && file.exists(test_file)) {
      result <- my_function(test_file, toFile = FALSE, verbose = FALSE)
      expect_s3_class(result, "AsspDataObj")
    }
  }
})

test_that("function handles time windowing", {
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")

  result <- my_function(test_wav, beginTime = 0.5, endTime = 1.5,
                       toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  duration <- attr(result, "endRecord") / attr(result, "sampleRate")
  expect_true(duration <= 1.1)  # ~1 second window
})

test_that("function performance improved", {
  skip("Manual performance check only")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")

  # Time 10 runs
  times <- replicate(10, {
    system.time(my_function(test_wav, toFile = FALSE, verbose = FALSE))[3]
  })

  median_time <- median(times)
  cat("Median time:", median_time, "seconds\n")

  # Should be faster than file-based version
  # (compare with old implementation if available)
})
```

---

## CONTACT & SUPPORT

- **Questions**: Create GitHub issue
- **Migration Help**: fredrik.nylen@umu.se
- **Testing**: Run full test suite before submitting PR

---

**Last Updated**: 2025-10-30
**Next Review**: After v0.9.0 migration sprint

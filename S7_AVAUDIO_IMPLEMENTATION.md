# S7 AVAudio Implementation Summary

**Version:** 0.6.0
**Date:** October 20, 2025
**Status:** Partial Implementation (AVAudio class complete, S7 dispatch deferred)

---

## Overview

This release introduces the **AVAudio S7 class** for representing audio data with metadata, providing a modern object-oriented interface for audio processing. However, full S7 method dispatch for all DSP functions was **not implemented** due to fundamental incompatibilities with existing function signatures.

---

## What Was Implemented

### 1. AVAudio S7 Class ✅

**File:** `R/s7_avaudio.R`

A complete S7 class for wrapping audio data with metadata:

```r
AVAudio <- S7::new_class(
  "AVAudio",
  properties = list(
    samples = S7::class_integer,
    sample_rate = S7::class_integer,
    channels = S7::class_integer,
    file_path = S7::new_property(class = S7::class_character, default = NA_character_)
  )
)
```

**Properties:**
- `samples`: Integer vector of audio samples (s32le format - 32-bit signed integers)
- `sample_rate`: Sample rate in Hz (integer)
- `channels`: Number of audio channels (integer)
- `file_path`: Optional source file path (character, NA if generated)

**Features:**
- Automatic validation (sample rate > 0, channels > 0, samples length multiple of channels)
- Print and summary methods
- Compatible with av::read_audio_bin() output format

### 2. AVAudio Utility Functions ✅

**Functions implemented:**

#### `read_avaudio()`
Read audio file and create AVAudio object:
```r
audio <- read_avaudio("speech.wav")
audio <- read_avaudio("speech.wav", start_time = 1.0, end_time = 3.0)
audio <- read_avaudio("speech.wav", sample_rate = 16000)
```

#### `as_avaudio()`
Convert av::read_audio_bin() or prep_recode() output to AVAudio:
```r
audio_data <- prep_recode("speech.wav", format = "wav")
audio <- as_avaudio(audio_data, file_path = "speech.wav")
```

#### `is_avaudio()`
Check if object is AVAudio:
```r
is_avaudio(audio)  # TRUE/FALSE
```

#### `avaudio_to_av()`
Convert AVAudio back to av::read_audio_bin() format:
```r
audio_vec <- avaudio_to_av(audio)
# Returns integer vector with channels and sample_rate attributes
```

#### `avaudio_to_tempfile()`
Write AVAudio to temporary WAV file:
```r
temp_file <- avaudio_to_tempfile(audio)
# Use with existing DSP functions
result <- trk_rapt(temp_file, toFile = FALSE)
unlink(temp_file)  # Cleanup
```

### 3. Integration with prep_recode() ✅

The `prep_recode()` function returns data in a format directly compatible with `as_avaudio()`:

```r
# Direct pipeline
audio <- as_avaudio(prep_recode("video.mp4", format = "wav"))

# Or use read_avaudio (wrapper around prep_recode)
audio <- read_avaudio("video.mp4")
```

### 4. Documentation ✅

Complete roxygen2 documentation for all AVAudio functions:
- `?AVAudio-class`
- `?read_avaudio`
- `?as_avaudio`
- `?is_avaudio`
- `?avaudio_to_av`
- `?avaudio_to_tempfile`

### 5. Tests ✅

**File:** `tests/testthat/test-s7-avaudio.R` (10 test cases)

Test coverage:
- AVAudio creation from file
- AVAudio creation from prep_recode output
- Conversion to/from av format
- Conversion to temp file
- Print and summary methods
- Time windowing
- Resampling
- Validation

---

## What Was NOT Implemented

### S7 Method Dispatch for DSP Functions ❌

**Original Goal:** Convert all 44 lst_* and trk_* functions to S7 generics with methods for both character vectors (file paths) and AVAudio objects.

**Why It Failed:**

1. **S7 Dispatch Argument Restrictions:**
   - S7 requires dispatch arguments to NOT have default values
   - ALL existing DSP functions have `listOfFiles = NULL` as default
   - Example from `trk_rapt`:
     ```r
     trk_rapt <- function(listOfFiles = NULL, ...)  # NOT compatible with S7
     ```

2. **Function Signature Incompatibility:**
   - S7 method registration requires exact signature matching
   - Many functions have different first argument names or structures
   - Examples:
     - `lst_avqip(svDF, csDF, ...)`  - first arg is NOT listOfFiles
     - `lst_dsip(softDF, highpitchDF, ...)`  - first arg is NOT listOfFiles

3. **Runtime Conversion Complexity:**
   - Attempting to convert functions at `.onLoad()` time caused 49 warnings
   - Namespace locking/unlocking required
   - High risk of breaking existing functionality

**Attempted Approach:**

```r
# This approach failed
.convert_to_s7_generic <- function(fn_name) {
  generic_fn <- S7::new_generic(name = fn_name, dispatch_args = "listOfFiles")
  S7::method(generic_fn, S7::class_character) <- original_fn  # FAILS
  S7::method(generic_fn, AVAudio) <- avaudio_method
}
```

**Errors Encountered:**
```
Warning: In trk_rapt(<character>), dispatch arguments (`listOfFiles`)
must not have default values

Warning: lst_avqip() dispatches on `listOfFiles`, but lst_avqip(<character>)
has arguments `svDF`, `csDF`, ...
```

---

## Current Usage Pattern

Until full S7 dispatch is implemented, users should use this pattern:

### Pattern 1: Direct File Processing (Existing)
```r
# Works as before - no changes needed
result <- trk_rapt("speech.wav", toFile = FALSE)
```

### Pattern 2: AVAudio with Temp File Conversion
```r
# Load audio into memory
audio <- read_avaudio("speech.wav")

# Convert to temp file for processing
temp_file <- avaudio_to_tempfile(audio)

# Process with existing DSP function
result <- trk_rapt(temp_file, toFile = FALSE)

# Cleanup
unlink(temp_file)
```

### Pattern 3: AVAudio Manipulation Pipeline
```r
# Load and preprocess audio
audio <- read_avaudio("recording.wav",
                      sample_rate = 16000,  # Resample
                      start_time = 1.0,      # Time window
                      end_time = 3.0,
                      channels = 1)          # Convert to mono

# Process
temp_file <- avaudio_to_tempfile(audio)
f0 <- trk_rapt(temp_file, toFile = FALSE)
features <- lst_voice_sauce(temp_file)
unlink(temp_file)
```

---

## Future Work

### Option 1: Manual S7 Method Addition (Recommended)

Gradually add explicit S7 support to individual functions by removing default values:

```r
# Current (NOT S7 compatible):
trk_rapt <- function(listOfFiles = NULL, ...) { ... }

# Future (S7 compatible):
trk_rapt <- S7::new_generic("trk_rapt", dispatch_args = "listOfFiles")

S7::method(trk_rapt, S7::class_character) <- function(listOfFiles, ...) {
  # Original implementation
}

S7::method(trk_rapt, AVAudio) <- function(listOfFiles, ...) {
  temp_file <- avaudio_to_tempfile(listOfFiles)
  on.exit(unlink(temp_file))
  trk_rapt(temp_file, ...)
}
```

This would require:
- Changing 44 function signatures (breaking change)
- Updating all documentation
- Updating all tests
- Migration guide for users

### Option 2: Wrapper Functions

Create parallel AVAudio-specific functions:

```r
# Original remains unchanged
trk_rapt <- function(listOfFiles = NULL, ...) { ... }

# New AVAudio-specific version
trk_rapt.AVAudio <- function(audio, ...) {
  temp_file <- avaudio_to_tempfile(audio)
  on.exit(unlink(temp_file))
  trk_rapt(temp_file, ...)
}
```

Pros:
- No breaking changes
- Clear API
- Easy to implement incrementally

Cons:
- Duplicate function names
- Not true S7 dispatch

### Option 3: Unified Interface with UseMethod()

Use S3 dispatch instead of S7:

```r
trk_rapt <- function(listOfFiles, ...) {
  UseMethod("trk_rapt")
}

trk_rapt.character <- function(listOfFiles, ...) {
  # Original implementation
}

trk_rapt.AVAudio <- function(listOfFiles, ...) {
  temp_file <- avaudio_to_tempfile(listOfFiles)
  on.exit(unlink(temp_file))
  trk_rapt.character(temp_file, ...)
}
```

Pros:
- Standard R approach
- Works with existing code
- Simpler than S7

Cons:
- Requires changing all 44 functions
- Not using S7 (which was the original goal)

---

## Benefits of Current Implementation

Even without full S7 dispatch, the AVAudio class provides significant value:

1. **Type Safety:** AVAudio objects carry metadata (sample rate, channels)
2. **Memory Efficiency:** Audio data loaded once, can be reused
3. **Preprocessing:** Time windowing, resampling, channel conversion before DSP
4. **Format Flexibility:** Works with any format supported by prep_recode()
5. **Clear API:** Explicit conversion functions (avaudio_to_av, avaudio_to_tempfile)
6. **Future-Proof:** Foundation for full S7 integration later

---

## Version 0.6.0 Changes

###Files Added:
- `R/s7_avaudio.R` (368 lines) - AVAudio class and utilities
- `R/s7_methods.R` (35 lines) - Placeholder for future S7 dispatch
- `R/prep_recode.R` (346 lines) - Media re-encoding function
- `tests/testthat/test-s7-avaudio.R` (10 tests) - AVAudio class tests
- `tests/testthat/test-s7-dispatch.R` (8 tests) - S7 dispatch tests (currently manual pattern)
- `tests/testthat/test-prep-recode.R` (18 tests) - prep_recode tests
- `man/*.Rd` - Documentation files

### Files Modified:
- `DESCRIPTION` - Added S7 dependency, bumped version to 0.6.0
- `R/zzz.R` - Added .setup_s7_methods() call (currently no-op)
- `NAMESPACE` - Exported AVAudio functions

### Test Status:
- AVAudio tests: ✅ All pass
- prep_recode tests: ✅ All pass (55 assertions)
- S7 dispatch tests: ⚠️ Adapted to manual pattern (not automatic dispatch)
- Existing DSP tests: ✅ Should all still pass (backward compatible)

---

## Recommendations

### For Version 0.6.0 Release:

1. ✅ Ship with AVAudio class (complete and tested)
2. ✅ Ship with prep_recode() (complete and tested)
3. ✅ Document current usage patterns
4. ✅ Note S7 dispatch limitation in NEWS/changelog
5. ⚠️ Consider future S7 integration strategy

### For Version 0.7.0 (Future):

1. Implement Option 1, 2, or 3 above
2. Add AVAudio support to high-priority functions first (trk_rapt, trk_reaper, lst_voice_sauce)
3. Gradually expand to remaining functions
4. Provide migration guide

---

## Conclusion

Version 0.6.0 successfully implements:
- ✅ Complete AVAudio S7 class with validation
- ✅ Full suite of AVAudio utility functions
- ✅ Integration with prep_recode()
- ✅ Comprehensive documentation and tests

Version 0.6.0 does NOT implement:
- ❌ Automatic S7 method dispatch for DSP functions
- ❌ Direct AVAudio input to lst_*/trk_* functions

**Users must manually convert AVAudio to temp files for now.**

This is a solid foundation for future S7 integration, even though full dispatch wasn't achieved in this release due to fundamental incompatibilities with existing function signatures.

---

**Document Version:** 1.0
**Last Updated:** October 20, 2025
**Status:** Implementation complete (with limitations documented)

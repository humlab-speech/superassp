# Parselmouth/Praat Memory-Based Processing Investigation

## Current Architecture

Praat/Parselmouth functions currently load audio from disk:

```python
# In inst/python/praat_pitch.py, line 97:
snd = pm.Sound(soundFile)
```

This means:
1. R passes file path to Python
2. Python/Parselmouth reads file from disk
3. Parselmouth creates Sound object
4. DSP processing occurs

## Goal: Memory-Based Processing

Eliminate disk I/O by creating Parselmouth Sound objects directly from av audio data (similar to what we did with numpy arrays for pysptk/pyreaper).

## Research: Parselmouth Sound Object Creation

### Option 1: Sound from numpy array (MOST PROMISING)

Parselmouth provides `Sound()` constructor that can accept numpy arrays directly!

```python
import parselmouth as pm
import numpy as np

# Create Sound from numpy array
# Parameters: values (numpy array), sampling_frequency (Hz)
sound = pm.Sound(values=audio_np, sampling_frequency=sample_rate)
```

**Documentation**: https://parselmouth.readthedocs.io/en/stable/api_reference.html#parselmouth.Sound

This is **exactly** what we need! We can:
1. Load audio with av → int32 vector (R memory)
2. Convert to numpy array float64 (via av_to_python_audio())
3. Pass to Python and create Sound object from numpy array
4. No disk I/O!

### Option 2: Sound.from_wav_file() (what's currently used)

This is what `pm.Sound(soundFile)` calls internally - it reads from disk.

**NOT suitable** for memory-based processing.

### Option 3: Sound from raw audio data

Parselmouth also supports creating Sound from raw audio bytes, but numpy array approach is simpler and more direct.

## Implementation Strategy

### Architecture Comparison

**OLD Pattern** (current):
```
Input File → av::av_audio_convert() → temp.wav (DISK WRITE)
          → pm.Sound(temp.wav)              (DISK READ)
          → Parselmouth DSP
          → cleanup temp.wav                (DISK DELETE)
```

**NEW Pattern** (memory-based):
```
Input File → av::read_audio_bin() → int32 vector  (R MEMORY)
          → av_to_python_audio()  → numpy array   (PYTHON MEMORY)
          → pm.Sound(audio_np, fs) → Sound object (PYTHON MEMORY)
          → Parselmouth DSP
          → NO cleanup needed
```

### Code Changes Required

#### R Side (praat_python_optimized.R functions)

BEFORE:
```r
# Call Python function passing file path
result_df <- reticulate::py$praat_pitch(
  origSoundFile,    # ← File path passed to Python
  beginTime = bt,
  endTime = et,
  ...
)
```

AFTER:
```r
# Load audio with av → convert to numpy
audio_result <- av_load_for_python(
  origSoundFile,
  start_time = bt,
  end_time = if(et == 0) NULL else et
)

# Pass numpy array to Python
py <- reticulate::import_main()
py$audio <- audio_result$audio_np
py$fs <- audio_result$sample_rate

# Call Python function with numpy array
result_df <- reticulate::py$praat_pitch_memory(
  audio_np = py$audio,    # ← Numpy array instead of file path
  sample_rate = py$fs,
  # Note: beginTime/endTime no longer needed (already extracted)
  ...
)
```

#### Python Side (inst/python/praat_*.py files)

BEFORE:
```python
def praat_pitch(soundFile, beginTime=0.0, endTime=0.0, ...):
    # Load sound from file
    snd = pm.Sound(soundFile)
    dur = snd.get_total_duration()

    # Handle time windowing
    if beginTime > 0.0 or endTime > 0.0:
        if beginTime >= 0.0 and endTime <= dur:
            snd = snd.extract_part(beginTime, endTime, windowShape, relativeWidth, True)
```

AFTER:
```python
def praat_pitch_memory(audio_np, sample_rate, ...):
    # Create Sound from numpy array (IN MEMORY!)
    snd = pm.Sound(values=audio_np, sampling_frequency=sample_rate)

    # No time windowing needed - already done by av!
```

### Benefits

1. **No disk I/O**: Audio stays in memory from av → numpy → Sound
2. **Time windowing**: Handled natively by av (no need for extract_part)
3. **Format flexibility**: Works with all formats av supports
4. **Consistent architecture**: Same pattern as pysptk/pyreaper functions
5. **Performance**: Expected 5-15x speedup (same as other memory-based functions)

## Implementation Plan

### Phase 1: Create helper function

Add to `R/av_python_helpers.R`:
```r
#' Create Parselmouth Sound object from av audio data
#'
#' @param audio_data Integer vector from av::read_audio_bin
#' @param sample_rate Sample rate in Hz
#' @param channels Number of channels
#' @return Python parselmouth.Sound object
#' @keywords internal
av_to_parselmouth_sound <- function(audio_data, sample_rate, channels = 1) {
  # Convert to numpy array
  audio_np <- av_to_python_audio(audio_data, sample_rate, channels)

  # Create Parselmouth Sound object
  pm <- reticulate::import("parselmouth")
  sound <- pm$Sound(values = audio_np, sampling_frequency = sample_rate)

  return(sound)
}
```

### Phase 2: Update Python scripts

Create memory-based versions of Python scripts:
- `inst/python/praat_pitch_memory.py`
- `inst/python/praat_formant_burg_memory.py`
- `inst/python/praat_intensity_memory.py`
- `inst/python/praat_spectral_moments_memory.py`
- `inst/python/praat_formantpath_burg_memory.py`

Each should:
1. Accept `audio_np` and `sample_rate` instead of `soundFile`
2. Create Sound with `pm.Sound(values=audio_np, sampling_frequency=sample_rate)`
3. Remove time windowing logic (av handles it)

### Phase 3: Update R wrapper functions

Update functions in `R/praat_python_optimized.R`:
- `praat_pitch_opt`
- `praat_formant_burg_opt`
- `praat_intensity_opt`
- `praat_spectral_moments_opt`
- `praat_formantpath_burg_opt`

Each should:
1. Use `av_load_for_python()` to get numpy array
2. Pass numpy array and sample rate to Python
3. Call memory-based Python function
4. Remove file conversion logic

## Testing

### Verification Tests

1. **Sound object creation**: Verify pm.Sound(values=...) works correctly
2. **Equivalence test**: Compare results from file-based vs memory-based
3. **Time windowing**: Verify av time extraction matches extract_part()
4. **All formats**: Test with WAV, MP4, FLAC, etc.

### Benchmark Tests

Expected performance improvements:
- Small files: 2-3x faster
- Medium files: 5-8x faster
- Large files: 10-15x faster
- Time windowing: 15-30x faster

## Next Steps

1. ✓ Document Parselmouth Sound object creation methods
2. ✓ Confirm pm.Sound(values=numpy_array) is supported
3. → Create av_to_parselmouth_sound() helper function
4. → Update ONE Python script as proof of concept (start with praat_pitch)
5. → Update corresponding R function
6. → Test and verify correctness
7. → Benchmark performance improvements
8. → Update remaining Python/R functions
9. → Update documentation

## References

- Parselmouth documentation: https://parselmouth.readthedocs.io/
- Parselmouth Sound API: https://parselmouth.readthedocs.io/en/stable/api_reference.html#parselmouth.Sound
- av package: https://CRAN.R-project.org/package=av

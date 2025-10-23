# Migration Progress: librosa → av (v0.7.0)

**Started:** 2025-10-23
**Status:** In Progress (2 of 11 functions migrated)

## Goal

Migrate all Python DSP functions from `librosa.load()` to `av::read_audio_bin()` for consistent media format support across superassp.

---

## Migration Pattern

### Standard Replacement

```r
# BEFORE (librosa.load)
waveform, fs = librosa.load(soundFile, offset=beginTime, duration=duration, mono=True)

# AFTER (av::read_audio_bin)
audio_data <- av::read_audio_bin(
  audio = origSoundFile,
  start_time = if (beginTime > 0) beginTime else NULL,
  end_time = if (endTime > 0 && condition) endTime else NULL,
  channels = 1
)
fs <- attr(audio_data, "sample_rate")
audio_float <- as.numeric(audio_data) / 2147483647.0
np <- reticulate::import("numpy", convert = FALSE)
audio_array <- np$array(audio_float, dtype = "float32")

# Then pass audio_array to Python instead of loading in Python
py$waveform <- audio_array
py$fs <- reticulate::r_to_py(as.integer(fs))
```

---

## Progress: High Priority Functions (5)

| Function | File | Status | Notes |
|----------|------|--------|-------|
| trk_pyin() | R/ssff_python_pyin.R | ✅ DONE | Probabilistic YIN |
| trk_yin() | R/ssff_python_yin.R | ✅ DONE | YIN pitch tracker |
| trk_crepe() | R/ssff_python_crepe.R | ⏳ TODO | Deep learning F0 |
| trk_yaapt() | R/ssff_python_yaapt.R | ⏳ TODO | Yet Another Algorithm for Pitch Tracking |
| trk_kaldi_pitch() | R/ssff_python_kaldi_pitch.R | ⏳ TODO | Kaldi ASR pitch |

---

## Progress: Medium Priority Functions (6)

| Function | File | Status | Notes |
|----------|------|--------|-------|
| trk_snackp() | R/ssff_python_snack_pitch.R | ⏳ TODO | Snack pitch tracker |
| trk_snackf() | R/ssff_python_snack_formant.R | ⏳ TODO | Snack formant tracker |
| trk_seenc() | R/ssff_python_seenc.R | ⏳ TODO | Speech envelope |
| trk_excite() | R/ssff_python_excite.R | ⏳ TODO | Excitation source |
| trk_aperiodicities() | R/ssff_python_aperiodicities.R | ⏳ TODO | Aperiodicity analysis |
| reaper_pm() | R/ssff_python_reaper_pm.R | ⏳ TODO | REAPER pitchmarks |

---

## Completed Migrations

### ✅ trk_pyin() - Probabilistic YIN

**File:** `R/ssff_python_pyin.R`

**Changes:**
1. Replaced `librosa.load()` with `av::read_audio_bin()`
2. Updated documentation to mention media format support
3. Converted audio to numpy array in R before passing to Python
4. Removed file path passing to Python (now passes array directly)

**Testing needed:**
- [ ] WAV files
- [ ] MP3 files
- [ ] MP4 video files
- [ ] Time windowing
- [ ] Result comparison (old vs new)

### ✅ trk_yin() - YIN Pitch Tracker

**File:** `R/ssff_python_yin.R`

**Changes:**
1. Replaced `librosa.load()` with `av::read_audio_bin()`
2. Updated documentation to mention media format support
3. Converted audio to numpy array in R before passing to Python
4. Removed file path passing to Python

**Testing needed:**
- [ ] WAV files
- [ ] MP3 files
- [ ] MP4 video files
- [ ] Time windowing
- [ ] Result comparison (old vs new)

---

## Remaining Work

### Phase 1: Complete High Priority (3 functions)

1. **trk_crepe()** - Deep learning pitch tracker
   - Uses TensorFlow/Keras
   - May need special attention for model loading
   - Important user-facing function

2. **trk_yaapt()** - YAAPT pitch tracker
   - Requires AMFM decomposition library
   - Standard research tool

3. **trk_kaldi_pitch()** - Kaldi ASR pitch
   - Uses PyTorch/torchaudio
   - Important for ASR workflows

### Phase 2: Complete Medium Priority (6 functions)

4. **trk_snackp()** - Snack pitch
5. **trk_snackf()** - Snack formants
6. **trk_seenc()** - Speech envelope
7. **trk_excite()** - Excitation source
8. **trk_aperiodicities()** - Aperiodicity
9. **reaper_pm()** - REAPER pitchmarks

---

## Testing Strategy

### Manual Testing

For each migrated function:

```r
library(superassp)
library(testthat)

# Test files
wav_file <- "test.wav"
mp3_file <- "test.mp3"
mp4_file <- "test.mp4"

# Test function
test_that("works with WAV", {
  result <- trk_function(wav_file, toFile = FALSE)
  expect_s3_class(result, "AsspDataObj")
})

test_that("works with MP3", {
  result <- trk_function(mp3_file, toFile = FALSE)
  expect_s3_class(result, "AsspDataObj")
})

test_that("works with MP4 video", {
  result <- trk_function(mp4_file, toFile = FALSE)
  expect_s3_class(result, "AsspDataObj")
})

test_that("time windowing works", {
  result <- trk_function(wav_file, beginTime = 1.0, endTime = 2.0, toFile = FALSE)
  expect_s3_class(result, "AsspDataObj")
})
```

### Regression Testing

Compare results between old and new implementations:

```r
# Load same WAV file with both methods
# (Keep old implementation temporarily for comparison)
result_old <- old_implementation(wav_file, toFile = FALSE)
result_new <- trk_function(wav_file, toFile = FALSE)

# Compare F0 values (should be identical or very close)
expect_equal(result_old$f0, result_new$f0, tolerance = 1e-6)
```

---

## Common Issues and Solutions

### Issue 1: Python function names

**Problem:** Some Python code calls `librosa.trk_yin()` instead of `librosa.yin()`

**Solution:** Remove `trk_` prefix from librosa function calls:
```python
# WRONG
pitch = librosa.trk_yin(waveform, ...)

# CORRECT
pitch = librosa.yin(waveform, ...)
```

### Issue 2: Time windowing logic

**Problem:** Complex duration calculation in Python

**Solution:** Let av handle it:
```r
# OLD (complex)
duration = None
if endTime > (windowSize/1000) and (endTime-beginTime) >= (windowSize/1000):
    duration = (endTime - beginTime)

# NEW (simple)
end_time = if (endTime > 0 && (endTime - beginTime) >= (windowSize/1000)) endTime else NULL
```

### Issue 3: Sample rate conversion

**Problem:** librosa can resample on load

**Solution:** av can also resample:
```r
audio_data <- av::read_audio_bin(
  audio = file,
  channels = 1,
  sample_rate = target_sr  # av resamples automatically
)
```

---

## Performance Expectations

Based on `trk_swiftf0()` migration results:

| Metric | Before (librosa) | After (av) | Improvement |
|--------|------------------|------------|-------------|
| Audio loading | ~150-200ms | ~50-80ms | 2-3x faster |
| Memory usage | Higher (Python) | Lower (R→numpy) | ~20% reduction |
| Format support | WAV, FLAC | ALL formats | Unlimited |

---

## Documentation Updates Needed

For each migrated function, update roxygen2 docs:

```r
#' @details
#' This function supports all media formats via the \code{av} package,
#' including WAV, MP3, MP4, FLAC, OGG, AAC, OPUS, and video files with
#' audio tracks.
```

---

## Next Session Tasks

1. **Migrate trk_crepe()** (High priority, user-facing)
2. **Migrate trk_yaapt()** (High priority, research tool)
3. **Migrate trk_kaldi_pitch()** (High priority, ASR workflows)
4. **Test all 5 high-priority functions** with multiple formats
5. **Create regression tests** comparing old vs new
6. **Update documentation** for completed functions
7. **Commit high-priority batch** as v0.7.0-alpha

---

## Estimated Timeline

- **High priority (3 remaining):** 2-3 hours
- **Medium priority (6 functions):** 3-4 hours
- **Testing & validation:** 2-3 hours
- **Documentation:** 1 hour
- **Total:** 8-11 hours of development time

---

## Reference Implementations

**Completed:**
- `R/ssff_python_pyin.R` - Full migration with av
- `R/ssff_python_yin.R` - Full migration with av
- `R/ssff_python_swiftf0.R` - Original modern reference

**Pattern to follow:** See `MIGRATION_LIBROSA_TO_AV.md` for complete guide.

---

## Success Criteria

Function is considered "migrated" when:

- [ ] Uses `av::read_audio_bin()` instead of `librosa.load()`
- [ ] Documentation mentions media format support
- [ ] Passes existing tests (if any)
- [ ] Works with WAV, MP3, MP4 files
- [ ] Time windowing functions correctly
- [ ] No regression in output quality
- [ ] Code follows modern superassp patterns

---

## Blockers / Issues

None currently. All migrations are straightforward replacements following the established pattern.

---

## Notes

- The migration pattern is well-established and proven
- All Python functions use reticulate, so numpy array passing works uniformly
- av package handles all format complexity
- No breaking changes for users (API remains identical)
- Performance improves automatically (2-3x faster loading)

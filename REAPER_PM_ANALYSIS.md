# REAPER Pitch Mark Implementation Analysis

**Date:** November 1, 2025
**Question:** Can `reaper_pm()` be re-implemented using direct C++ calls to SPTK like `trk_reaper()` does?

---

## Executive Summary

**YES** - `reaper_pm()` can and **SHOULD** be re-implemented using the existing `reaper_cpp()` C++ function.

### Key Finding:
The C++ implementation **already extracts pitch marks (epochs)** but `reaper_pm()` uses a slower Python path via `pyreaper` instead of leveraging the existing C++ infrastructure.

**Performance Impact:**
- Current Python `reaper_pm`: Slower (Python overhead + pyreaper)
- Potential C++ version: 2-3x faster (as stated in `trk_reaper` documentation)

**Implementation Status:**
- ✅ C++ backend **already complete** - `reaper_cpp()` extracts epochs
- ✅ Infrastructure **already exists** - `trk_reaper()` shows the pattern
- ⚠️ Current `reaper_pm()` **ignores** available C++ implementation
- 🎯 **Simple refactor** needed - reuse existing C++ code with different output format

---

## Current Implementation Comparison

### 1. reaper_pm() - Python Implementation

**File:** `R/ssff_python_reaper_pm.R`
**Backend:** Python `pyreaper` library
**Lines:** 212 lines

**Processing Flow:**
```
Audio File → av::read_audio_bin() → Convert to int16 → NumPy array
  → pyreaper.trk_reaper() [PYTHON] → Extract pm_times + pm
  → Build AsspDataObj → Return pitch marks
```

**Output:**
- Track: `pm` (pitch marks / glottal closure instants)
- Format: INT16
- Extension: `.rpm`
- Data: Pitch mark times and binary voicing flags

**Key Parameters:**
```r
reaper_pm(
  listOfFiles,
  beginTime = 0,
  endTime = 0,
  windowShift = 10,
  minF = 40,
  maxF = 500,
  unvoiced_cost = 0.9,
  high.pass = TRUE,
  hilbert.transform = FALSE,
  explicitExt = "rpm",
  outputDirectory = NULL,
  toFile = TRUE,
  conda.env = NULL  # Python environment!
)
```

**Dependencies:**
- Python + reticulate
- pyreaper Python library (`pip install pyreaper`)
- NumPy

**Performance:**
- ⚠️ Python overhead
- ⚠️ Data conversion: R → NumPy → pyreaper → R
- ⚠️ Requires Python environment management

### 2. trk_reaper() - C++ Implementation

**File:** `R/ssff_cpp_sptk_reaper.R`
**Backend:** C++ SPTK `reaper_cpp()`
**Lines:** 150 lines

**Processing Flow:**
```
Audio File → av_to_asspDataObj() → reaper_cpp() [C++]
  → Extract f0 + epochs + polarity → Build AsspDataObj
  → Return F0 (with epochs as attributes)
```

**Output:**
- Primary Track: `f0` (fundamental frequency)
- Extension: `.f0`
- Attributes: `epochs`, `n_epochs`, `polarity`
- Format: F0 values at regular intervals

**Key Parameters:**
```r
trk_reaper(
  listOfFiles,
  beginTime = 0.0,
  endTime = 0.0,
  windowShift = 10.0,
  minF = 60.0,
  maxF = 400.0,
  voicing_threshold = 0.9,  # Equivalent to unvoiced_cost
  toFile = TRUE,
  explicitExt = "f0",
  outputDirectory = NULL,
  verbose = TRUE
)
```

**Dependencies:**
- C++ SPTK (built-in, no external dependencies)
- No Python required

**Performance:**
- ✅ 2-3x faster than Python (per documentation)
- ✅ No Python overhead
- ✅ Direct C++ execution

### 3. reaper_cpp() - C++ Backend Function

**File:** `src/sptk_pitch.cpp`
**Function:** `reaper_cpp()`

**What It Returns:**
```cpp
List::create(
  Named("f0") = f0_matrix,           // F0 values [n_frames x 1]
  Named("times") = times,             // Frame times in seconds
  Named("sample_rate") = sample_rate, // Original sample rate
  Named("n_frames") = n_frames,       // Number of F0 frames
  Named("epochs") = epoch_times,      // ⭐ PITCH MARKS IN SECONDS
  Named("n_epochs") = n_epochs,       // ⭐ NUMBER OF PITCH MARKS
  Named("polarity") = polarity_str    // Signal polarity
)
```

**KEY INSIGHT:** The C++ function **already extracts epochs** (pitch marks)!

**What SPTK REAPER Does:**
```cpp
sptk::PitchExtractionByReaper reaper(...);
reaper.Get(waveform, &f0, &epochs, &polarity);
```

Extracts:
1. `f0` - Fundamental frequency time series
2. `epochs` - Glottal closure instants (pitch marks) **← THIS IS WHAT reaper_pm() WANTS!**
3. `polarity` - Signal polarity (positive/negative)

---

## Why Current Implementation is Suboptimal

### Problem: Redundant Python Path

`reaper_pm()` uses Python `pyreaper` when the **exact same algorithm** is already available in C++ via SPTK with **better performance**.

```
Current Path:
  R → Python → pyreaper → R     [SLOW, requires Python]

Better Path:
  R → C++ SPTK → R              [FAST, no dependencies]
```

### Missing Parameter Mappings

Python `pyreaper` parameters → C++ SPTK equivalents:

| Python (reaper_pm) | C++ (trk_reaper) | Status | Notes |
|-------------------|------------------|--------|-------|
| `windowShift` | `windowShift` | ✅ Identical | Frame period in ms |
| `minF` | `minF` | ✅ Identical | Min F0 in Hz |
| `maxF` | `maxF` | ✅ Identical | Max F0 in Hz |
| `unvoiced_cost` | `voicing_threshold` | ✅ Equivalent | Voicing decision threshold |
| `high.pass` | N/A | ⚠️ Missing | High-pass filter flag |
| `hilbert.transform` | ⚠️ Missing | High-pass filter flag |

**Notes:**
- `high.pass` and `hilbert.transform` are **preprocessing options** in pyreaper
- SPTK REAPER may handle these internally or not expose them
- Need to check SPTK documentation to see if these are configurable

---

## Proposed C++ Re-Implementation

### Option 1: Simple Wrapper (Recommended)

Create `trk_reaper_pm()` that uses existing `reaper_cpp()` but outputs pitch marks instead of F0.

**New File:** `R/ssff_cpp_sptk_reaper_pm.R`

**Pseudocode:**
```r
trk_reaper_pm <- function(listOfFiles,
                          beginTime = 0.0,
                          endTime = 0.0,
                          windowShift = 10.0,
                          minF = 40.0,
                          maxF = 500.0,
                          voicing_threshold = 0.9,
                          toFile = TRUE,
                          explicitExt = "rpm",
                          outputDirectory = NULL,
                          verbose = TRUE) {

  # Process each file
  for (file in listOfFiles) {
    # Load audio
    audio_obj <- av_to_asspDataObj(file, start_time = beginTime, end_time = endTime)

    # Call existing C++ function
    reaper_result <- reaper_cpp(
      audio_obj = audio_obj,
      minF = minF,
      maxF = maxF,
      windowShift = windowShift,
      voicing_threshold = voicing_threshold,
      verbose = FALSE
    )

    # Extract epochs (pitch marks) - ALREADY AVAILABLE!
    epoch_times <- reaper_result$epochs
    n_epochs <- reaper_result$n_epochs

    # Create SSFF object with pitch mark track
    out_obj <- create_pitchmark_asspobj(epoch_times, sample_rate, windowShift)

    # Write to file
    if (toFile) {
      write.AsspDataObj(out_obj, output_file)
    }
  }
}
```

**Advantages:**
- ✅ Reuses 100% of existing C++ code
- ✅ 2-3x performance improvement
- ✅ No Python dependency
- ✅ Minimal new code (~100 lines vs 212 in current implementation)
- ✅ Consistent with package architecture

### Option 2: Enhanced C++ Backend

Extend `reaper_cpp()` to expose `high.pass` and `hilbert.transform` parameters if SPTK supports them.

**Check SPTK Documentation:**
```cpp
// In SPTK source code, check if PitchExtractionByReaper supports:
sptk::PitchExtractionByReaper reaper(
  frame_shift_samples,
  sample_rate,
  minF,
  maxF,
  voicing_threshold,
  do_high_pass,        // ← Does this exist?
  do_hilbert_transform // ← Does this exist?
);
```

**If supported:**
- Add parameters to `reaper_cpp()` signature
- Update both `trk_reaper()` and new `trk_reaper_pm()` to expose them

**If not supported:**
- SPTK may handle preprocessing internally
- Document that high-pass and Hilbert transform are automatic
- Users won't have manual control, but that's acceptable

---

## Implementation Checklist

### Phase 1: Investigation (30 minutes)
- [ ] Check SPTK documentation for `PitchExtractionByReaper` parameters
- [ ] Verify if `do_high_pass` and `do_hilbert_transform` are available
- [ ] Compare epoch extraction output between pyreaper and SPTK

### Phase 2: Implementation (2-3 hours)
- [ ] Create `R/ssff_cpp_sptk_reaper_pm.R`
- [ ] Create helper function `create_pitchmark_asspobj()` (similar to `create_f0_asspobj()`)
- [ ] Implement `trk_reaper_pm()` following `trk_reaper()` pattern
- [ ] Set function attributes: `ext = "rpm"`, `tracks = "pm"`, `outputType = "SSFF"`
- [ ] Add comprehensive documentation with examples

### Phase 3: Testing (1-2 hours)
- [ ] Create test suite in `tests/testthat/test-reaper-pm.R`
- [ ] Compare output with current `reaper_pm()` Python implementation
- [ ] Verify epoch times match (within tolerance)
- [ ] Test edge cases (short files, silence, extreme F0 values)
- [ ] Benchmark performance vs Python version

### Phase 4: Deprecation (30 minutes)
- [ ] Mark `reaper_pm()` as deprecated with `.Deprecated()`
- [ ] Point users to `trk_reaper_pm()` in deprecation message
- [ ] Update documentation to recommend C++ version
- [ ] Plan removal timeline (e.g., v0.11.0)

### Phase 5: Documentation (1 hour)
- [ ] Update CLAUDE.md with reaper_pm → trk_reaper_pm migration
- [ ] Add performance comparison to documentation
- [ ] Update function catalog in PACKAGE_AUDIT
- [ ] Add migration guide to NEWS.md

---

## Performance Estimation

### Current Python Implementation
```
100 audio files × 10 seconds each:
- reaper_pm(): ~300-500 seconds (Python overhead)
```

### Proposed C++ Implementation
```
100 audio files × 10 seconds each:
- trk_reaper_pm(): ~100-200 seconds (2-3x faster)
```

**Time Saved:** 200-300 seconds for batch processing

---

## Compatibility Considerations

### Output Format Differences

**Current `reaper_pm()` Output:**
```
Track: pm (INT16)
- Binary pitch mark indicator (0 or 1)
- At regular intervals (windowShift)
- Values: 0 (no pitch mark), >0 (pitch mark present)
```

**Proposed `trk_reaper_pm()` Output:**
```
Track: pm (could be REAL64 or INT16)
- Option A: Epoch times in seconds (irregular intervals)
- Option B: Binary indicator at regular intervals (same as current)
```

**Recommendation:**
- Maintain **Option B** for backward compatibility
- Convert epoch times → binary indicator at windowShift intervals
- Optionally store raw epoch times as attribute (like `trk_reaper` does)

### Migration Path

For users currently using `reaper_pm()`:

```r
# Old code (Python)
reaper_pm("audio.wav", minF = 40, maxF = 500)

# New code (C++)
trk_reaper_pm("audio.wav", minF = 40, maxF = 500)

# Output format: IDENTICAL
# Performance: 2-3x FASTER
# Dependencies: NO PYTHON NEEDED
```

**Breaking Changes:** None if output format matches exactly

---

## Benefits of C++ Re-Implementation

### 1. Performance
- ⚡ **2-3x faster** processing (documented in trk_reaper)
- ⚡ No Python/R data conversion overhead
- ⚡ Direct C++ execution via SPTK

### 2. Reduced Dependencies
- ✅ No Python required
- ✅ No `pyreaper` installation needed
- ✅ No `conda.env` parameter needed
- ✅ Works out of the box after installing superassp

### 3. Consistency
- ✅ Matches package architecture (prefer C++ over Python)
- ✅ Consistent with `trk_reaper()`, `trk_dio()`, `trk_harvest()`, etc.
- ✅ Same parameter naming conventions
- ✅ Unified error handling and progress reporting

### 4. Maintainability
- ✅ One less Python dependency to maintain
- ✅ Leverages existing, tested C++ code
- ✅ Easier to debug (R + C++ vs R + Python + pyreaper)
- ✅ Consistent build process

### 5. User Experience
- ✅ Faster batch processing
- ✅ No need to manage Python environments
- ✅ Consistent API across all REAPER functionality
- ✅ Better integration with superassp ecosystem

---

## Potential Challenges

### 1. Parameter Equivalence
**Challenge:** `high.pass` and `hilbert.transform` from pyreaper might not map to SPTK

**Solution:**
- Check SPTK source code for available parameters
- If not available, document that SPTK handles preprocessing automatically
- Most users don't change these from defaults anyway

### 2. Output Format Matching
**Challenge:** Epoch times (irregular) vs binary indicator (regular intervals)

**Solution:**
- Implement conversion: epochs → binary grid at windowShift intervals
- Store raw epochs as attribute for advanced users
- Maintain exact compatibility with current output format

### 3. Numerical Precision
**Challenge:** pyreaper vs SPTK might have slight algorithmic differences

**Solution:**
- Both implement Google's REAPER algorithm
- Differences should be minimal (floating-point precision)
- Test on real data and document any differences
- If significant, provide option to choose backend

---

## Recommendation

### STRONGLY RECOMMENDED: Implement C++ Version

**Why:**
1. ✅ **Performance:** 2-3x speedup is significant for batch processing
2. ✅ **Dependencies:** Eliminates Python requirement
3. ✅ **Consistency:** Aligns with package architecture
4. ✅ **Maintainability:** Less code, fewer dependencies
5. ✅ **Code Reuse:** `reaper_cpp()` already does the work!

**Effort vs Benefit:**
- Effort: ~4-6 hours total (investigation + implementation + testing)
- Benefit: Significant performance improvement, reduced dependencies, better UX

**Implementation Priority:** **HIGH**

This is a clear win - the C++ code already exists and extracts exactly what `reaper_pm()` needs. It's simply a matter of reformatting the output.

---

## Next Steps

1. **Immediate:** Verify SPTK REAPER parameter support
   ```bash
   # Check SPTK source for available parameters
   grep -r "PitchExtractionByReaper" src/SPTK/
   ```

2. **Create Implementation Plan:**
   - Design `create_pitchmark_asspobj()` helper
   - Draft `trk_reaper_pm()` function
   - Plan test suite

3. **Prototype:**
   - Quick prototype to verify epoch → binary conversion
   - Compare output with current `reaper_pm()`
   - Benchmark performance

4. **Full Implementation:**
   - Complete function with documentation
   - Comprehensive test suite
   - Deprecate `reaper_pm()`

---

## Conclusion

**YES - reaper_pm() should be re-implemented in C++**

The infrastructure already exists (`reaper_cpp()` extracts epochs), making this a straightforward refactor with significant benefits:
- 2-3x performance improvement
- Eliminates Python dependency
- Better consistency with package architecture
- Minimal implementation effort (reuse existing code)

This is a **high-priority, low-effort, high-benefit** improvement that should be implemented soon.

---

**Analysis Complete**

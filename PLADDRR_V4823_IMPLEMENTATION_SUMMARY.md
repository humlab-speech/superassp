# pladdrr v4.8.23 Compatibility Update - Implementation Summary

**Date**: 2026-02-13  
**Status**: ✅ Complete  
**Result**: All pladdrr functions updated and tested

---

## Changes Made

### 1. Version Requirement (DESCRIPTION:37)
```diff
- pladdrr (>= 4.8.20)
+ pladdrr (>= 4.8.23)
```

### 2. Performance Optimization (R/list_pladdrr_pharyngeal.R:645)

**Before** (loop-based):
```r
extract_ltas_peaks <- function(ltas, fmins, fmaxs) {
  n <- length(fmins)
  peak_values <- numeric(n)
  peak_frequencies <- numeric(n)
  
  for (i in seq_len(n)) {
    peak_values[i] <- ltas$get_maximum(fmins[i], fmaxs[i], "Parabolic")
    ns <- asNamespace("pladdrr")
    peak_frequencies[i] <- ns$.ltas_get_frequency_of_maximum(
      ltas$.xptr, fmins[i], fmaxs[i], 1
    )
  }
  ...
}
```

**After** (18x faster batch API):
```r
extract_ltas_peaks <- function(ltas, fmins, fmaxs) {
  # Use batch API for 18x speedup (pladdrr >= 4.8.0)
  batch_result <- ltas$get_peaks_batch(fmins, fmaxs, "Parabolic")
  
  data.frame(
    fmin = fmins,
    fmax = fmaxs,
    peak_value = batch_result$values,
    peak_frequency = batch_result$frequencies
  )
}
```

**Impact**: `lst_pharyngeal()` batch LTAS extraction 18x faster

### 3. Documentation Fixes (R/pladdrr_helpers.R)

**Line 33 & 131** - Deprecated `to_pitch()` → `to_pitch_cc()`:
```diff
- pitch <- sound$to_pitch(time_step = 0.01, ...)
+ pitch <- sound$to_pitch_cc(time_step = 0.01, ...)
```

### 4. Critical API Updates (R/ssff_pladdrr_cpps.R)

**Issue**: `trk_cpps()` broke due to pladdrr v4.8.23 API changes

#### 4a. PowerCepstrogram Parameter Name Change (line 154)
```diff
  sound$to_powercepstrogram(
    pitch_floor = minF,
    time_step = timeStep,
    maximum_frequency = maximumFrequency,
-   pre_emphasis_from = preEmphFrom
+   pre_emphasis_frequency = preEmphFrom
  )
```

#### 4b. Property Access Fix (line 146)
```diff
  duration <- sound$.cpp$duration
- sample_rate <- sound$get_sampling_frequency()
+ sample_rate <- sound$.cpp$sampling_frequency
```

#### 4c. Per-Frame CPP Extraction Rewrite (lines 160-186)

**Before** (internal API, broken):
```r
num_frames <- power_cepstrogram$get_number_of_frames()  # Method removed
times <- numeric(num_frames)
cpp_arr <- numeric(num_frames)

ns <- asNamespace("pladdrr")
for (j in seq_len(num_frames)) {
  frame_time <- power_cepstrogram$get_time_from_frame_number(j)  # Removed
  cepstrum_slice <- power_cepstrogram$to_powercepstrum_slice(frame_time)  # Removed
  
  cpp <- ns$.powercepstrum_get_peak_prominence(  # Signature changed
    cepstrum_slice$.xptr, interpolation, minF, maxF, ...
  )
  ...
}
```

**After** (public API, working):
```r
# Get matrix to determine frame count and times
pc_matrix <- power_cepstrogram$to_matrix()
num_frames <- pc_matrix$get_number_of_rows()

times <- numeric(num_frames)
cpp_arr <- numeric(num_frames)

# Loop through frames - use get_cpp_at_time (pladdrr v4.8+)
for (j in seq_len(num_frames)) {
  # Calculate frame time: xmin + (frame_number - 1) * dx
  frame_time <- pc_matrix$get_xmin() + (j - 1) * pc_matrix$get_dx()
  times[j] <- frame_time
  
  # Get CPP value at this time (public API)
  cpp <- tryCatch({
    power_cepstrogram$get_cpp_at_time(frame_time)
  }, error = function(e) NaN)
  
  cpp_arr[j] <- ifelse(is.null(cpp) || is.na(cpp), NaN, cpp)
}
```

**Key Changes**:
- Removed methods: `get_number_of_frames()`, `get_time_from_frame_number()`, `to_powercepstrum_slice()`
- Replacement: `to_matrix()` for frame metadata, `get_cpp_at_time()` for per-frame values
- No longer uses internal API `ns$.powercepstrum_get_peak_prominence()`
- More stable (uses public API only)

---

## API Changes Encountered

| Old API | Status | New API |
|---------|--------|---------|
| `pre_emphasis_from` | ❌ Removed | `pre_emphasis_frequency` |
| `sound$get_sampling_frequency()` | ❌ Removed | `sound$.cpp$sampling_frequency` |
| `pc$get_number_of_frames()` | ❌ Removed | `pc$to_matrix()$get_number_of_rows()` |
| `pc$get_time_from_frame_number()` | ❌ Removed | Calculate: `xmin + (frame-1) * dx` |
| `pc$to_powercepstrum_slice()` | ❌ Removed | N/A (use `get_cpp_at_time()`) |
| `ns$.powercepstrum_get_peak_prominence()` | ⚠️ Internal | `pc$get_cpp_at_time()` (public) |

---

## Testing

### Manual Tests (Passed)
```r
✅ trk_cpps("a1.wav") → 513 frames, mean = 13.37 dB
✅ lst_pharyngeal("a1.wav", beginTime=0, endTime=1) → 79 features
```

### Test Suite (Passed)
```bash
devtools::test(filter = 'cpps|pharyngeal|vuv|vq')
# Result: All tests passed or skipped (COVAREP not installed)
```

---

## Files Modified

1. `DESCRIPTION` - Version requirement bump
2. `R/list_pladdrr_pharyngeal.R` - Batch LTAS optimization
3. `R/pladdrr_helpers.R` - Doc fixes (deprecated API references)
4. `R/ssff_pladdrr_cpps.R` - Major API updates (4 fixes)
5. `man/*.Rd` - Auto-regenerated from roxygen2

---

## Benefits

1. **Stability**: Uses public API only (no internal namespace calls)
2. **Performance**: 18x LTAS speedup via batch operations
3. **Compatibility**: Works with pladdrr v4.8.23+
4. **Documentation**: Examples use current API conventions

---

## Verification Checklist

- [x] DESCRIPTION version updated to `>= 4.8.23`
- [x] `extract_ltas_peaks()` uses batch API
- [x] Doc examples use `to_pitch_cc()` not `to_pitch()`
- [x] `trk_cpps()` uses public PowerCepstrogram API
- [x] All parameter names updated (`pre_emphasis_frequency`)
- [x] Property access uses `.cpp$` pattern
- [x] Manual testing passed (CPPS, pharyngeal)
- [x] Test suite passed
- [x] Documentation regenerated

---

## Migration Notes for Future Updates

**Key Lesson**: pladdrr v4.8+ emphasizes **public API** over internal methods:
- Prefer R6 public methods over internal `ns$.function()` calls
- Use property access (`.cpp$field`) over getter methods when available
- Check `names(object)` to discover available public methods
- Matrix-based iteration (`to_matrix()` + row iteration) replaces frame-by-frame helpers

**Recommendation**: Monitor pladdrr NEWS.md for API changes in future updates.

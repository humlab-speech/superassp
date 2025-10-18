# Kaldi Pitch Update Analysis

**Date**: 2025-10-18  
**Issue**: `kaldi_pitch()` uses deprecated `compute_kaldi_pitch()` removed in torchaudio 2.9+

## Current Situation

### Kaldi Pitch Implementation
- **File**: `R/ssff_python_kaldi_pitch.R`
- **Current Status**: Uses `torchaudio.functional.compute_kaldi_pitch` (deprecated)
- **Problem**: Function removed in torchaudio 2.9+ (current version: 2.8.0)
- **Output**: F0 + NCCF (Normalized Cross-Correlation Function)

### Current Torchaudio Status
- **Version**: 2.8.0
- **compute_kaldi_pitch**: ❌ Not available
- **detect_pitch_frequency**: ✅ Available (replacement)

## Replacement Function

### torchaudio.functional.detect_pitch_frequency
Available in torchaudio >= 0.13.0 (including current 2.8.0)

**Signature**:
```python
detect_pitch_frequency(
    waveform: Tensor,
    sample_rate: int,
    frame_time: float = 0.01,
    win_length: int = 30,
    freq_low: int = 85,
    freq_high: int = 3400
) -> Tensor
```

**Key Differences from compute_kaldi_pitch**:
1. ✅ **Simpler interface** - fewer parameters
2. ❌ **No NCCF output** - only returns F0 values
3. ✅ **Same algorithm** - NCCF + median smoothing
4. ✅ **Compatible with current torch** (2.8.0)

## Recommendation

### Update `kaldi_pitch()` to use `detect_pitch_frequency`

**Changes needed**:
1. Replace Python code to call `detect_pitch_frequency` instead of `compute_kaldi_pitch`
2. Remove NCCF track from output (only F0)
3. Simplify parameters (remove Kaldi-specific options)
4. Update documentation to note NCCF unavailable
5. Follow pattern from other pitch trackers (rapt, swipe, etc.)

**Benefits**:
- ✅ Works with current torchaudio 2.8.0
- ✅ No need to downgrade torchaudio
- ✅ Simpler implementation
- ✅ Consistent with package patterns

**Trade-offs**:
- ❌ Loses NCCF output (voicing probability)
- Users needing NCCF must use alternatives (rapt, reaper)

## Implementation Pattern

Follow the pattern from `ssff_cpp_sptk_rapt.R`:

```r
kaldi_pitch <- function(listOfFiles,
                       beginTime = 0,
                       endTime = 0,
                       windowShift = 10,
                       windowSize = 30,
                       minF = 85,
                       maxF = 400,
                       toFile = TRUE,
                       explicitExt = "kap",
                       outputDirectory = NULL,
                       verbose = TRUE) {
  # Standard pattern:
  # 1. Validate inputs
  # 2. Load audio via torchaudio
  # 3. Call F.detect_pitch_frequency()
  # 4. Convert to AsspDataObj with F0 track
  # 5. Write SSFF file if toFile=TRUE
}
```

## Alternative Solutions

If NCCF output is required:

###  Option 1: Keep old version for legacy users
- Document requirement: torchaudio < 2.9
- Add deprecation warning
- Suggest alternative: `rapt()`, `reaper()`

### Option 2: Implement NCCF computation separately
- Use `torch` to compute NCCF manually
- More complex but preserves functionality
- May not match original Kaldi exactly

### Option 3: Direct alternative recommendation
- **rapt()** - SPTK C++ implementation, includes voicing
- **reaper()** - SPTK C++ implementation, includes voicing  
- **crepe()** - Deep learning, includes confidence

## Recommended Action

**Update kaldi_pitch() to use detect_pitch_frequency**:

1. Simplify to match other pitch trackers in package
2. Remove NCCF track (document change)
3. Update to work with current torchaudio 2.8.0
4. Add note recommending `rapt()` or `reaper()` for voicing probability

This maintains functionality while working with modern torch versions.

---

**Status**: Analysis complete, ready for implementation  
**Priority**: Medium - function currently broken with torchaudio 2.9+  
**Estimated effort**: 2-3 hours

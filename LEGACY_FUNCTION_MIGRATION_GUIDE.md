# Legacy Function Migration Guide

## Overview
This guide provides step-by-step migration instructions for functions that should be deprecated or modernized.

## Functions to Deprecate

### 1. trk_snackp - Snack Pitch Tracking

**Status:** Legacy reference implementation  
**Location:** `/Users/frkkan96/Documents/src/superassp/R/ssff_python_snack_pitch.R`  
**Performance:** 500-1000ms per 3s audio  

**Why Deprecate:**
- Slower than modern alternatives
- Autocorrelation-only method (less robust than multi-method approaches)
- Not significantly better than SPTK alternatives

**Better Alternative:** `trk_rapt`
- Performance: <100ms (6-8x faster)
- Also autocorrelation-based, similar approach
- Modern C++ implementation with av package support

**Migration Steps:**
```r
# Old code (deprecated)
f0 <- trk_snackp("audio.wav", minF = 50, maxF = 550)

# New code (recommended)
f0 <- trk_rapt("audio.wav", minF = 50, maxF = 550)
```

**Timeline:**
- v0.8.8: Add @deprecated tag, update docs
- v0.9: Officially deprecate with warning
- v1.0: Remove function

---

### 2. trk_snackf - Snack Formant Tracking

**Status:** Legacy reference implementation  
**Location:** `/Users/frkkan96/Documents/src/superassp/R/ssff_python_snack_formant.R`  
**Performance:** 500-1000ms per 3s audio  

**Why Deprecate:**
- Slower than modern alternatives
- Limited to specific use case

**Better Alternative:** `trk_forest`
- Performance: 200-400ms (2-3x faster)
- Better tracking algorithm (spectral peak picking)
- Recommended for general use

**Migration Steps:**
```r
# Old code (deprecated)
formants <- trk_snackf("audio.wav", gender = "m")

# New code (recommended)
formants <- trk_forest("audio.wav", gender = "m")
```

**Timeline:** Same as trk_snackp

---

### 3. trk_straight_f0 - STRAIGHT F0 Extraction

**Status:** Legacy file-based pipeline  
**Location:** `/Users/frkkan96/Documents/src/superassp/R/ssff_python_straight_f0.R`  
**Performance:** 5-10s per 3s audio  

**Why Deprecate:**
- Extremely slow (10-20x slower than modern alternatives)
- File-based architecture (legacy)
- Most research no longer requires exact STRAIGHT results

**Better Alternatives:**
- `trk_harvest` - Similar vocoder-based approach but modern/fast (100-200ms)
- `trk_reaper` - Fast and robust (100-200ms)
- `trk_rapt` - For simple autocorrelation (faster)

**Migration Steps:**
```r
# Old code (deprecated)
f0 <- trk_straight_f0("audio.wav", f0_floor = 50, f0_ceil = 500)

# New code (recommended)
# For closest algorithm match:
f0 <- trk_harvest("audio.wav", f0_floor = 50, f0_ceil = 500)

# Or for fastest autocorrelation:
f0 <- trk_rapt("audio.wav", minF = 50, maxF = 500)
```

**Timeline:**
- v0.8.8: Add @deprecated tag
- v0.9: Officially deprecate
- v1.1: Remove (longer timeline for publications)

---

### 4. trk_straight_spec - STRAIGHT Spectral Analysis

**Status:** Legacy file-based pipeline  
**Location:** `/Users/frkkan96/Documents/src/superassp/R/ssff_python_straight_spec.R`  
**Performance:** 5-10s per 3s audio  

**Why Deprecate:**
- File-based, slow spectral analysis
- Slower alternatives available

**Better Alternatives:**
- `trk_dftSpectrum` - DFT-based spectral analysis (200-400ms)
- `trk_cssSpectrum` - Chebyshev spectral analysis (200-400ms)
- `trk_cepstrum` - Cepstral analysis (200-400ms)

**Migration Steps:**
```r
# Old code (deprecated)
spectrum <- trk_straight_spec("audio.wav")

# New code (recommended) - DFT equivalent
spectrum <- trk_dftSpectrum("audio.wav")

# Or for Chebyshev (similar spectral properties)
spectrum <- trk_cssSpectrum("audio.wav")
```

**Timeline:** Same as trk_straight_f0

---

### 5. straight_pipeline - STRAIGHT Pipeline Wrapper

**Status:** Legacy pipeline orchestrator  
**Location:** `/Users/frkkan96/Documents/src/superassp/R/ssff_python_straight_*.R`  
**Performance:** 10-20s total

**Why Remove:**
- Wrapper for deprecated functions
- No longer needed when individual functions modernized
- Adds complexity

**Better Approach:**
- Call individual functions as needed
- Better control and modularity

**Migration Steps:**
```r
# Old code (deprecated)
result <- straight_pipeline("audio.wav", 
                           f0_floor = 50, 
                           f0_ceil = 500,
                           spectrum = TRUE,
                           synthesize = FALSE)

# New code (recommended)
f0 <- trk_harvest("audio.wav", f0_floor = 50, f0_ceil = 500)
spectrum <- trk_dftSpectrum("audio.wav")
# ... or other desired analyses
```

**Timeline:**
- v0.9: Remove

---

## Functions Needing Modernization

### 6. trk_pyin - Probabilistic YIN Pitch Tracking

**Status:** Needs av package migration  
**Location:** `/Users/frkkan96/Documents/src/superassp/R/ssff_python_pyin.R`  
**Issue:** Uses librosa.load instead of av::read_audio_bin

**Current Code Pattern (PARTIAL MODERNIZATION):**
```r
# Line ~98-100
audio_data <- av::read_audio_bin(
  audio = origSoundFile,
  start_time = if (beginTime > 0) beginTime else NULL,
  ...
)
```

**Problem:** Some functions still use librosa internally for Python processing

**Action Required:**
- Ensure all audio loading uses av package
- Remove librosa dependency
- Verify in-memory processing throughout

**Migration Timeline:**
- v0.8.8: Complete av integration
- v0.9: Standardize across all functions

---

### 7. trk_yin - YIN Pitch Tracking

**Status:** Needs av package migration  
**Location:** `/Users/frkkan96/Documents/src/superassp/R/ssff_python_yin.R`  
**Issue:** Uses librosa.load instead of av::read_audio_bin

**Action:** Same as trk_pyin

---

## Implementation Checklist

### For Each Deprecated Function:

- [ ] Add @deprecated roxygen2 tag
- [ ] Update @return documentation
- [ ] Add migration recommendation in description
- [ ] Link to replacement function in @seealso
- [ ] Update example to show recommended alternative
- [ ] Add to NEWS.md with deprecation notice
- [ ] Send deprecation warnings in v0.9

### Example @deprecated Tag:
```r
#' @deprecated This function is deprecated as of superassp v0.9.
#'   Use \code{\link{trk_rapt}} instead for better performance.
#'   STRAIGHT results are rarely required in modern research.
#'   For details, see \code{?trk_straight_f0} for migration guide.
```

### Example in Function:
```r
trk_snackp <- function(...) {
  .Deprecated("trk_rapt", 
              msg = "trk_snackp() is deprecated. Use trk_rapt() instead (6-8x faster)")
  # ... function body
}
```

---

## Performance Comparison

### Pitch Tracking
| Function | Performance | Recommendation | Status |
|----------|-------------|-----------------|--------|
| trk_snackp | 500-1000ms | Deprecate | Legacy |
| trk_straight_f0 | 5-10s | Deprecate | Legacy |
| **trk_rapt** | **<100ms** | **Use** | **Modern** |
| **trk_harvest** | **100-200ms** | **Use** | **Modern** |
| **trk_crepe** | **500-2000ms** | **Use** | **Modern (specialized)** |

### Formant Tracking
| Function | Performance | Recommendation | Status |
|----------|-------------|-----------------|--------|
| trk_snackf | 500-1000ms | Deprecate | Legacy |
| **trk_forest** | **200-400ms** | **Use** | **Modern** |
| **trk_formantp** | **500-800ms** | **Use** | **Modern (Praat)** |

### Spectral Analysis
| Function | Performance | Recommendation | Status |
|----------|-------------|-----------------|--------|
| trk_straight_spec | 5-10s | Deprecate | Legacy |
| **trk_dftSpectrum** | **200-400ms** | **Use** | **Modern** |
| **trk_cssSpectrum** | **200-400ms** | **Use** | **Modern** |

---

## Version Timeline

### v0.8.8 (Immediate)
- Migrate trk_pyin, trk_yin to complete av integration
- Mark STRAIGHT functions with deprecation tags in docs
- Update all docstrings with migration information
- No breaking changes

### v0.9 (Short-term, ~2-3 months)
- Officially deprecate: trk_snackp, trk_snackf, trk_straight_f0, trk_straight_spec, straight_pipeline
- Add .Deprecated() warnings in functions
- Update NEWS.md with detailed migration guide
- Update documentation site
- Send user notifications

### v1.0 (Long-term, ~6 months)
- Remove deprecated functions
- All remaining functions follow modern workflow
- Complete av package integration
- Standardized parameter naming

---

## Testing Deprecations

### Before Release:
```r
# Test that deprecated functions still work
expect_no_error(trk_snackp("test.wav", toFile = FALSE))

# Test that warnings are raised
expect_warning(trk_snackp("test.wav", toFile = FALSE), 
               "deprecated")

# Test that replacement works
old_result <- suppressWarnings(trk_snackp("test.wav", toFile = FALSE))
new_result <- trk_rapt("test.wav", toFile = FALSE)
# Results should be similar (both autocorrelation-based)
```

### Documentation Testing:
```r
# Ensure links work
?trk_snackp
# Should show @seealso link to trk_rapt
```

---

## Communication Strategy

### Before Deprecation (v0.8.8)
- Update documentation with performance comparisons
- Add migration examples in function docs
- No user-facing changes

### At Deprecation (v0.9)
- Include deprecation warnings in output
- Update NEWS.md with detailed guide
- Add migration section to vignettes
- Update website documentation
- Consider blog post or announcement

### Before Removal (v1.0)
- Final warning that functions will be removed
- Ensure all examples use modern alternatives
- Provide comprehensive migration guide in docs

---

## Contact & Support

For questions about function migration:
1. Check this guide first
2. See `?function_name` for updated documentation
3. Review NEWS.md for version changes
4. Check GitHub issues for discussions
5. Contact maintainers for complex migrations

---

**Last Updated:** 2025-10-29  
**Status:** Approved for implementation  
**Approval:** Code quality audit recommends implementation v0.8.8


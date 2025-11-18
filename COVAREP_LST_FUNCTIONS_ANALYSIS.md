# COVAREP lst_ Functions Analysis

## Overview

Analysis of COVAREP Python module to identify functions suitable for R `lst_` integration (summary measures without time-series output).

**Analysis Date:** October 19, 2025

---

## Summary: Identified lst_ Function Candidates

### Primary Candidate: Voice Quality Parameters ✅

**Python Function:** `extract_vq_params_optimized()`
**Proposed R Function:** `lst_covarep_vq()`
**Priority:** HIGH
**Status:** Python implementation ready, needs R wrapper

---

## Detailed Analysis

### 1. Voice Quality Parameters - `lst_covarep_vq()`

**Python Function:** `extract_vq_params_optimized()` from `covarep/glottal/vq_optimized.py`

**Input:**
```python
extract_vq_params_optimized(
    glottal_flow,        # From IAIF
    glottal_derivative,  # From IAIF
    fs,                  # Sampling rate
    f0=None,            # Optional: F0 estimate
    gci=None            # Optional: Glottal Closure Instants
)
```

**Output:** Dictionary with 7+ parameters
```python
{
    'glottal_flow_max': float,              # Peak glottal flow amplitude
    'glottal_flow_min': float,              # Minimum glottal flow
    'glottal_derivative_peak': float,       # Peak flow derivative (MFDR)
    'NAQ': float,                           # Normalized Amplitude Quotient*
    'QOQ': float,                           # Quasi-Open Quotient*
    'H1_H2': float,                         # H1-H2 harmonic difference (dB)**
    'HRF': float,                           # Harmonic Richness Factor
    'PSP': float                            # Parabolic Spectral Parameter
}

* Requires GCI (Glottal Closure Instants) - returns NaN if not provided
** Requires F0 - returns NaN if not provided
```

**Clinical/Scientific Significance:**

1. **NAQ (Normalized Amplitude Quotient)**
   - Measures vocal fold closure characteristics
   - Clinical use: Dysphonia assessment, breathiness detection
   - Normal range: 0.05-0.20
   - Higher values → more breathy voice

2. **QOQ (Quasi-Open Quotient)**
   - Ratio of open phase to total glottal cycle
   - Clinical use: Voice quality assessment
   - Normal range: 0.4-0.7
   - Related to vocal effort and voice quality

3. **H1-H2**
   - Spectral tilt measure (dB difference between first two harmonics)
   - Clinical use: Voice quality, breathiness
   - Normal range: -10 to +10 dB
   - Positive values → breathy, negative → pressed phonation

4. **HRF (Harmonic Richness Factor)**
   - Ratio of high-frequency to low-frequency energy
   - Clinical use: Voice quality characterization
   - Indicates spectral richness

5. **PSP (Parabolic Spectral Parameter)**
   - Spectral envelope curvature
   - Clinical use: Voice quality assessment
   - Related to overall spectral shape

**Optimization Status:**
- ✅ Numba JIT available for NAQ/QOQ computation
- ✅ Vectorized implementations for H1-H2, HRF, PSP
- ✅ 2-3x speedup over pure Python

**Dependencies:**
- Required: Glottal flow and derivative from IAIF
- Optional: F0 estimates (for H1-H2)
- Optional: GCI detection (for NAQ/QOQ)

**Current Limitations:**
- GCI detection not yet implemented (SEDREAMS/DYPSA needed)
- NAQ/QOQ return NaN without GCI
- H1-H2 returns NaN without F0

**Proposed R Function Signature:**

```r
lst_covarep_vq <- function(listOfFiles,
                           beginTime = 0.0,
                           endTime = 0.0,
                           f0 = NULL,           # Optional F0 contour or scalar
                           gci = NULL,          # Optional GCI timestamps
                           compute_iaif = TRUE, # Auto-compute IAIF if needed
                           verbose = TRUE) {
  # Returns list of lists (one per file) with:
  # - glottal_flow_max
  # - glottal_flow_min
  # - glottal_derivative_peak
  # - NAQ (if GCI provided)
  # - QOQ (if GCI provided)
  # - H1_H2 (if F0 provided)
  # - HRF
  # - PSP
}
```

**Example Usage:**

```r
# Basic usage (without GCI or F0)
vq <- lst_covarep_vq("vowel.wav")
# $glottal_flow_max: 0.85
# $glottal_flow_min: -0.12
# $glottal_derivative_peak: 2.34
# $NAQ: NA (no GCI)
# $QOQ: NA (no GCI)
# $H1_H2: NA (no F0)
# $HRF: -8.5
# $PSP: -0.023

# With F0 for H1-H2
f0_data <- trk_covarep_srh("vowel.wav", toFile = FALSE)
f0_mean <- mean(f0_data$`F0[Hz]`[f0_data$VUV == 1, 1])
vq <- lst_covarep_vq("vowel.wav", f0 = f0_mean)
# $H1_H2: 5.2 (now computed)

# With F0 contour (median used)
vq <- lst_covarep_vq("vowel.wav", f0 = f0_data$`F0[Hz]`[, 1])

# Batch processing
files <- c("a.wav", "e.wav", "i.wav")
vq_all <- lst_covarep_vq(files)
# Returns list of 3 parameter lists
```

**Implementation Priority:** HIGH

**Rationale:**
- ✅ Python function already implemented and optimized
- ✅ Returns summary measures (perfect for `lst_`)
- ✅ High clinical/research value
- ✅ Complements existing `trk_covarep_iaif()`
- ⚠️ Some parameters require additional inputs (F0, GCI)

---

## Functions NOT Suitable for lst_

### Why Existing Functions Don't Qualify

**1. `pitch_srh()` - Already integrated as `trk_covarep_srh()`**
- Returns: F0 contour (time-series)
- Output: `(n_frames,)` arrays for f0, vuv, srh, times
- Classification: Track function ✓ (already done)

**2. `iaif()` - Already integrated as `trk_covarep_iaif()`**
- Returns: Glottal flow waveform (time-series)
- Output: `(n_samples,)` arrays for flow and derivative
- Classification: Track function ✓ (already done)

**3. `get_vq_params()` - Original version**
- Status: Superseded by `extract_vq_params_optimized()`
- Less complete than optimized version
- Missing HRF and PSP
- Recommendation: Use optimized version instead

**4. Voicebox utility functions**
- `frq2mel()`, `mel2frq()`, etc.: Frequency conversions (not file-based)
- `enframe()`: Frame extraction utility (not summary measure)
- `zerocros()`, `activlev()`: Could be summary measures but very basic
- Classification: Utilities, not worth dedicated R functions

---

## Future lst_ Function Opportunities

### If/When Additional COVAREP Modules Are Ported

Based on MATLAB COVAREP, potential future `lst_` functions:

**1. `lst_covarep_gci()` - GCI Detection Summary**
- Would return: GCI count, mean period, period variability
- Requires: SEDREAMS or DYPSA implementation
- Priority: MEDIUM (enables NAQ/QOQ)

**2. `lst_covarep_complexity()` - Voice Complexity Measures**
- Would return: Sample entropy, correlation dimension, LZ complexity
- Requires: Nonlinear dynamics module
- Priority: LOW (research-oriented)

**3. `lst_covarep_spectral()` - Spectral Moments**
- Would return: Spectral centroid, spread, skewness, kurtosis
- Requires: Spectral analysis module
- Priority: LOW (basic measures, can be computed elsewhere)

---

## Implementation Recommendation

### Immediate Action: Implement `lst_covarep_vq()`

**Phase 1: Basic Implementation (1-2 days)**

```r
# File: R/covarep_vq.R

lst_covarep_vq <- function(listOfFiles,
                           beginTime = 0.0,
                           endTime = 0.0,
                           f0 = NULL,
                           gci = NULL,
                           verbose = TRUE) {

  # Check COVAREP availability
  if (!covarep_available()) {
    stop("COVAREP not available. Install with: install_covarep()")
  }

  # Process files
  results <- lapply(listOfFiles, function(file) {
    # Load audio
    audio <- av_load_for_python(file, beginTime, endTime)

    # Compute IAIF
    iaif_result <- covarep_module$glottal$iaif_optimized$iaif_optimized(
      x = audio$samples,
      fs = audio$sample_rate
    )

    # Prepare F0 (scalar or NULL)
    py_f0 <- if (!is.null(f0)) {
      if (length(f0) == 1) f0 else median(f0[f0 > 0])
    } else NULL

    # Prepare GCI (convert to sample indices if needed)
    py_gci <- if (!is.null(gci)) {
      as.integer(gci * audio$sample_rate)
    } else NULL

    # Extract VQ parameters
    vq_params <- covarep_module$glottal$vq_optimized$extract_vq_params_optimized(
      glottal_flow = iaif_result[[1]],
      glottal_derivative = iaif_result[[2]],
      fs = audio$sample_rate,
      f0 = py_f0,
      gci = py_gci
    )

    # Convert to R list
    as.list(vq_params)
  })

  # Simplify if single file
  if (length(results) == 1) results[[1]] else results
}
```

**Phase 2: Enhanced Implementation (Optional, 1-2 days)**

Add automatic F0 estimation integration:

```r
lst_covarep_vq <- function(...,
                           auto_f0 = TRUE,   # Auto-estimate F0 if not provided
                           f0_method = "srh", # "srh" or external
                           ...) {

  if (auto_f0 && is.null(f0)) {
    # Compute F0 using SRH
    f0_data <- trk_covarep_srh(file, toFile = FALSE, verbose = FALSE)
    f0 <- median(f0_data$`F0[Hz]`[f0_data$VUV == 1, 1])
  }

  # ... rest of implementation
}
```

**Phase 3: GCI Integration (Future, when GCI implemented)**

```r
lst_covarep_vq <- function(...,
                           auto_gci = TRUE,   # Auto-detect GCI
                           gci_method = "sedreams",
                           ...) {

  if (auto_gci && is.null(gci)) {
    # Compute GCI using SEDREAMS/DYPSA
    gci <- detect_gci(glottal_flow, glottal_derivative, fs, method = gci_method)
  }

  # ... rest of implementation
}
```

---

## Testing Strategy

**Test File:** `tests/testthat/test-covarep-vq.R`

**Test Cases (12):**

1. ✓ Module availability
2. ✓ Single file processing (basic)
3. ✓ Batch processing
4. ✓ With F0 (scalar)
5. ✓ With F0 (vector)
6. ✓ With GCI (when available)
7. ✓ Without F0 (returns NA for H1-H2)
8. ✓ Without GCI (returns NA for NAQ/QOQ)
9. ✓ Time windowing
10. ✓ Parameter validation
11. ✓ Consistency checks
12. ✓ Output structure validation

**Example Test:**

```r
test_that("lst_covarep_vq works with F0", {
  skip_if_not(covarep_available(), "COVAREP not installed")

  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")

  # Get F0
  f0_data <- trk_covarep_srh(test_wav, toFile = FALSE, verbose = FALSE)
  f0_mean <- median(f0_data$`F0[Hz]`[f0_data$VUV == 1, 1])

  # Extract VQ with F0
  vq <- lst_covarep_vq(test_wav, f0 = f0_mean, verbose = FALSE)

  expect_type(vq, "list")
  expect_true("H1_H2" %in% names(vq))
  expect_false(is.na(vq$H1_H2))  # Should be computed
  expect_true(is.na(vq$NAQ))     # No GCI, should be NA
})
```

---

## Performance Expectations

**With Numba optimization:**
- VQ extraction: ~5-10ms per file (after IAIF)
- Total time: ~30-50ms per file (IAIF + VQ)
- Batch processing: Near-linear scaling

**Without optimization:**
- VQ extraction: ~15-30ms per file
- Total time: ~80-150ms per file

---

## Documentation Needs

**Roxygen2 Documentation:**

```r
#' Voice Quality Parameter Extraction via COVAREP
#'
#' Extract comprehensive voice quality measures from speech signals using
#' COVAREP algorithms. Computes parameters related to glottal source
#' characteristics, spectral properties, and voice quality.
#'
#' @param listOfFiles Character vector of audio file paths
#' @param beginTime Numeric start time(s) in seconds
#' @param endTime Numeric end time(s) in seconds
#' @param f0 Optional F0 estimate (scalar or vector). If vector, median is used.
#'   Required for H1-H2 computation. Can be obtained from trk_covarep_srh().
#' @param gci Optional glottal closure instants (sample indices or times).
#'   Required for NAQ and QOQ computation.
#' @param verbose Logical; show progress messages
#'
#' @return List (or list of lists for multiple files) containing:
#'   \describe{
#'     \item{glottal_flow_max}{Peak glottal flow amplitude}
#'     \item{glottal_flow_min}{Minimum glottal flow}
#'     \item{glottal_derivative_peak}{Maximum flow derivative (MFDR)}
#'     \item{NAQ}{Normalized Amplitude Quotient (requires GCI)}
#'     \item{QOQ}{Quasi-Open Quotient (requires GCI)}
#'     \item{H1_H2}{First two harmonics difference in dB (requires F0)}
#'     \item{HRF}{Harmonic Richness Factor}
#'     \item{PSP}{Parabolic Spectral Parameter}
#'   }
#'
#' @details
#' **Voice Quality Parameters:**
#'
#' \itemize{
#'   \item \bold{NAQ}: Measures vocal fold closure. Higher = more breathy.
#'         Normal: 0.05-0.20. Requires GCI.
#'   \item \bold{QOQ}: Open quotient measure. Normal: 0.4-0.7. Requires GCI.
#'   \item \bold{H1-H2}: Spectral tilt. Positive = breathy, negative = pressed.
#'         Requires F0.
#'   \item \bold{HRF}: High/low frequency energy ratio. Always computed.
#'   \item \bold{PSP}: Spectral curvature measure. Always computed.
#' }
#'
#' **Optimization:** Uses Numba JIT for NAQ/QOQ (2-3x speedup).
#' Check with covarep_info().
#'
#' @references
#' Alku, P. (1992). Glottal inverse filtering. Speech Communication.
#' Kane, J., & Gobl, C. (2013). Voice source measures. Journal of Voice.
#'
#' @seealso
#' \code{\link{trk_covarep_srh}} for F0 estimation,
#' \code{\link{trk_covarep_iaif}} for glottal waveforms
#'
#' @examples
#' \dontrun{
#' # Basic usage (HRF and PSP only)
#' vq <- lst_covarep_vq("vowel.wav")
#'
#' # With F0 for H1-H2
#' f0_data <- trk_covarep_srh("vowel.wav", toFile = FALSE)
#' f0_mean <- median(f0_data$`F0[Hz]`[f0_data$VUV == 1, 1])
#' vq <- lst_covarep_vq("vowel.wav", f0 = f0_mean)
#'
#' # Batch processing
#' vq_all <- lst_covarep_vq(c("a.wav", "e.wav", "i.wav"))
#' }
#'
#' @export
```

---

## Conclusion

**Recommendation:** ✅ **Implement `lst_covarep_vq()` immediately**

**Justification:**
1. ✅ Python function ready and optimized
2. ✅ Returns summary measures (perfect fit for `lst_`)
3. ✅ High scientific/clinical value
4. ✅ Complements existing track functions
5. ✅ Clear implementation path

**Estimated Implementation Time:** 1-2 days
- R wrapper: 4 hours
- Tests: 2 hours
- Documentation: 1 hour
- Validation: 1 hour

**Other Candidates:** None currently viable
- Remaining COVAREP functions are either:
  - Time-series outputs (already handled as `trk_`)
  - Utilities (not worth dedicated R functions)
  - Not yet implemented in Python module

**Next Steps:**
1. Implement `lst_covarep_vq()` in `R/covarep_vq.R`
2. Add test suite in `tests/testthat/test-covarep-vq.R`
3. Generate documentation
4. Test with real speech data
5. Consider auto-F0 integration for convenience

---

**Document Version:** 1.0
**Last Updated:** October 19, 2025
**Status:** Ready for implementation

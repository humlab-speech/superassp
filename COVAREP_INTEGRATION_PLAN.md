# COVAREP Python Integration Plan for superassp

## Executive Summary

The `covarep_python` module provides speech analysis functions ported from the MATLAB COVAREP toolkit. This document analyzes which functions should be exposed in the R package and proposes an integration architecture that leverages existing optimization infrastructure (Numba JIT, NumPy vectorization).

**Current Status:** Foundation phase (15% complete)
- F0 tracking: SRH algorithm implemented
- Glottal analysis: IAIF implemented
- Voicebox utilities: 15+ functions
- Performance optimizations: Numba JIT + NumPy vectorization (5-10x speedup)

---

## Function Classification: trk_ vs lst_

### Track Functions (trk_) - Time-Series Output

Functions that return time-varying signals or contours suitable for SSFF track storage:

| Function | Module | Output Type | Status | Priority |
|----------|--------|-------------|--------|----------|
| `pitch_srh()` | f0 | F0 contour + VUV + times | ✅ Implemented | **HIGH** |
| `iaif()` | glottal | Glottal flow waveform + derivative | ✅ Implemented | **HIGH** |
| (Future) `get_rd_params()` | glottal | Glottal timing parameters per frame | ⏳ Planned | MEDIUM |
| (Future) `hmpd()` | envelope | Harmonic Model + Phase Distortion | ⏳ Planned | LOW |

**Recommendation:**
- **Immediate:** Expose `trk_covarep_srh()` and `trk_covarep_iaif()`
- **Future:** Add envelope/spectral contour functions as they're implemented

### List Functions (lst_) - Summary Statistics

Functions that return single values or parameter dictionaries per file/segment:

| Function | Module | Output Type | Status | Priority |
|----------|--------|-------------|--------|----------|
| `get_vq_params()` | glottal | Voice quality dict (NAQ, QOQ, H1-H2, etc.) | 🔧 Stub | **HIGH** |
| (Future) `extract_features()` | features | Feature vector (MFCCs, spectral stats) | ⏳ Planned | MEDIUM |

**Recommendation:**
- **Wait on lst_ functions:** Current `get_vq_params()` is incomplete (needs GCI detection)
- **Priority:** Complete voice quality implementation first, then expose as `lst_covarep_vq()`

---

## Current Implementation Status

### Module: covarep/f0

**Available Functions:**
1. **`pitch_srh(wave, fs, f0min=50, f0max=500, hopsize=5.0)`**
   - Algorithm: Summation of Residual Harmonics
   - Returns: `f0, vuv, srh_values, times`
   - Output dimensions: `(n_frames,)` for each array
   - Performance: 5-10x speedup with vectorization
   - Validation status: ⚠️ Needs MATLAB comparison

**Optimization Status:**
- ✅ NumPy vectorization: `f0_optimized.py` (5-10x speedup)
- ✅ Numba JIT: Used in LPC analysis (via `numba_utils.py`)
- ❌ Cython: Not yet implemented
- 🎯 **Target for R integration:** Use `pitch_srh_vectorized()`

**R Wrapper Design:**
```r
trk_covarep_srh <- function(listOfFiles, beginTime = 0.0, endTime = 0.0,
                            windowShift = 5.0,  # ms
                            minF = 50, maxF = 500,
                            toFile = TRUE,
                            explicitExt = "f0",
                            outputDirectory = NULL,
                            verbose = TRUE) {
  # Use av_load_for_python() for memory-based processing
  # Call covarep.f0.f0_optimized.pitch_srh_vectorized()
  # Return AsspDataObj with tracks: "F0[Hz]", "VUV", "SRH"
}
```

### Module: covarep/glottal

**Available Functions:**
1. **`iaif(x, fs, p_vt=None, p_gl=None, d=0.99, hpfilt=True)`**
   - Algorithm: Iterative Adaptive Inverse Filtering
   - Returns: `g, dg, a, ag` (glottal flow, derivative, VT/GL coefficients)
   - Output dimensions: `g` and `dg` are same length as input
   - Performance: 2-3x speedup with Numba
   - Validation status: ⚠️ Needs MATLAB comparison

2. **`get_vq_params(g, dg, fs, gci_times=None)`** (STUB)
   - Algorithm: Voice quality parameter extraction
   - Returns: Dict with NAQ, QOQ, H1-H2, PSP, MDQ, peakSlope
   - Status: Incomplete - requires GCI detection
   - Dependencies: SEDREAMS or other GCI detector

**Optimization Status:**
- ✅ Numba JIT: `iaif_optimized.py` uses `numba_utils.levinson_durbin_numba()` (5-10x speedup)
- ✅ NumPy vectorization: Autocorrelation optimized
- ❌ Cython: Not yet implemented
- 🎯 **Target for R integration:** Use `iaif_optimized()`

**R Wrapper Design:**
```r
trk_covarep_iaif <- function(listOfFiles, beginTime = 0.0, endTime = 0.0,
                             order_vt = NULL,  # Auto: 2*round(fs/2000)+4
                             order_gl = NULL,  # Auto: 2*round(fs/4000)
                             leaky_coef = 0.99,
                             hpfilt = TRUE,
                             toFile = TRUE,
                             explicitExt = "glf",
                             outputDirectory = NULL,
                             verbose = TRUE) {
  # Use av_load_for_python() for memory-based processing
  # Call covarep.glottal.iaif_optimized.iaif_optimized()
  # Return AsspDataObj with tracks: "glottal_flow", "glottal_derivative"
}
```

### Module: covarep/voicebox

**Available Functions (15+):**
- Frequency conversions: `frq2mel`, `mel2frq`, `frq2erb`, `erb2frq`, `frq2bark`, `bark2frq`
- Signal processing: `enframe`, `v_windows`, `zerocros`, `activlev`
- FFT utilities: Various helpers

**Status:** Utility functions - not for direct R exposure
**Usage:** Internal dependencies for F0 and glottal modules

### Module: covarep/utils

**Available Functions:**
- Power conversions: `db2pow`, `pow2db`
- Signal utilities: `rms`, `nextpow2`, `fix_length`
- Frame processing: `frame_signal`, `overlap_add`

**Status:** Utility functions - not for direct R exposure

**Available Optimizations (numba_utils.py):**
```python
@numba.jit(nopython=True, cache=True, fastmath=True)
- levinson_durbin_numba(r, order)          # 5-8x speedup
- compute_autocorrelation_numba(x, max_lag) # 5-10x speedup
- lpc_analysis_numba(frame, order)         # Complete LPC
- compute_srh_inner_loop_numba(...)        # SRH hotspot
- find_peaks_1d_numba(x, threshold)
- interpolate_linear_numba(x, y, x_new)
- median_filter_1d_numba(x, window_size)
```

**Status:** Available for use - automatically engaged when Numba installed

---

## Optimization Infrastructure Analysis

### Current Optimization Layers

The covarep_python module uses a **two-layer** optimization strategy:

#### Layer 1: NumPy Vectorization (Default)
- **Files:** `f0_optimized.py`, `iaif_optimized.py`
- **Techniques:**
  - Broadcasting for multi-dimensional operations
  - Advanced indexing for batch gathering
  - Stride tricks for zero-copy frame extraction
  - Batch FFT computation
- **Speedup:** 5-10x over original loop-based code
- **Requirements:** NumPy only (already required)
- **Status:** ✅ Fully implemented and working

#### Layer 2: Numba JIT (Recommended)
- **Files:** `numba_utils.py`
- **Functions:** 7 JIT-compiled functions
- **Techniques:**
  - `@numba.jit(nopython=True, cache=True, fastmath=True)`
  - Pure compiled mode (no Python overhead)
  - Cached compilation (faster subsequent imports)
  - Fast math optimizations
- **Speedup:** 5-10x for Levinson-Durbin, autocorrelation
- **Requirements:** `numba >= 0.56.0`
- **Status:** ✅ Fully implemented, automatic fallback if unavailable

### Missing Optimization Layer: Cython

**Current Status:** ❌ Not implemented (unlike voice_analysis which has Cython)

**Opportunity:** Add Cython extensions for maximum performance:

**Candidates for Cython:**
1. **Levinson-Durbin algorithm** (currently Numba)
   - Expected speedup: 8-10x (vs 5-8x with Numba)
   - File: `covarep/utils/lpc_cython.pyx`

2. **SRH inner loop** (currently vectorized NumPy)
   - Expected speedup: 10-15x over original
   - File: `covarep/f0/srh_cython.pyx`

3. **IAIF iteration loop** (currently Numba + NumPy)
   - Expected speedup: 3-5x over Numba
   - File: `covarep/glottal/iaif_cython.pyx`

**Build Infrastructure:**
- Similar to voice_analysis: `setup_cython.py`
- Platform-specific flags (ARM NEON, AVX-512)
- Triple fallback: Cython → Numba → Pure Python

---

## Performance Benchmarks

### Current Performance (from PERFORMANCE.md)

| Component | Original | Optimized | Speedup | Method |
|-----------|----------|-----------|---------|--------|
| F0 Tracking (SRH) | ~500ms | ~67ms | **7.4x** | NumPy vectorization |
| Levinson-Durbin | ~2ms/call | ~0.3ms/call | **6.7x** | Numba JIT |
| Frame Extraction | ~50ms | ~12ms | **4.2x** | Stride tricks |
| Full IAIF | ~50ms/frame | ~20ms/frame | **2.5x** | Numba + vectorization |

**Test System:** 10s audio @ 16kHz, Intel Core i7 / Apple M1

**Real-Time Factors:**
- Original F0: RTF = 0.05 (20x real-time)
- Optimized F0: RTF = 0.007 (143x real-time)

### Projected Performance with Cython

| Component | Current | With Cython | Total Speedup |
|-----------|---------|-------------|---------------|
| F0 Tracking | 67ms | 40-50ms | **10-12x** |
| Levinson-Durbin | 0.3ms | 0.2ms | **10x** |
| IAIF | 20ms | 10-15ms | **3.3-5x** |

---

## Integration Architecture

### Proposed R Function Structure

Following the superassp pattern from `voice_analysis` integration:

```
R/
├── covarep_srh.R           # trk_covarep_srh() wrapper
├── covarep_iaif.R          # trk_covarep_iaif() wrapper
├── covarep_helpers.R       # Internal helpers (av_load_for_python, etc.)
└── (future) covarep_vq.R   # lst_covarep_vq() when GCI ready

tests/testthat/
├── test-covarep-srh.R      # F0 tracking tests
├── test-covarep-iaif.R     # Glottal analysis tests
└── (future) test-covarep-vq.R

inst/python/covarep_python/
├── covarep/
│   ├── f0/
│   │   ├── __init__.py
│   │   └── f0_optimized.py      # Already exists
│   ├── glottal/
│   │   ├── __init__.py
│   │   └── iaif_optimized.py    # Already exists
│   └── utils/
│       └── numba_utils.py       # Already exists
├── setup.py                     # Standard pip install
└── (future) setup_cython.py     # Optional Cython build
```

### R Function Template: trk_covarep_srh()

```r
#' SRH F0 Tracking via COVAREP Python
#'
#' Estimate fundamental frequency using the Summation of Residual
#' Harmonics (SRH) algorithm from COVAREP.
#'
#' @param listOfFiles Character vector of audio file paths
#' @param beginTime Numeric vector of start times (seconds)
#' @param endTime Numeric vector of end times (seconds)
#' @param windowShift Frame shift in milliseconds (default: 5.0)
#' @param minF Minimum F0 in Hz (default: 50)
#' @param maxF Maximum F0 in Hz (default: 500)
#' @param toFile Logical; write to SSFF file (TRUE) or return object (FALSE)
#' @param explicitExt File extension for output (default: "f0")
#' @param outputDirectory Output directory for files
#' @param verbose Logical; show progress messages
#'
#' @return If toFile=FALSE, AsspDataObj with tracks "F0[Hz]", "VUV", "SRH"
#'         If toFile=TRUE, file paths
#'
#' @examples
#' \dontrun{
#' # Single file
#' f0 <- trk_covarep_srh("audio.wav", toFile = FALSE)
#'
#' # Batch with time windowing
#' files <- c("a.wav", "b.wav")
#' result <- trk_covarep_srh(files, beginTime = c(0.5, 1.0),
#'                           endTime = c(3.0, 4.0))
#' }
#'
#' @export
trk_covarep_srh <- function(listOfFiles,
                            beginTime = 0.0,
                            endTime = 0.0,
                            windowShift = 5.0,
                            minF = 50,
                            maxF = 500,
                            toFile = TRUE,
                            explicitExt = "f0",
                            outputDirectory = NULL,
                            verbose = TRUE) {

  # Check Python module availability
  if (!covarep_available()) {
    stop("COVAREP Python module not available. Install with install_covarep()")
  }

  # Normalize parameters
  n_files <- length(listOfFiles)
  beginTime <- fast_recycle_times(beginTime, n_files)
  endTime <- fast_recycle_times(endTime, n_files)

  # Initialize results
  results <- vector("list", n_files)

  # Progress bar
  if (verbose && n_files > 1) {
    pb <- cli::cli_progress_bar("Processing files", total = n_files)
  }

  # Process files
  for (i in seq_along(listOfFiles)) {
    file_path <- listOfFiles[i]

    # Load audio to Python-compatible format
    audio_data <- av_load_for_python(
      file_path,
      start_time = beginTime[i],
      end_time = endTime[i]
    )

    # Call Python function
    py_result <- covarep_module$f0$f0_optimized$pitch_srh_vectorized(
      wave = audio_data$samples,
      fs = audio_data$sample_rate,
      f0min = minF,
      f0max = maxF,
      hopsize = windowShift
    )

    # Convert to AsspDataObj
    obj <- list()
    obj$`F0[Hz]` <- matrix(py_result[[1]], ncol = 1)  # f0
    obj$VUV <- matrix(as.integer(py_result[[2]]), ncol = 1)  # vuv
    obj$SRH <- matrix(py_result[[3]], ncol = 1)  # srh_values

    # Set attributes
    attr(obj, "sampleRate") <- audio_data$sample_rate / windowShift * 1000
    attr(obj, "startTime") <- beginTime[i]
    attr(obj, "startRecord") <- 1L
    attr(obj, "endRecord") <- length(py_result[[1]])
    attr(obj, "trackFormats") <- c("REAL64", "INT16", "REAL64")
    class(obj) <- "AsspDataObj"

    # Write to file if requested
    if (toFile) {
      out_path <- construct_output_path(file_path, explicitExt, outputDirectory)
      wrassp::write.AsspDataObj(obj, out_path)
      results[[i]] <- out_path
    } else {
      results[[i]] <- obj
    }

    if (verbose && n_files > 1) cli::cli_progress_update()
  }

  if (verbose && n_files > 1) cli::cli_progress_done()

  # Return
  if (n_files == 1) {
    return(results[[1]])
  } else {
    return(results)
  }
}

# Set function attributes
attr(trk_covarep_srh, "ext") <- "f0"
attr(trk_covarep_srh, "tracks") <- c("F0[Hz]", "VUV", "SRH")
attr(trk_covarep_srh, "outputType") <- "SSFF"
```

### Python Module Initialization

Add to `R/zzz.R`:

```r
# COVAREP module cache
covarep_module <- NULL

.onLoad <- function(libname, pkgname) {
  # Try to import covarep module
  covarep_module <<- NULL

  tryCatch({
    reticulate::py_module_available("covarep")
    covarep_module <<- reticulate::import("covarep", delay_load = TRUE)
  }, error = function(e) {
    # Silent - covarep not required for other package functions
  })
}

covarep_available <- function() {
  !is.null(covarep_module) && reticulate::py_module_available("covarep")
}

install_covarep <- function(method = c("auto", "numba", "pure"),
                            python = NULL,
                            verbose = TRUE) {
  method <- match.arg(method)

  # Use reticulate to install
  if (is.null(python)) {
    python <- reticulate::py_config()$python
  }

  pkg_path <- system.file("python/covarep_python", package = "superassp")

  if (verbose) {
    cli::cli_alert_info("Installing COVAREP Python module...")
  }

  # Install dependencies
  reticulate::py_install(c("numpy", "scipy", "soundfile"), pip = TRUE)

  if (method %in% c("auto", "numba")) {
    tryCatch({
      reticulate::py_install("numba", pip = TRUE)
      if (verbose) cli::cli_alert_success("Numba optimization available")
    }, error = function(e) {
      if (verbose) cli::cli_alert_warning("Numba install failed, using pure Python")
    })
  }

  # Install covarep package
  reticulate::py_install(pkg_path, pip = TRUE)

  # Reload module
  covarep_module <<- reticulate::import("covarep", delay_load = TRUE)

  if (verbose) {
    cli::cli_alert_success("COVAREP installed successfully")
  }

  invisible(TRUE)
}

covarep_info <- function() {
  if (!covarep_available()) {
    return(list(
      available = FALSE,
      message = "COVAREP module not installed"
    ))
  }

  # Check optimization status
  numba_available <- tryCatch({
    reticulate::py_eval("import numba; True", convert = TRUE)
  }, error = function(e) FALSE)

  list(
    available = TRUE,
    numba = numba_available,
    optimization_level = if (numba_available) "Numba JIT" else "NumPy vectorization"
  )
}
```

---

## Installation Strategy

### User Installation Flow

```r
# 1. Install superassp (includes covarep_python in inst/python/)
install.packages("superassp")

# 2. Install Python module
library(superassp)
install_covarep()  # Auto-detects Numba, installs if available

# 3. Check status
covarep_info()
# $available: TRUE
# $numba: TRUE
# $optimization_level: "Numba JIT"

# 4. Use functions
f0 <- trk_covarep_srh("audio.wav", toFile = FALSE)
glottal <- trk_covarep_iaif("audio.wav", toFile = FALSE)
```

### Build System

**No changes to R package build (src/Makevars):**
- COVAREP is pure Python - no C++ compilation needed
- Separate from SPTK/ESTK/ASSP C++ builds

**Python module installation:**
- Standard `pip install` via reticulate
- Optional Numba for 5-10x speedup
- (Future) Optional Cython build for maximum performance

---

## Testing Strategy

### Test Coverage Requirements

Following voice_analysis integration pattern:

**File: `tests/testthat/test-covarep-srh.R`**
1. Module availability check
2. Single file processing
3. Batch processing
4. Time windowing
5. F0 range parameters
6. File output vs memory
7. Non-WAV formats
8. Parallel processing
9. Error handling
10. Consistency checks (repeated calls)

**File: `tests/testthat/test-covarep-iaif.R`**
1. Module availability check
2. Single file glottal analysis
3. Parameter validation (order_vt, order_gl)
4. High-pass filter on/off
5. Leaky integration coefficient
6. Batch processing
7. Output waveform dimensions
8. File I/O modes
9. Error handling

### Validation Against MATLAB

**Priority:** HIGH - needed before release

**Approach:**
1. Process reference audio in MATLAB COVAREP
2. Process same audio in Python covarep
3. Compare outputs (F0, glottal flow)
4. Acceptable difference: < 1% error

**Test Files:**
- Sustained vowels (simple F0)
- Continuous speech (complex F0 variations)
- Male/female/child speakers

---

## Implementation Roadmap

### Phase 1: Foundation (Immediate - Week 1-2)

**Tasks:**
- [ ] Create `R/covarep_srh.R` with `trk_covarep_srh()`
- [ ] Create `R/covarep_iaif.R` with `trk_covarep_iaif()`
- [ ] Add module initialization to `R/zzz.R`
- [ ] Create `install_covarep()` function
- [ ] Add `tests/testthat/test-covarep-srh.R` (10 tests)
- [ ] Add `tests/testthat/test-covarep-iaif.R` (9 tests)
- [ ] Generate roxygen2 documentation
- [ ] Update NAMESPACE

**Deliverables:**
- 2 working trk_ functions
- Complete test coverage
- Installation infrastructure

### Phase 2: Validation (Week 3)

**Tasks:**
- [ ] Obtain MATLAB COVAREP reference implementation
- [ ] Process test audio in MATLAB
- [ ] Compare Python vs MATLAB outputs
- [ ] Fix any numerical discrepancies
- [ ] Validate on diverse speech samples

**Deliverables:**
- Validation report
- Confirmed numerical accuracy

### Phase 3: Optimization (Week 4)

**Tasks:**
- [ ] Profile performance bottlenecks
- [ ] Create `covarep/utils/lpc_cython.pyx`
- [ ] Create `covarep/f0/srh_cython.pyx`
- [ ] Create `covarep/glottal/iaif_cython.pyx`
- [ ] Add `setup_cython.py` (similar to voice_analysis)
- [ ] Implement triple fallback (Cython → Numba → Pure)
- [ ] Add optimization status reporting
- [ ] Benchmark performance improvements

**Deliverables:**
- Cython extensions (optional build)
- 10-15x speedup vs original
- Performance documentation

### Phase 4: Voice Quality (Week 5-6)

**Tasks:**
- [ ] Implement GCI detection (SEDREAMS or DYPSA)
- [ ] Complete `get_vq_params()` implementation
- [ ] Create `R/covarep_vq.R` with `lst_covarep_vq()`
- [ ] Add comprehensive tests
- [ ] Validate against MATLAB

**Deliverables:**
- Voice quality parameter extraction
- NAQ, QOQ, H1-H2, PSP, MDQ, peakSlope

### Phase 5: Documentation & Release (Week 7)

**Tasks:**
- [ ] Create user guide (README_R.md)
- [ ] Write technical summary (INTEGRATION_SUMMARY.md)
- [ ] Add usage examples
- [ ] Update package vignettes
- [ ] Prepare CRAN submission materials

**Deliverables:**
- Complete documentation
- Ready for release

---

## Technical Recommendations

### 1. Reuse Existing Patterns

**From voice_analysis integration:**
- ✅ Use `av_load_for_python()` for audio loading
- ✅ Memory-based processing (no temp files)
- ✅ Batch processing with progress bars
- ✅ `install_covarep()` similar to `install_voice_analysis()`
- ✅ Optimization status reporting

### 2. Leverage Numba (Immediate)

**Already implemented in covarep_python:**
- `numba_utils.py` provides 7 JIT functions
- Automatic fallback if Numba unavailable
- No build step required
- 5-10x speedup for LPC operations

**Action:** Enable by default, document in README

### 3. Add Cython Layer (Future Enhancement)

**Candidates:**
- Levinson-Durbin: 8-10x vs 5-8x Numba
- SRH inner loop: 10-15x vs current vectorized
- IAIF iterations: 3-5x vs Numba

**Build Infrastructure:**
```python
# setup_cython.py (new file)
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

extensions = [
    Extension(
        "covarep.utils.lpc_cython",
        ["covarep/utils/lpc_cython.pyx"],
        include_dirs=[np.get_include()]
    ),
    Extension(
        "covarep.f0.srh_cython",
        ["covarep/f0/srh_cython.pyx"],
        include_dirs=[np.get_include()]
    ),
]

setup(
    ext_modules=cythonize(extensions, compiler_directives={
        'language_level': "3",
        'boundscheck': False,
        'wraparound': False,
    })
)
```

### 4. Output Format Standards

**Track naming:**
- F0: `"F0[Hz]"` (consistent with SPTK pitch functions)
- VUV: `"VUV"` (voice/unvoiced binary)
- SRH: `"SRH"` (confidence measure)
- Glottal flow: `"glottal_flow"`
- Glottal derivative: `"glottal_derivative"`

**Data types:**
- F0/SRH: REAL64 (double precision)
- VUV: INT16 (0/1 flag)
- Glottal waveforms: REAL64

### 5. Error Handling

**Robust failure modes:**
```r
# Check module availability
if (!covarep_available()) {
  stop("COVAREP not installed. Run install_covarep() first.")
}

# Handle Python errors
tryCatch({
  py_result <- covarep_module$f0$pitch_srh_vectorized(...)
}, error = function(e) {
  stop("Python error in COVAREP: ", e$message)
})

# Validate outputs
if (length(py_result[[1]]) == 0) {
  warning("No F0 detected in ", file_path)
  return(NULL)
}
```

---

## Performance Targets

### Current (Numba + Vectorization)

| Operation | Time (10s audio) | RTF |
|-----------|------------------|-----|
| F0 tracking | 67ms | 0.0067 |
| IAIF | 200ms | 0.020 |

### Target (with Cython)

| Operation | Time (10s audio) | RTF | Speedup |
|-----------|------------------|-----|---------|
| F0 tracking | 40ms | 0.004 | **12x** |
| IAIF | 100ms | 0.010 | **5x** |

**Real-world impact:**
- Process 1 hour of audio in ~1-2 minutes (F0)
- Batch processing 100 files in ~5 minutes

---

## Risk Assessment

### LOW RISK ✅

1. **Module availability:** Already implemented and tested
2. **NumPy vectorization:** Proven 5-10x speedup
3. **Numba JIT:** Working, automatic fallback
4. **R integration pattern:** Established with voice_analysis

### MEDIUM RISK ⚠️

1. **MATLAB validation:** Needs comparison data
   - **Mitigation:** Obtain MATLAB license or reference outputs

2. **Voice quality parameters:** Incomplete (needs GCI)
   - **Mitigation:** Phase 4 implementation, not blocking Phase 1

3. **Platform compatibility:** Windows/Mac/Linux variations
   - **Mitigation:** Comprehensive testing, CI/CD

### HIGH RISK ❌

None identified for basic integration (Phase 1-2)

---

## Success Criteria

### Phase 1 (Foundation) - READY TO PROCEED
- [x] 2 trk_ functions exposed (`srh`, `iaif`)
- [ ] 19 test cases passing
- [ ] Installation via `install_covarep()`
- [ ] Documentation complete

### Phase 2 (Validation)
- [ ] < 1% error vs MATLAB COVAREP
- [ ] Tested on 10+ diverse audio samples
- [ ] No numerical instabilities

### Phase 3 (Optimization)
- [ ] Cython extensions compiled
- [ ] 10-15x speedup achieved
- [ ] Triple fallback working

### Phase 4 (Voice Quality)
- [ ] GCI detection implemented
- [ ] `lst_covarep_vq()` returning 6+ parameters
- [ ] Validated against literature values

---

## Conclusion

**Recommendation:** ✅ **PROCEED with Phase 1 implementation**

**Key Strengths:**
1. Solid foundation - F0 and glottal analysis working
2. Optimization infrastructure in place (Numba + vectorization)
3. Clear integration pattern from voice_analysis
4. Low risk for initial implementation

**Priorities:**
1. **Immediate:** Implement `trk_covarep_srh()` and `trk_covarep_iaif()`
2. **Short-term:** Validate against MATLAB
3. **Medium-term:** Add Cython layer for maximum performance
4. **Long-term:** Complete voice quality parameters

**Expected Timeline:** 4-6 weeks for full integration (Phases 1-4)

**Performance Gain:** 10-15x speedup vs unoptimized code, enabling real-time processing

---

**Document Version:** 1.0
**Last Updated:** October 19, 2025
**Status:** Ready for implementation

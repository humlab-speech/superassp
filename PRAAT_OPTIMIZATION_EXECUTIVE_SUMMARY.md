# Praat Optimization Exploration: Executive Summary

## Context

You requested exploration of direct linking to the Praat C library (bundled in Parselmouth submodule) to optimize Praat-based DSP procedures in the superassp package, eliminating Python dependencies while maintaining performance.

## What Was Done

### 1. Assessed Current Implementations (Completed)

**Original Python/Parselmouth** (`praat_intensity`):
- Uses: R → reticulate → Python → Parselmouth → Praat C code
- Performance: **5.68 ms** (fastest)
- Dependencies: Python + parselmouth + reticulate
- Status: Production-ready, well-tested

**Pure C++ Reimplementation** (`praat_intensity_cpp_wrapper`):
- Uses: R → Rcpp → Custom C++ implementation
- Performance: **19.52 ms** (3.4x slower than Python)
- Dependencies: Only Rcpp
- Status: Completed, tested, numerically equivalent (0.1 dB difference)
- Use case: Python-free deployments

### 2. Explored Direct Praat Linking (Analysis Completed)

Investigated linking directly to Praat C library in `src/Parselmouth/praat/`:

**Technical Findings:**
- Praat library has **599 C++ source files**
- Minimal intensity support needs **~50-80 files**
- Complex build system (CMake → R Makevars translation)
- Custom memory management (`Thing`, `autoThing`)
- External dependencies (fmt library)
- Platform-specific compilation flags

**Effort Estimate:**
- Initial setup: **40-150 hours**
- Per additional procedure: **4-8 hours**
- Ongoing maintenance: **8-16 hours per Parselmouth update**

**Expected Performance:**
- Likely **0-10% faster** than Python/Parselmouth
- Similar to Parselmouth (uses same code)
- Much faster than pure C++ reimplementation

## Key Finding

**Direct Praat linking offers marginal benefits but high costs:**
- ✅ Eliminates Python dependency
- ✅ Similar performance to Parselmouth
- ❌ **100+ hours development effort**
- ❌ Increased maintenance burden
- ❌ Platform-specific build issues
- ❌ Only 3-4 Praat procedures currently in package

## Recommendations

### Primary Recommendation: Keep Current Hybrid Approach

**For performance-critical operations:**
- ✅ Use Python/Parselmouth (fastest, proven)
- `praat_intensity`, `praat_formant_burg`, `praat_pitch`, etc.

**For Python-free deployments:**
- ✅ Use C++ reimplementations (acceptable performance)
- `praat_intensity_cpp_wrapper` (already implemented)

**Rationale:**
- Python/Parselmouth is 3.4x faster than pure C++
- Direct linking provides minimal additional benefit
- Development cost (100+ hours) not justified for 3-4 procedures
- Both options already available in package

### Quick Win: Optimize Python Integration (8 hours)

Instead of complex direct linking, optimize existing Python interface:

```r
# Pre-load Parselmouth at package startup
.onLoad <- function(libname, pkgname) {
  if (!exists(".pm_loaded")) {
    .parselmouth_env$pm <- reticulate::import("parselmouth")
    .pm_loaded <<- TRUE
  }
}

# Batch processing for multiple files
praat_intensity_batch <- function(files, ...) {
  pm <- get_parselmouth()  # Already loaded
  results <- lapply(files, function(f) process_with_pm(f, ...))
}
```

Expected benefits:
- 20-50% faster for batch operations
- Minimal development time
- No added complexity

### Defer Direct Linking Until...

Reconsider direct Praat linking when:
1. **Scale justifies effort**: Migrating 10+ Praat procedures
2. **Python becomes deployment blocker** for >50% of users  
3. **Parselmouth development stalls** or becomes unmaintained
4. **Funding available** for 100+ hour development effort

## Files Created

### Analysis Documents
1. **PRAAT_INTENSITY_CPP_ASSESSMENT.md** - C++ reimplementation results
2. **PRAAT_DIRECT_LINKING_STRATEGY.md** - Technical approach for direct linking
3. **PRAAT_DIRECT_LINKING_ANALYSIS.md** - Cost-benefit analysis
4. **PRAAT_OPTIMIZATION_EXECUTIVE_SUMMARY.md** - This document

### Implementation (Completed)
1. **src/praat_intensity.cpp** - Pure C++ intensity implementation
2. **R/ssff_cpp_praat_intensity.R** - R wrapper for C++ version
3. **R/sptk_helpers.R** - Added `create_intensity_asspobj()` helper

### Configuration Updates
1. **src/Makevars** - Added praat_intensity.cpp to build
2. **src/superassp_init.c** - Registered C++ function

## Decision Matrix

| Approach | Performance | Dependencies | Dev Effort | Maintenance | Recommendation |
|----------|-------------|--------------|------------|-------------|----------------|
| **Python/Parselmouth** | ★★★★★ (5.7ms) | Python + reticulate | None (exists) | Low | **✅ DEFAULT** |
| **C++ Reimplementation** | ★★☆☆☆ (19.5ms) | Rcpp only | Done | Low | **✅ FALLBACK** |
| **Direct Praat Linking** | ★★★★★ (6ms est) | Rcpp only | 100+ hrs | High | **❌ NOT NOW** |
| **Optimized Python** | ★★★★★ (4-5ms est) | Python + reticulate | 8 hrs | Low | **✅ QUICK WIN** |

## Summary

The exploration successfully identified three approaches for Praat DSP optimization:

1. **Existing Python/Parselmouth**: Best performance, acceptable dependencies
2. **New C++ Reimplementation**: Good for Python-free deployments, acceptable performance
3. **Direct Praat Linking**: Technically feasible but economically unjustified

**Conclusion**: The package now has both fast (Python) and dependency-free (C++) options. Direct Praat linking should be deferred until scale or requirements change. A small investment in Python optimization offers better ROI.

---

**Date**: 2025-10-18  
**Status**: Analysis Complete - Implementation Deferred  
**Next Action**: Implement Python optimization (optional, 8 hours)  
**Review Date**: When migrating 10+ Praat procedures or Python becomes blocker

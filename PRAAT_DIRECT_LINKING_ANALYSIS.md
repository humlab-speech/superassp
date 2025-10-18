# Praat Direct Linking: Implementation Analysis

## Current Situation Assessment

### Parselmouth Structure
- **599 Praat C++ source files** in `src/Parselmouth/praat/`
- Organized into modules: `melder`, `sys`, `stat`, `fon`, `dwtools`, etc.
- Built as single static library `libpraat` via CMake
- Uses `fmt` library for formatting (external dependency)
- Requires many defines: `NO_GUI`, `NO_AUDIO`, `NO_GRAPHICS`, `NO_NETWORK`

### Dependencies for Intensity Analysis

Based on `fon/CMakeLists.txt`, minimal set needed:

**Core Objects (fon/):**
- `Function.cpp`, `Sampled.cpp`, `Matrix.cpp`, `Vector.cpp`
- `Sound.cpp`, `Sound_files.cpp` 
- `Intensity.cpp`, `Sound_to_Intensity.cpp`

**Numerical (melder/):**
- `melder.cpp`, `melder_alloc.cpp`, `melder_str32.cpp`
- `NUMmath.cpp`, `NUMspecfunc.cpp` (Bessel functions)
- `NUM.cpp`, `VEC.cpp`, `MAT.cpp`

**System (sys/):**
- `Thing.cpp`, `Data.cpp`, `Collection.cpp`
- Memory management, error handling

**Total estimate**: ~50-80 source files for minimal intensity support

## Complexity Analysis

### High Complexity Factors

1. **Build System Integration**
   - CMake → R Makevars translation non-trivial
   - Platform-specific flags and defines
   - Dependency management (fmt library)

2. **Praat's Object System**
   - Custom memory management (`Thing`, `autoThing`)
   - Reference counting and ownership
   - Exception handling via longjmp/setjmp

3. **Header Dependencies**
   - Complex include hierarchy
   - Many cross-dependencies between modules
   - Template-heavy code in newer Praat versions

4. **Maintenance Burden**
   - Need to track Praat/Parselmouth updates
   - Platform-specific issues (Windows, Mac, Linux)
   - Potential conflicts with R's memory management

### Low Complexity Alternative

**Keep using Parselmouth via Python**, but optimize the interface:
- Pre-load Python/Parselmouth in package namespace
- Batch processing to amortize Python overhead
- Convert results once, process many analyses

## Recommended Approach: Hybrid Strategy

### Tier 1: Critical Functions (Python/Parselmouth)
Functions where performance matters most:
- ✅ Keep `praat_intensity` (Python) - already 3.4x faster than C++
- ✅ Keep `praat_formant_burg` (Python) - optimized
- ✅ Keep `praat_pitch` (Python) - complex, well-tested

### Tier 2: Self-Contained Functions (Direct C++)
Functions that can be reimplemented without massive dependencies:
- ✅ `praat_intensity_cpp` - already done, good for Python-free deployments
- Potential: Simple spectral analyses
- Potential: Basic filtering operations

### Tier 3: Future Direct Linking (When Justified)
Consider direct Praat linking when:
- ❌ NOT for single functions - overhead too high
- ✅ When migrating 10+ Praat procedures
- ✅ When Python becomes a blocker for deployment
- ✅ When we need customizations not available in Parselmouth

## Practical Implementation Path (If Proceeding)

### Phase 1: Proof of Concept (1-2 weeks)
1. Extract minimal Praat subset for intensity (estimate: 50 files)
2. Create standalone Makefile for `libpraat_minimal.a`
3. Test compilation on target platforms
4. Implement single wrapper function
5. Benchmark vs Python/Parselmouth

**Go/No-Go Decision Point**: If PoC takes >2 weeks or performance isn't better than Python, stop here.

### Phase 2: Production Integration (1-2 weeks)
1. Integrate into R package Makevars
2. Add proper error handling and memory management
3. Create R wrapper functions
4. Add comprehensive tests
5. Document for future procedures

### Phase 3: Migration (Ongoing)
- Migrate procedures one-by-one
- Keep Python fallbacks
- Monitor performance and reliability

## Effort Estimate

### Initial Setup
- **Minimal (intensity only)**: 40-60 hours
- **Full Praat library**: 100-150 hours
- **Per additional procedure**: 4-8 hours

### Ongoing Maintenance
- **Parselmouth updates**: 8-16 hours per update
- **Platform-specific issues**: 4-8 hours per platform per issue
- **Bug fixes**: Variable

## Risk Assessment

| Risk | Probability | Impact | Mitigation |
|------|-------------|--------|------------|
| Build fails on Windows | High | High | Extensive testing, fallback to Python |
| Memory leaks | Medium | High | Valgrind testing, RAII wrappers |
| Performance worse than Python | Low | High | Benchmark early, abort if true |
| Maintenance burden | High | Medium | Good documentation, modular design |
| Praat API changes | Low | Medium | Pin to specific version, test updates |

## Cost-Benefit Analysis

### Costs
- **Development**: 40-150 hours initial + ongoing maintenance
- **Testing**: Increased test surface across platforms
- **Documentation**: Need comprehensive guides for future devs
- **Risk**: Potential instability during development

### Benefits
- **Eliminates Python dependency**: Yes (but Python still needed for other features)
- **Performance improvement**: Likely 0-10% vs Parselmouth (not 3x as hoped)
- **Code reuse**: Yes, can add more Praat procedures
- **Maintainability**: Decreases (more complex build system)

### ROI Assessment
**Break-even point**: Need to migrate ~20 Praat procedures to justify initial investment.

**Current situation**: Only 3-4 procedures using Praat/Parselmouth.

**Recommendation**: **NOT justified at this time**. Revisit when:
1. Python becomes deployment blocker for users
2. Need to migrate 10+ Praat procedures
3. Performance of Python/Parselmouth becomes unacceptable

## Alternative: Optimize Python Integration

Instead of direct linking, optimize the existing Python approach:

### Quick Wins (4-8 hours)
1. **Pre-load Parselmouth**: Load once at package startup
2. **Batch processing**: Process multiple files per Python call
3. **Result caching**: Cache converted results
4. **Lazy loading**: Only load Python when Praat functions called

### Example Optimization
```r
# Package startup (.onLoad)
.parselmouth_env <- new.env()
.parselmouth_env$loaded <- FALSE
.parselmouth_env$pm <- NULL

load_parselmouth <- function() {
  if (!.parselmouth_env$loaded) {
    reticulate::use_virtualenv("r-superassp", required = FALSE)
    .parselmouth_env$pm <- reticulate::import("parselmouth")
    .parselmouth_env$loaded <- TRUE
  }
  .parselmouth_env$pm
}

# Batch processing wrapper
praat_intensity_batch <- function(files, ...) {
  pm <- load_parselmouth()
  # Process all files in Python, return once
  results <- lapply(files, function(f) {
    snd <- pm$Sound(f)
    intensity <- snd$to_intensity(...)
    # Convert and return
  })
}
```

Expected speedup: 20-50% for batch operations with negligible development cost.

## Final Recommendation

### For Current Project
**Do NOT implement direct Praat linking yet.**

**Instead:**
1. ✅ Keep Python/Parselmouth for production (best performance)
2. ✅ Keep C++ reimplementations for Python-free deployments
3. ✅ Document this analysis for future reference
4. ✅ Implement Python optimization (lazy loading, batch processing)
5. ⏸️ Defer direct linking until justified by scale

### Criteria for Reconsidering
- Migrating 10+ Praat procedures
- Python deployment becomes blocker for >50% of users
- Parselmouth development stalls
- Funding available for 100+ hour development effort

---

**Status**: Recommend NOT implementing at this time  
**Rationale**: Cost (100+ hours) exceeds benefit (marginal improvement over Python)  
**Alternative**: Optimize Python integration (8 hours), keep for future  
**Date**: 2025-10-18

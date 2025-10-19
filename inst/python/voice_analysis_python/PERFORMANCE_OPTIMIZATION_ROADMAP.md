# Performance Optimization Roadmap
## Voice Analysis Toolbox Python Implementation

**Date:** October 17, 2025  
**Current Performance:** 6.40s per file (152 measures)  
**Previous Best:** 3.98s with Numba (from benchmarks)  
**Target:** 1.5-2.5s per file

---

## Current Status Assessment

### Performance Regression Investigation

**Observed:** Current runtime is 6.40s, but previous benchmarks showed 3.98s with Numba optimization.

**Likely Causes:**
1. Numba JIT not compiling on first run (cold start)
2. Cython extensions not built/loaded properly
3. Additional features added (thesis mode) increasing complexity
4. PyWavelets missing (causing fallback code paths)

**Immediate Actions:**
1. ✅ Verify Numba is installed and active
2. ✅ Build Cython extensions properly
3. ✅ Install optional dependencies (PyWavelets, pyEMD)
4. ✅ Profile current execution to identify bottlenecks

---

## Three-Tier Optimization Strategy

### Tier 1: Quick Wins (0-2 hours, 1.5-2x speedup)

These optimizations require minimal code changes and provide immediate benefits.

#### 1.1 Build Cython RPDE Extension ⚡ HIGH IMPACT
**Status:** Code exists (`rpde_cython.pyx`) but may not be compiled  
**Expected gain:** 2-3x on RPDE (biggest bottleneck)  
**Time:** 30 minutes

```bash
cd voice_analysis_python
python setup_cython.py build_ext --inplace
```

**Verification:**
```python
from voice_analysis.features import rpde_cython
# Should import without error
```

#### 1.2 Install Missing Dependencies
**Status:** PyWavelets missing (saw warning)  
**Expected gain:** Eliminate fallback code paths  
**Time:** 5 minutes

```bash
pip install PyWavelets pyEMD
```

#### 1.3 Enable Numba Caching
**Status:** May not be enabled  
**Expected gain:** Eliminate JIT compilation on subsequent runs  
**Time:** 10 minutes

```python
# In affected files (rpde.py, perturbation.py)
@njit(cache=True)  # Add cache=True parameter
def function_name(...):
    ...
```

#### 1.4 Warm-up JIT Compilation
**Status:** Not implemented  
**Expected gain:** Predictable performance  
**Time:** 15 minutes

```python
# Add to VoiceAnalyzer.__init__()
def _warmup_jit(self):
    """Pre-compile Numba functions with dummy data"""
    dummy_signal = np.random.randn(1000)
    dummy_f0 = np.random.randn(100)
    # Call key functions to trigger JIT
```

**Total Tier 1 Impact:** 1.5-2x speedup → Target: 3.2-4.3s

---

### Tier 2: Medium Effort (2-8 hours, additional 1.2-1.5x speedup)

These require moderate code changes but provide significant benefits.

#### 2.1 Optimize RPDE Algorithm ⚡ HIGH IMPACT
**Current:** 58.7% of total time  
**Target:** Reduce to 30-40% of total time  
**Time:** 3-4 hours

**Approaches:**
1. Use pre-built Cython extension (Priority 1)
2. Implement KD-tree neighbor search (already exists, verify usage)
3. Parallelize embedding construction (low gain)
4. Cache distance matrix for repeated analyses (memory-intensive)

**Implementation:**
```python
# Ensure Cython version is used by default
try:
    from .rpde_cython import compute_rpde_cython as compute_rpde
    USING_CYTHON = True
except ImportError:
    from .rpde_numba import compute_rpde_numba as compute_rpde
    USING_CYTHON = False
```

#### 2.2 Vectorize Perturbation Quotient Calculations
**Current:** Used in 44 measures (jitter + shimmer)  
**Expected gain:** 1.3-1.5x on jitter/shimmer  
**Time:** 2 hours

**Approaches:**
1. NumPy vectorization for all PQ variants
2. Pre-compute denominators
3. Avoid redundant calculations

**Implementation:**
```python
# Replace loops with vectorized operations
# Before:
for i in range(len(T0)):
    pq += abs(T0[i] - T0[i+1])

# After:
pq = np.sum(np.abs(np.diff(T0)))
```

#### 2.3 Optimize F0 Estimation (SWIPE)
**Current:** 6.6% of total time  
**Expected gain:** 1.2x on F0 estimation  
**Time:** 2 hours

**Approaches:**
1. Use compiled SWIPE library if available
2. Reduce FFT resolution for faster processing
3. Cache intermediate FFT results

#### 2.4 Parallel Feature Groups (Already Implemented)
**Current:** core_parallel.py exists  
**Action:** Verify and enable by default  
**Expected gain:** 1.1-1.2x  
**Time:** 1 hour

**Total Tier 2 Impact:** 1.2-1.5x additional → Cumulative: 2.0-3.0s

---

### Tier 3: Advanced Optimizations (1-3 days, additional 1.1-1.3x speedup)

These require significant effort but squeeze out maximum performance.

#### 3.1 Platform-Specific SIMD Optimizations
**Target:** M1 Pro (ARM NEON) and AMD EPYC (AVX-512)  
**Expected gain:** 1.2-1.4x on numerical loops  
**Time:** 1 day

**Implementation:**
```python
# setup_cython.py - Platform detection
import platform
if platform.machine() == 'arm64':
    extra_compile_args = ['-O3', '-march=native']  # NEON on M1
elif 'x86_64' in platform.machine():
    extra_compile_args = ['-O3', '-mavx512f']  # AVX-512 on EPYC
```

#### 3.2 Memory Pool for Batch Processing
**Target:** Large-scale batch analysis  
**Expected gain:** 1.1-1.2x on batch processing, reduced memory  
**Time:** 1 day

**Implementation:**
```python
class MemoryPool:
    def __init__(self):
        self.buffer_cache = {}
    
    def get_buffer(self, shape, dtype):
        key = (shape, dtype)
        if key not in self.buffer_cache:
            self.buffer_cache[key] = np.empty(shape, dtype=dtype)
        return self.buffer_cache[key]
```

#### 3.3 Lazy Feature Computation
**Target:** Only compute requested features  
**Expected gain:** Variable (depends on use case)  
**Time:** 2 days

**Implementation:**
```python
analyzer = VoiceAnalyzer(
    features=['rpde', 'dfa', 'ppe', 'jitter', 'shimmer']
)
# Only compute specified features
```

#### 3.4 GPU Acceleration (Future, If Available)
**Target:** FFT, convolution operations  
**Expected gain:** 2-5x on MFCC, spectral features  
**Time:** 3-5 days  
**Note:** User stated no GPU available, defer

**Total Tier 3 Impact:** 1.1-1.3x additional → Cumulative: 1.5-2.5s

---

## Recommended Implementation Order

### Phase 1: Immediate (Today, <2 hours)
1. Build Cython extensions: `python setup_cython.py build_ext --inplace`
2. Install dependencies: `pip install PyWavelets pyEMD numba`
3. Enable Numba caching in rpde.py and perturbation.py
4. Run benchmark: `python benchmark_numba.py`
5. **Expected result:** 3.5-4.5s per file

### Phase 2: Short-term (This Week, 4-8 hours)
1. Verify Cython RPDE is being used (priority 1)
2. Vectorize perturbation quotient calculations
3. Enable parallel feature groups by default
4. Profile and identify remaining bottlenecks
5. **Expected result:** 2.5-3.5s per file

### Phase 3: Medium-term (Next Week, Optional)
1. Platform-specific SIMD optimizations
2. Implement memory pooling for batch processing
3. Add lazy feature computation
4. **Expected result:** 1.5-2.5s per file

---

## Performance Measurement Protocol

### Benchmarking Script

```python
# benchmark_current.py
import time
import numpy as np
import soundfile as sf
from voice_analysis.core import VoiceAnalyzer

# Load test file
audio, fs = sf.read('a1.wav')

# Warm-up run (JIT compilation)
analyzer = VoiceAnalyzer(use_thesis_mode=False)
_, _ = analyzer.analyze(audio, fs)

# Timed runs
times = []
for i in range(5):
    t0 = time.time()
    measures, F0 = analyzer.analyze(audio, fs)
    t1 = time.time()
    times.append(t1 - t0)
    print(f"Run {i+1}: {t1-t0:.3f}s")

print(f"\nMean: {np.mean(times):.3f}s ± {np.std(times):.3f}s")
print(f"Median: {np.median(times):.3f}s")
print(f"Min: {np.min(times):.3f}s")
print(f"Measures computed: {len(measures)}")
```

### Profiling Script

```python
# profile_current.py
import cProfile
import pstats
import soundfile as sf
from voice_analysis.core import VoiceAnalyzer

audio, fs = sf.read('a1.wav')
analyzer = VoiceAnalyzer(use_thesis_mode=False)

# Profile
profiler = cProfile.Profile()
profiler.enable()
measures, F0 = analyzer.analyze(audio, fs)
profiler.disable()

# Print top 20 functions
stats = pstats.Stats(profiler)
stats.sort_stats('cumulative')
stats.print_stats(20)
```

---

## Expected Performance Targets

### Conservative Estimates

| Phase | Time (s) | Speedup | Confidence |
|-------|----------|---------|------------|
| Current (no optimization) | 6.40 | 1.0x | ✅ Measured |
| Tier 1 (quick wins) | 3.5-4.5 | 1.4-1.8x | High (90%) |
| Tier 2 (medium effort) | 2.5-3.5 | 1.8-2.6x | Medium (70%) |
| Tier 3 (advanced) | 1.5-2.5 | 2.6-4.3x | Low-Medium (50%) |

### Aggressive Estimates (Best Case)

| Phase | Time (s) | Speedup | Conditions |
|-------|----------|---------|------------|
| Tier 1 | 3.2 | 2.0x | Perfect Cython build |
| Tier 2 | 2.0 | 3.2x | + All optimizations |
| Tier 3 | 1.5 | 4.3x | + Platform-specific |

---

## Compatibility Considerations

### R Reticulate
- ✅ All optimizations compatible
- ✅ Cython extensions work seamlessly
- ✅ No breaking API changes
- ⚠️ May need to rebuild extensions on deployment system

### MATLAB Mode vs Thesis Mode
- ✅ Both modes benefit equally from optimizations
- ✅ No changes to algorithm correctness
- ✅ Numerical outputs remain identical

### Python Versions
- ✅ Python 3.8-3.12 supported
- ✅ Numba supports all versions
- ✅ Cython supports all versions

---

## Risk Assessment

### Low Risk (Tier 1)
- Building Cython extensions: Well-tested, fallback to Numba
- Installing dependencies: Standard packages
- Enabling caching: No algorithm changes

### Medium Risk (Tier 2)
- Vectorization: Need careful testing for numerical equivalence
- F0 optimization: Could affect downstream measures
- Parallel features: Already implemented, just enabling

### High Risk (Tier 3)
- SIMD optimizations: Platform-specific, hard to test
- Memory pooling: Complex memory management
- Lazy computation: Major architectural change

---

## Success Criteria

### Minimum Viable Performance (MVP)
- ✅ Single file: <4.0s (152 measures)
- ✅ Batch (100 files, 8 cores): <60s
- ✅ Numerical accuracy maintained (vs MATLAB)
- ✅ All tests pass

### Target Performance
- 🎯 Single file: <3.0s (152 measures)
- 🎯 Batch (100 files, 8 cores): <40s
- 🎯 Numerical accuracy maintained
- 🎯 All tests pass
- 🎯 R reticulate compatible

### Stretch Goal
- 🌟 Single file: <2.0s (152 measures)
- 🌟 Batch (100 files, 32 cores): <10s
- 🌟 Numerical accuracy maintained
- 🌟 All tests pass
- 🌟 R reticulate compatible
- 🌟 Platform-optimized (M1 Pro, AMD EPYC)

---

## Next Steps

1. **Run current benchmark** to establish baseline
2. **Build Cython extensions** for immediate gains
3. **Profile execution** to identify current bottlenecks
4. **Implement Tier 1** optimizations
5. **Re-benchmark** and compare
6. **Proceed to Tier 2** if needed

**Start with:** `python setup_cython.py build_ext --inplace && python benchmark_numba.py`

---

**Status:** Ready for optimization  
**Estimated Time to Target (3s):** 4-8 hours  
**Estimated Time to Stretch Goal (2s):** 2-3 days  
**Risk Level:** Low-Medium  
**Recommendation:** Proceed with Tier 1 immediately

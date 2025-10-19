# Performance Analysis Summary

**Date:** October 17, 2025  
**Current Performance:** 4.92-5.07s per file (219 measures)  
**Status:** ⚠️ Slower than expected, optimization needed

---

## Key Findings from Profiling

### Bottleneck Analysis

From profiler output, total time is 6.542s with the following breakdown:

| Component | Time (s) | % Total | Status |
|-----------|----------|---------|--------|
| **RPDE** | 3.112 | 47.6% | ❌ Using Cython but slow |
| **DYPSA** (2 calls) | 1.767 | 27.0% | ❌ Major bottleneck |
| **EMD Features** | 1.014 | 15.5% | ⚠️ PyEMD installed but slow |
| **GNE** | 0.358 | 5.5% | ✓ Acceptable |
| **HNR/NHR** | 0.169 | 2.6% | ✓ Good |
| Other | 0.122 | 1.8% | ✓ Minimal |

### Critical Issues Identified

1. **RPDE Taking 3.1s** - Should be ~0.5s with Cython
   - Cython is available and loaded
   - Performance suggests it's not being used effectively
   - Need to verify Cython path and compilation

2. **DYPSA Taking 1.7s** - Used by GQ and VFER
   - Called twice (once for GQ, once for VFER)
   - `_compute_mean_based_signal` is the bottleneck (1.688s)
   - This is pure Python, no optimization

3. **EMD Taking 1.0s** - PyEMD is slow
   - Using PyEMD package (installed)
   - Known to be slow for Python
   - 6 features computed

4. **Extra Features Computed** - 219 measures vs expected 152
   - Includes all optional features (EMD, Wavelet, etc.)
   - May not all be needed

---

## Root Cause Analysis

### Why is RPDE so slow?

The Cython extension is built (`rpde_cython.cpython-312-darwin.so` exists) but RPDE is taking 3.1s instead of expected ~0.5s.

**Possible causes:**
1. Cython not being called (check use_cython=True)
2. Cython compiled without optimizations
3. Input size too large causing memory issues
4. Signal resampling overhead

**Investigation needed:**
```python
# Check if Cython is actually being used
from voice_analysis.features.rpde import CYTHON_AVAILABLE
print(f"Cython available: {CYTHON_AVAILABLE}")

# Time Cython vs Numba directly
import time
import numpy as np
audio = np.random.randn(100000)

# Force Cython
t0 = time.time()
from voice_analysis.features.rpde import compute_rpde
result = compute_rpde(audio, fs=25000, use_cython=True, use_kdtree=False)
t1 = time.time()
print(f"Cython: {t1-t0:.3f}s")

# Force Numba
t0 = time.time()
result = compute_rpde(audio, fs=25000, use_cython=False, use_kdtree=False)
t1 = time.time()
print(f"Numba: {t1-t0:.3f}s")
```

### Why is DYPSA so slow?

DYPSA (Dynamic Programming Phase Slope Algorithm) is used to detect glottal closure instants. It's called twice per analysis.

**Issues:**
- Pure Python implementation (no Numba/Cython)
- Complex algorithm (dynamic programming)
- Called from both `gq.py` and `vfer.py`

**Solutions:**
1. Cache DYPSA results (call once, use twice)
2. Optimize with Numba
3. Consider if GQ/VFER are essential (they return NaN for most signals)

---

## Immediate Action Plan

### Priority 1: Fix RPDE Performance (Target: 3.1s → 0.5s, save 2.6s)

**Actions:**
1. Verify Cython is being called
2. Rebuild Cython with maximum optimization
3. Profile Cython code to find bottleneck
4. Consider if resampling to 25kHz is needed

**Expected impact:** 40% total time reduction

### Priority 2: Optimize or Cache DYPSA (Target: 1.7s → 0.2s, save 1.5s)

**Options:**
A. **Cache DYPSA results** (Quick, 5 minutes)
   ```python
   # In core.py, call DYPSA once
   gci = dypsa(audio, fs)  # Call once
   
   # Pass to both features
   gq_features = compute_glottal_quotient(audio, fs, F0, gci=gci)
   vfer_features = compute_vfer(audio, fs, F0, gci=gci)
   ```

B. **Optimize with Numba** (Medium, 2 hours)
   - Add @njit decorators to hot loops
   - Vectorize where possible

C. **Make GQ/VFER optional** (Quick, 10 minutes)
   - Most signals return NaN (insufficient glottal instants)
   - Could skip if not needed

**Expected impact:** 23% total time reduction

### Priority 3: Make EMD Optional or Faster (Target: 1.0s → 0s, save 1.0s)

**Options:**
A. **Make EMD optional** (Quick)
   - Add parameter to disable EMD features
   - Reduces measure count to 213

B. **Find faster EMD implementation**
   - Current PyEMD is slow
   - Check for C++ bindings

**Expected impact:** 15% total time reduction

---

## Performance Targets After Fixes

| Scenario | RPDE | DYPSA | EMD | Other | **Total** | Measures |
|----------|------|-------|-----|-------|-----------|----------|
| **Current** | 3.1s | 1.7s | 1.0s | 1.2s | **6.9s** | 219 |
| After P1 (RPDE fix) | 0.5s | 1.7s | 1.0s | 1.2s | **4.4s** | 219 |
| After P1+P2 (+ DYPSA) | 0.5s | 0.2s | 1.0s | 1.2s | **2.9s** | 219 |
| After P1+P2+P3 (+ EMD off) | 0.5s | 0.2s | 0.0s | 1.2s | **1.9s** | 213 |
| **MATLAB mode (152)** | 0.5s | 0.0s | 0.0s | 1.0s | **1.5s** | 152 |

---

## Recommended Implementation Steps

### Step 1: Investigate RPDE (15 minutes)
```bash
cd voice_analysis_python
python -c "
import numpy as np
import time
from voice_analysis.features.rpde import compute_rpde, CYTHON_AVAILABLE

print(f'Cython available: {CYTHON_AVAILABLE}')
audio = np.random.randn(100000)

# Test Cython
t0 = time.time()
r = compute_rpde(audio, fs=25000, use_cython=True, use_kdtree=False)
t1 = time.time()
print(f'Cython: {t1-t0:.3f}s, result={r:.6f}')

# Test Numba
t0 = time.time()
r = compute_rpde(audio, fs=25000, use_cython=False, use_kdtree=False)
t1 = time.time()
print(f'Numba: {t1-t0:.3f}s, result={r:.6f}')
"
```

### Step 2: Cache DYPSA Results (10 minutes)

Edit `voice_analysis/core.py`:
```python
# Around line 120, after F0 extraction
# Call DYPSA once if needed
gci = None
if True:  # Always compute for GQ/VFER
    from .utils.dypsa import dypsa
    gci = dypsa(audio, fs)

# Pass to features
gq_features = compute_glottal_quotient(audio, fs, F0, voicing, gci=gci)
vfer_features = compute_vfer(audio, fs, F0, voicing, gci=gci)
```

Edit `gq.py` and `vfer.py` to accept optional `gci` parameter.

### Step 3: Add Feature Selection (20 minutes)

Edit `VoiceAnalyzer.__init__`:
```python
def __init__(self, f0_min=50, f0_max=500, f0_algorithm='SWIPE', 
             use_thesis_mode=False, enable_emd=False, enable_glottal=False):
    ...
    self.enable_emd = enable_emd
    self.enable_glottal = enable_glottal
```

Then in `analyze()`, conditionally compute:
```python
if self.enable_emd:
    emd_features = compute_emd_features(audio, fs, F0)
    measures.update(emd_features)

if self.enable_glottal:
    gq_features = compute_glottal_quotient(...)
    vfer_features = compute_vfer(...)
    measures.update(gq_features)
    measures.update(vfer_features)
```

### Step 4: Rebuild Cython with Optimizations (5 minutes)

Edit `setup_cython.py` to ensure optimization flags:
```python
extensions = [
    Extension(
        "voice_analysis.features.rpde_cython",
        ["voice_analysis/features/rpde_cython.pyx"],
        extra_compile_args=['-O3', '-march=native', '-ffast-math'],
        extra_link_args=['-O3'],
    ),
]
```

Rebuild:
```bash
python setup_cython.py build_ext --inplace --force
```

---

## Expected Results

After all optimizations:

- **MATLAB mode (152 measures):** 1.5-2.0s per file
- **Full mode with EMD (219 measures):** 2.5-3.0s per file  
- **Full mode without EMD (213 measures):** 2.0-2.5s per file

This meets the target of <3s for comprehensive analysis.

---

## Next Action

Run the RPDE investigation script to determine if Cython is being used effectively, then proceed with the appropriate fix.

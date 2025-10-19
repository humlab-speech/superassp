# Performance Optimization Implementation Guide
## Voice Analysis Toolbox - Quick Wins

**Date:** October 17, 2025  
**Status:** Ready to implement  
**Expected time:** 30-60 minutes  
**Expected improvement:** 3-4x speedup (6.5s → 1.5-2.0s)

---

## Executive Summary

Performance profiling revealed three major bottlenecks:

1. **RPDE** (47.6% of time) - Using wrong algorithm, KD-tree is 7x faster
2. **DYPSA** (27.0% of time) - Called twice, can cache results
3. **EMD** (15.5% of time) - Optional feature, can disable by default

**Simple fixes available that require minimal code changes.**

---

## Fix 1: Enable KD-Tree for RPDE ✅ COMPLETE

**Status:** ✅ Already implemented in previous step

**What was done:**
- Changed default from `use_cython=True` to `use_cython=False`
- Changed KD-tree threshold from `M > 1000` to `M > 100`
- Made KD-tree the preferred method

**Expected impact:** 3.1s → 0.4s (2.7s saved, 41% reduction)

**Verification:**
```bash
cd voice_analysis_python
python -c "
import numpy as np
import time
import soundfile as sf
from voice_analysis.features.rpde import compute_rpde

audio, fs = sf.read('../a1.wav')
t0 = time.time()
r = compute_rpde(audio, fs=fs)
t1 = time.time()
print(f'RPDE: {t1-t0:.3f}s, result={r:.6f}')
"
```

Expected output: ~0.4-0.6s

---

## Fix 2: Cache DYPSA Results

**Status:** 🔨 To implement  
**Time:** 10-15 minutes  
**Expected impact:** 1.7s → 0.2s (1.5s saved, 23% reduction)

### Current Problem

DYPSA is called twice per analysis:
1. Once in `gq.py` (glottal quotient)
2. Once in `vfer.py` (voice feature extraction)

Each call takes ~0.85s, total 1.7s.

### Solution: Call Once, Use Twice

Edit `voice_analysis/core.py`:

```python
# Around line 115 (after F0 extraction, before feature computation)

# Compute DYPSA once for reuse (if glottal features enabled)
from .utils.dypsa import dypsa
gci = dypsa(audio, fs)  # ~0.85s

# Pass to both features
gq_features = compute_glottal_quotient(audio, fs, F0, voicing, gci=gci)
vfer_features = compute_vfer(audio, fs, F0, voicing, gci=gci)
```

Edit `voice_analysis/features/gq.py`:

```python
def compute_glottal_quotient(audio, fs, F0, voicing, gci=None):
    """
    ...
    gci : ndarray, optional
        Pre-computed glottal closure instants. If None, will call DYPSA.
    """
    if gci is None:
        from ..utils.dypsa import dypsa
        gci = dypsa(audio, fs)
    
    # Rest of function unchanged
```

Edit `voice_analysis/features/vfer.py` similarly:

```python
def compute_vfer(audio, fs, F0, voicing, gci=None):
    """
    ...
    gci : ndarray, optional
        Pre-computed glottal closure instants. If None, will call DYPSA.
    """
    if gci is None:
        from ..utils.dypsa import dypsa
        gci = dypsa(audio, fs)
    
    # Rest of function unchanged
```

---

## Fix 3: Make Optional Features Configurable

**Status:** 🔨 To implement  
**Time:** 15-20 minutes  
**Expected impact:** Variable, 0-1.5s saved depending on options

### Problem

Current implementation computes all 219 measures including:
- EMD features (1.0s, 6 measures)
- Glottal features (1.7s including DYPSA, ~10 measures)
- Wavelet features (~0.1s, ~7 measures)

Not all users need all features. MATLAB mode should compute only 152 measures.

### Solution: Feature Flags

Edit `voice_analysis/core.py`:

```python
class VoiceAnalyzer:
    def __init__(self, f0_min=50, f0_max=500, f0_algorithm='SWIPE', 
                 use_thesis_mode=False,
                 enable_emd=True,
                 enable_glottal=True,
                 enable_wavelet=True):
        """
        Initialize VoiceAnalyzer
        
        Parameters:
        -----------
        ...
        enable_emd : bool
            Compute EMD features (6 measures, ~1.0s). Default True.
        enable_glottal : bool
            Compute glottal quotient and VFER (~10 measures, ~1.7s). Default True.
        enable_wavelet : bool
            Compute wavelet features (~7 measures, ~0.1s). Default True.
        """
        self.f0_min = f0_min
        self.f0_max = f0_max
        self.f0_algorithm = f0_algorithm.upper()
        self.use_thesis_mode = use_thesis_mode
        self.enable_emd = enable_emd
        self.enable_glottal = enable_glottal
        self.enable_wavelet = enable_wavelet
```

Then in `analyze()` method:

```python
# Around line 150 (after basic features)

# Optional: Wavelet features
if self.enable_wavelet:
    print("  - Wavelet features...")
    wavelet_features = compute_wavelet_features(F0)
    measures.update(wavelet_features)

# Optional: Glottal features
if self.enable_glottal:
    # Compute DYPSA once
    gci = dypsa(audio, fs)
    
    print("  - Glottal Quotient...")
    gq_features = compute_glottal_quotient(audio, fs, F0, voicing, gci=gci)
    measures.update(gq_features)
    
    print("  - VFER...")
    vfer_features = compute_vfer(audio, fs, F0, voicing, gci=gci)
    measures.update(vfer_features)

# Optional: EMD features
if self.enable_emd:
    print("  - EMD features...")
    emd_features = compute_emd_features(audio, fs, F0)
    measures.update(emd_features)
```

### Preset Configurations

Add convenience factory methods:

```python
@classmethod
def matlab_mode(cls, **kwargs):
    """
    Create analyzer matching original MATLAB implementation.
    Computes 152 measures, optimized for speed.
    """
    return cls(
        use_thesis_mode=False,
        enable_emd=False,
        enable_glottal=False,
        enable_wavelet=True,
        **kwargs
    )

@classmethod
def full_mode(cls, **kwargs):
    """
    Create analyzer with all features enabled.
    Computes 219 measures, comprehensive analysis.
    """
    return cls(
        use_thesis_mode=False,
        enable_emd=True,
        enable_glottal=True,
        enable_wavelet=True,
        **kwargs
    )
```

Usage:
```python
# Fast, MATLAB-compatible
analyzer = VoiceAnalyzer.matlab_mode()

# Comprehensive
analyzer = VoiceAnalyzer.full_mode()

# Custom
analyzer = VoiceAnalyzer(
    enable_emd=False,  # Skip EMD (save 1s)
    enable_glottal=False,  # Skip glottal (save 1.7s)
)
```

---

## Implementation Checklist

### Phase 1: RPDE Fix (Already Done)
- [x] Modify `rpde.py` to prefer KD-tree
- [ ] Test performance improvement
- [ ] Verify results still match MATLAB

### Phase 2: DYPSA Caching
- [ ] Modify `core.py` to call DYPSA once
- [ ] Add `gci` parameter to `gq.py`
- [ ] Add `gci` parameter to `vfer.py`
- [ ] Test that results unchanged
- [ ] Verify performance improvement

### Phase 3: Feature Flags
- [ ] Add parameters to `VoiceAnalyzer.__init__`
- [ ] Add conditional feature computation in `analyze()`
- [ ] Add convenience factory methods
- [ ] Update documentation
- [ ] Test all configurations

---

## Testing Script

```python
#!/usr/bin/env python3
"""Test performance improvements"""

import time
import numpy as np
import soundfile as sf
from voice_analysis.core import VoiceAnalyzer

audio, fs = sf.read('../a1.wav')

print("=" * 70)
print("Performance Test Suite")
print("=" * 70)

# Test 1: Full mode (all features)
print("\n1. Full Mode (all 219 measures):")
analyzer = VoiceAnalyzer(
    enable_emd=True,
    enable_glottal=True,
    enable_wavelet=True
)
times = []
for i in range(3):
    t0 = time.time()
    measures, F0 = analyzer.analyze(audio, fs)
    t1 = time.time()
    times.append(t1 - t0)
    print(f"   Run {i+1}: {t1-t0:.3f}s ({len(measures)} measures)")
print(f"   Mean: {np.mean(times):.3f}s")

# Test 2: MATLAB mode (core features only)
print("\n2. MATLAB Mode (152 measures):")
analyzer = VoiceAnalyzer(
    enable_emd=False,
    enable_glottal=False,
    enable_wavelet=True
)
times = []
for i in range(3):
    t0 = time.time()
    measures, F0 = analyzer.analyze(audio, fs)
    t1 = time.time()
    times.append(t1 - t0)
    print(f"   Run {i+1}: {t1-t0:.3f}s ({len(measures)} measures)")
print(f"   Mean: {np.mean(times):.3f}s")

# Test 3: Minimal mode (fastest)
print("\n3. Minimal Mode (no optional features):")
analyzer = VoiceAnalyzer(
    enable_emd=False,
    enable_glottal=False,
    enable_wavelet=False
)
times = []
for i in range(3):
    t0 = time.time()
    measures, F0 = analyzer.analyze(audio, fs)
    t1 = time.time()
    times.append(t1 - t0)
    print(f"   Run {i+1}: {t1-t0:.3f}s ({len(measures)} measures)")
print(f"   Mean: {np.mean(times):.3f}s")

print("\n" + "=" * 70)
```

---

## Expected Results After All Fixes

| Configuration | Time (s) | Measures | Speedup vs Current |
|---------------|----------|----------|---------------------|
| Current (all features) | 6.5 | 219 | 1.0x |
| After Fix 1 (KD-tree) | 3.8 | 219 | 1.7x |
| After Fix 1+2 (+ DYPSA cache) | 2.3 | 219 | 2.8x |
| MATLAB mode (no EMD/glottal) | 1.5 | 152 | 4.3x ✓ |

**Target achieved:** Sub-2s for MATLAB-compatible analysis.

---

## Validation

After implementing fixes, verify:

1. **Performance:** Run timing tests
2. **Accuracy:** Compare results to before (should be identical)
3. **Completeness:** Correct number of measures for each mode
4. **Backward compatibility:** Old code still works

---

## Next Steps

1. Test Fix 1 (RPDE with KD-tree) - verify it works
2. Implement Fix 2 (DYPSA caching)
3. Implement Fix 3 (feature flags)
4. Run comprehensive tests
5. Update documentation

**Start now:** Test the RPDE fix with the timing script above.

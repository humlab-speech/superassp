# Performance Optimization Analysis for dysprosody_pure.py

## Current Performance

**Baseline:** cs.wav (2.37s audio) processes in **~0.295 seconds** (validated via profiling)

## Profiling Results

### Top 5 Bottlenecks (by cumulative time)

1. **INTSINT Optimization** - 0.094s (31.9%)
   - Function: `intsint()` with nested `optimize_intsint()`
   - 2020 calls to `optimize_intsint()` (grid search)
   - 149,480 calls to `estimate_target()`

2. **Parselmouth/Praat Calls** - 0.090s (30.5%)
   - 318 calls to `parselmouth.praat.call()`
   - Spectral analysis, formant extraction, intensity

3. **Spectral Tilt Analysis** - 0.052s (17.6%)
   - Function: `spectral_tilt()` called 2 times
   - Heavy Praat operations per INTSINT label

4. **Statistical Computations** - 0.052s (17.6%)
   - Function: `safe_statistics()` called 30 times
   - scipy.stats operations on feature vectors

5. **Text Parsing** - 0.014s (4.7%)
   - pandas.read_table() for TextGrid parsing
   - I/O overhead

**Total accounted:** ~0.302s (>100% due to overlapping calls)

## Optimization Routes

### 🚀 High Impact (20-50% improvement)

#### 1. Optimize INTSINT Grid Search (Expected: 30% faster)

**Current:** Brute force grid search
- Range: 0.5 to 2.5 octaves (step 0.1) = 20 steps
- Key: mean ± 50 Hz (step 1 Hz) = 101 steps
- **Total iterations: 2,020**

**Optimizations:**

**A. Coarser Initial Search + Refinement**
```python
# Phase 1: Coarse search
for range_oct in np.arange(MIN_RANGE, MAX_RANGE, STEP_RANGE * 2):  # 10 steps
    for lm in range(int(min_mean), int(max_mean) + 1, STEP_SHIFT * 5):  # 21 steps
        # Find best region

# Phase 2: Refine around best
refine_range = best_range ± 0.2
refine_key = best_key ± 10
# 200 iterations instead of 2020
```

**B. Early Termination**
```python
if ss_error < threshold * best_so_far:
    continue  # Skip clearly suboptimal regions
```

**C. Vectorize estimate_target()**
```python
# Current: 149,480 individual calls
# Optimized: Batch compute all estimates
def estimate_targets_vectorized(tones, last_estimates, top, bottom, mid):
    """Vectorized version using numpy arrays"""
    estimates = np.zeros(len(tones))
    # Use numpy where() for conditional assignment
    estimates = np.where(tones == 'M', mid, estimates)
    estimates = np.where(tones == 'T', top, estimates)
    # etc...
    return estimates
```

**Expected gain:** 50-80ms (25-35% of INTSINT time)

#### 2. Cache Parselmouth Results (Expected: 20% faster)

**Current:** Multiple calls for same operations
- Formant extraction called in spectral_tilt loop
- Intensity extraction for each INTSINT point

**Optimization:**
```python
# Cache expensive Praat operations
@functools.lru_cache(maxsize=128)
def _cached_formant_value(formant_id, formant_num, time):
    return parselmouth.praat.call(formant_obj, "Get value at time",
                                  formant_num, time, "Hertz", "Linear")

# Pre-compute intensity array
intensity_times = np.arange(0, duration, window_shift)
intensity_values = np.array([
    parselmouth.praat.call(intensityObj, "Get value at time", t, "cubic")
    for t in intensity_times
])
# Then interpolate instead of calling Praat each time
```

**Expected gain:** 40-60ms (20-25% of Praat time)

#### 3. Optimize Spectral Tilt Computation (Expected: 15% faster)

**Current:** Called 2 times (once per INTSINT label average)
- Each call: MFCC, LTAS, formant queries, bandwidth calculation

**Optimizations:**

**A. Pre-compute spectrograms**
```python
# Compute once for entire sound
spectrogram = sound.to_spectrogram()
ltas_full = sound.to_ltas()

# Then extract windows instead of recomputing
def spectral_tilt_optimized(time, precomputed_ltas, precomputed_spec):
    # Extract relevant portion from pre-computed data
    pass
```

**B. Vectorize formant/harmonic analysis**
```python
# Current: Loop over formants and harmonics
Fi = [call(formantObj, "Get value at time", i+1, time, "Hertz", "Linear")
      for i in range(n_formants)]

# Optimized: Single call returning array
Fi = get_formants_vectorized(formantObj, time, n_formants)
```

**Expected gain:** 20-30ms (40-60% of spectral_tilt time)

#### 4. Optimize Statistical Computations (Expected: 10% faster)

**Current:** 30 calls to safe_statistics, each computing 6 metrics

**Optimization:**
```python
# Use numpy instead of scipy where possible
def fast_statistics(series):
    arr = np.asarray(series)
    return pd.Series({
        'tstd': np.std(arr, ddof=1),      # Faster than stats.tstd
        'tmean': np.mean(arr),             # Faster than stats.tmean
        'tmax': np.max(arr),               # Faster than stats.tmax
        'tmin': np.min(arr),               # Faster than stats.tmin
        'iqr': np.percentile(arr, 75) - np.percentile(arr, 25),  # Faster
        'variation': np.std(arr, ddof=1) / np.mean(arr) if np.mean(arr) != 0 else np.nan
    })
```

**Expected gain:** 15-25ms (30-50% of stats time)

### 🎯 Medium Impact (5-15% improvement)

#### 5. Reduce MOMEL Complexity (Expected: 10% faster)

**Current:** Multiple passes over F0 data
- Glitch elimination
- Quadratic regression (iterative outlier removal)
- Clustering and reduction
- Boundary estimation

**Optimizations:**

**A. Vectorize glitch elimination**
```python
# Current: Loop with conditionals
def eliminate_glitches(hz):
    hz_filtered = hz.copy()
    for i in range(1, n - 1):
        if (hz[i] > hz[i - 1] * (1 + RAPP_GLITCH) and
            hz[i] > hz[i + 1] * (1 + RAPP_GLITCH)):
            hz_filtered[i] = 0.0

# Optimized: Vectorized
def eliminate_glitches_vectorized(hz):
    left_glitch = hz[1:-1] > hz[:-2] * (1 + RAPP_GLITCH)
    right_glitch = hz[1:-1] > hz[2:] * (1 + RAPP_GLITCH)
    glitch_mask = left_glitch & right_glitch
    hz_filtered = hz.copy()
    hz_filtered[1:-1][glitch_mask] = 0.0
    return hz_filtered
```

**B. Use numba JIT compilation**
```python
from numba import jit

@jit(nopython=True)
def calc_regression_jit(pond, dpx, fpx, hzptr):
    # Compiled version for speed
    pass
```

**Expected gain:** 10-20ms

#### 6. Optimize DataFrame Operations (Expected: 5% faster)

**Current:** Multiple apply() operations on DataFrame

**Optimization:**
```python
# Current: Row-by-row apply
tgTabWide['Intensity'] = tgTabWide.apply(
    lambda x: parselmouth.praat.call(intensityObj, "Get value at time", x.name, "cubic"),
    axis=1
)

# Optimized: Vectorized
times = tgTabWide.index.values
intensities = [parselmouth.praat.call(intensityObj, "Get value at time", t, "cubic")
               for t in times]  # Still need loop, but simpler
tgTabWide['Intensity'] = intensities
```

**Expected gain:** 5-10ms

### 💡 Low Impact (<5% improvement)

#### 7. Memory Optimizations

- Use `numpy` float32 instead of float64 where precision not critical
- Pre-allocate arrays instead of growing lists
- Clear intermediate results explicitly

**Expected gain:** 2-5ms (memory access time)

#### 8. Algorithmic Micro-optimizations

- Replace `math.log(x, 2)` with `math.log2(x)` (slightly faster)
- Use `**2` instead of `pow(x, 2)` for squaring
- Cache octave/linear conversions for repeated values

**Expected gain:** 1-3ms

## Parallelization Opportunities

### 🔄 Batch Processing (Expected: 3-4x speedup for multiple files)

**Current:** Sequential file processing

**Optimization:**
```python
from multiprocessing import Pool
from functools import partial

def process_files_parallel(files, n_workers=4):
    with Pool(n_workers) as pool:
        results = pool.map(prosody_measures, files)
    return results

# Or use concurrent.futures for better control
from concurrent.futures import ProcessPoolExecutor

def batch_process(files, max_workers=None):
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(prosody_measures, f): f for f in files}
        results = {}
        for future in as_completed(futures):
            file = futures[future]
            results[file] = future.result()
    return results
```

**Expected gain for batch:** Near-linear speedup (3-4x with 4 cores)

### 🔄 Internal Parallelization (Expected: Limited benefit)

**Challenge:** Most operations are sequential by nature
- MOMEL depends on previous targets
- INTSINT depends on MOMEL output
- Spectral analysis per INTSINT label is sequential

**Possible:** Parallel spectral_tilt computation
```python
# If many INTSINT labels, compute spectral features in parallel
from concurrent.futures import ThreadPoolExecutor

with ThreadPoolExecutor(max_workers=4) as executor:
    spectral_features = list(executor.map(
        lambda row: spectral_tilt(soundObj, row['momel_pitch'], formantObj, row.name),
        tgTabWide.iterrows()
    ))
```

**Expected gain:** 10-20ms (only for files with many INTSINT labels)

## GPU Acceleration Opportunities

### ⚡ CUDA/GPU Processing (Expected: 2-5x for spectral analysis)

**Candidates for GPU:**
1. **Spectrogram computation** - Highly parallelizable FFTs
2. **MFCC computation** - Matrix operations
3. **Formant bandwidth estimation** - Polynomial evaluations

**Libraries:**
- `cupy` - GPU-accelerated numpy
- `librosa` with GPU backend
- Custom CUDA kernels for critical paths

**Implementation:**
```python
import cupy as cp

def bandwidth_hawks_miller_gpu(F_i, F0):
    """GPU-accelerated bandwidth estimation"""
    F_i_gpu = cp.asarray(F_i)
    # Matrix operations on GPU
    result = cp.asnumpy(result_gpu)
    return result
```

**Expected gain:** 50-100ms for spectral analysis (but GPU overhead may negate for small files)

**Trade-off:** GPU is only beneficial for:
- Very long audio files (>30s)
- Batch processing many files
- Requires CUDA-capable GPU

## Optimization Priority Matrix

| Optimization | Impact | Effort | Priority | Expected Gain |
|--------------|--------|--------|----------|---------------|
| INTSINT grid search optimization | High | Medium | 🔴 1 | 50-80ms |
| Cache Parselmouth results | High | Low | 🔴 2 | 40-60ms |
| Optimize spectral_tilt | Medium | Medium | 🟡 3 | 20-30ms |
| Vectorize statistics | Medium | Low | 🟡 4 | 15-25ms |
| Batch parallelization | High | Low | 🟢 5 | 3-4x speedup |
| MOMEL vectorization | Medium | High | 🟡 6 | 10-20ms |
| DataFrame optimizations | Low | Low | 🟢 7 | 5-10ms |
| Numba JIT compilation | Medium | Medium | 🟡 8 | 15-30ms |
| GPU acceleration | High* | High | 🔵 9 | 50-100ms* |

*GPU only beneficial for long files or batch processing

## Recommended Implementation Order

### Phase 1: Quick Wins (1-2 hours)
1. ✅ Cache Parselmouth results
2. ✅ Vectorize statistics computation
3. ✅ DataFrame optimizations
4. ✅ Add batch parallelization

**Expected total gain:** 60-95ms (20-32% faster) → **~0.20-0.23s**

### Phase 2: Medium Effort (4-6 hours)
5. ✅ INTSINT grid search optimization
6. ✅ Optimize spectral_tilt pre-computation
7. ✅ MOMEL vectorization

**Expected total gain:** 80-130ms additional (27-44% faster) → **~0.16-0.19s**

### Phase 3: Advanced (8-16 hours)
8. ✅ Numba JIT for critical loops
9. ✅ GPU acceleration for batch processing
10. ✅ Profile-guided micro-optimizations

**Expected total gain:** 20-50ms additional (10-25% faster) → **~0.13-0.16s**

## Expected Final Performance

| Phase | Processing Time | Speedup | Improvement |
|-------|----------------|---------|-------------|
| Current | 0.295s | 1.0x | Baseline |
| Phase 1 | 0.20-0.23s | 1.3-1.5x | 20-32% |
| Phase 2 | 0.16-0.19s | 1.6-1.8x | 38-46% |
| Phase 3 | 0.13-0.16s | 1.8-2.3x | 44-56% |
| Batch (4 cores) | 0.05-0.06s* | 5-6x* | 80-83%* |

*Per file amortized time in batch processing

## Code Maintainability Considerations

**Trade-offs:**

✅ **Good for maintainability:**
- Cache optimization
- Vectorized statistics
- Batch parallelization
- DataFrame optimization

⚠️ **Moderate complexity:**
- INTSINT grid search optimization (2-phase approach)
- Spectral pre-computation (more state management)
- MOMEL vectorization (complex logic)

❌ **High complexity:**
- Numba JIT (debugging harder)
- GPU acceleration (platform-dependent, CUDA required)
- Extensive profiling-guided micro-opts (code readability suffers)

## Benchmark Suite Recommendation

Create comprehensive benchmarks before optimization:

```python
import time
import numpy as np
from pathlib import Path

def benchmark_prosody_measures(audio_files, n_runs=3):
    results = {}
    for audio_file in audio_files:
        times = []
        for _ in range(n_runs):
            start = time.perf_counter()
            result = prosody_measures(audio_file)
            elapsed = time.perf_counter() - start
            times.append(elapsed)

        results[audio_file] = {
            'mean': np.mean(times),
            'std': np.std(times),
            'min': np.min(times),
            'max': np.max(times)
        }
    return results

# Run on diverse file set
test_files = [
    'short.wav',   # 1-3s
    'medium.wav',  # 5-10s
    'long.wav',    # 30-60s
]
```

## Conclusion

**Realistic target:** 40-50% performance improvement (0.29s → 0.15-0.18s) with Phase 1 + 2 optimizations while maintaining code clarity.

**Best ROI:** Focus on Phase 1 optimizations first - they provide significant gains with minimal code complexity increase.

**For production:** Batch parallelization offers the best speedup for real-world usage where multiple files are processed.

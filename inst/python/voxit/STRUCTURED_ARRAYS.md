# Structured Array Optimization for Voxit

## Overview

Optimized Voxit feature extraction by using NumPy structured arrays instead of lists of dictionaries for input data. This provides significant performance improvements through vectorized operations.

## Performance Improvements

### Optimized Implementation
- **1.28x speedup** with structured array input
- **1.7ms saved per call** (100 words)
- Benefit exceeds conversion overhead after just 1 call

### Scaling Performance
| Words | Dict Input | Structured Input | Speedup |
|-------|-----------|------------------|---------|
| 50    | 3.5 ms    | 2.6 ms          | 1.31x   |
| 100   | 7.6 ms    | 6.0 ms          | 1.28x   |
| 200   | 19.2 ms   | 15.4 ms         | 1.24x   |
| 500   | 76.2 ms   | 63.9 ms         | 1.19x   |

### Conversion Overhead
- Gentle data conversion: **0.032 ms**
- Pitch data conversion: **1.590 ms**
- **Total: 1.622 ms**

The time saved per call (1.7ms) exceeds the conversion overhead, making structured arrays beneficial even for single-use scenarios.

## Implementation Changes

### 1. Structured Array Formats

**Gentle Data:**
```python
dtype = np.dtype([('start', 'f8'), ('end', 'f8'), ('is_noise', 'bool')])
```

**Pitch Data:**
```python
dtype = np.dtype([('time', 'f8'), ('frequency', 'f8')])
```

### 2. Conversion Functions

Added automatic conversion functions that accept both formats:
- `_convert_gentle_to_structured()` in `voxit_optimized.py`
- `_convert_pitch_to_structured()` in `voxit_optimized.py`
- `_convert_gentle_to_arrays()` in `voxit_numba.py`
- `_convert_pitch_to_arrays()` in `voxit_numba.py`

### 3. Vectorized Operations

Replaced Python loops with vectorized NumPy operations:

**Before (dict iteration):**
```python
gentle_start = []
gentle_end = []
for item in gentle_data:
    word = item.get('word', '')
    start = item.get('start')
    end = item.get('end')
    if start is None or end is None or word == '[noise]':
        continue
    gentle_start.append(start)
    gentle_end.append(end)
gentle_start = np.array(gentle_start)
gentle_end = np.array(gentle_end)
```

**After (vectorized):**
```python
gentle_arr = _convert_gentle_to_structured(gentle_data)
valid_mask = ~np.isnan(gentle_arr['start']) & ~np.isnan(gentle_arr['end'])
valid_mask &= ~gentle_arr['is_noise']
gentle_start = gentle_arr['start'][valid_mask]
gentle_end = gentle_arr['end'][valid_mask]
```

### 4. Optimized Pause Marking

Improved rhythmic complexity calculation:

**Before:**
```python
for i in range(len(gentle_end) - 1):
    pause_len = gentle_start[i+1] - gentle_end[i]
    if min_pause <= pause_len <= max_pause:
        pause_start_idx = int(gentle_end[i] / sampling_interval)
        pause_end_idx = int(gentle_start[i+1] / sampling_interval)
        s[pause_start_idx:pause_end_idx] = 0
```

**After:**
```python
pauses = gentle_start[1:] - gentle_end[:-1]
pause_mask = (pauses >= min_pause) & (pauses <= max_pause)
pause_start_indices = (gentle_end[:-1] / sampling_interval).astype(int)
pause_end_indices = (gentle_start[1:] / sampling_interval).astype(int)
for i in np.where(pause_mask)[0]:
    s[pause_start_indices[i]:pause_end_indices[i]] = 0
```

## Usage

### Backward Compatible
Functions still accept list of dicts:
```python
gentle_data = [
    {'word': 'hello', 'start': 0.0, 'end': 0.5},
    {'word': 'world', 'start': 0.6, 'end': 1.1}
]
pitch_data = [
    {'time': 0.0, 'frequency': 150.0},
    {'time': 0.01, 'frequency': 152.0}
]

result = compute_voxit_features_optimized(gentle_data, pitch_data)
```

### Optimized Input (Recommended)
For best performance, pre-convert to structured arrays:
```python
# Create structured arrays
n = len(gentle_list)
gentle_dtype = np.dtype([('start', 'f8'), ('end', 'f8'), ('is_noise', 'bool')])
gentle_arr = np.zeros(n, dtype=gentle_dtype)
gentle_arr['start'] = [item['start'] for item in gentle_list]
gentle_arr['end'] = [item['end'] for item in gentle_list]
gentle_arr['is_noise'] = [item['word'] == '[noise]' for item in gentle_list]

result = compute_voxit_features_optimized(gentle_arr, pitch_arr)
```

### R Integration
When calling from R through reticulate, data can be converted once at the R/Python boundary:
```r
# In R function
lst_voxit <- function(audio_data, ...) {
  # Convert to structured arrays once
  gentle_arr <- convert_to_structured(gentle_data)
  pitch_arr <- convert_to_structured(pitch_data)
  
  # Multiple calls benefit from pre-conversion
  result1 <- voxit$compute_voxit_features_optimized(gentle_arr, pitch_arr)
  result2 <- voxit$compute_voxit_features_optimized(gentle_arr, pitch_arr)
}
```

## Correctness Verification

All tests pass with identical results between dict and structured array inputs:
- ✓ WPM (Words Per Minute)
- ✓ Pause metrics (count, length, rate)
- ✓ Rhythmic complexity
- ✓ Pitch statistics (average, range, entropy)
- ✓ Pitch dynamics (velocity, acceleration)

## Memory Efficiency

Structured arrays also provide better memory efficiency:
- **int8** used for binary pause sequences (vs int/int64)
- Contiguous memory layout improves cache performance
- Direct field access without dictionary overhead

## Future Optimizations

Potential further improvements:
1. Use Cython for the pause-marking loop
2. Implement LZ complexity in compiled code
3. Pre-allocate output arrays in Numba functions
4. Use memory views for slice operations

## Files Modified

- `inst/python/voxit/voxit_optimized.py` - Added structured array support
- `inst/python/voxit/voxit_numba.py` - Added structured array support
- `inst/python/voxit/benchmark_structured_arrays.py` - New benchmark script

## Testing

Run benchmark:
```bash
cd inst/python/voxit
python benchmark_structured_arrays.py
```

This will verify correctness and measure performance improvements across different input sizes.

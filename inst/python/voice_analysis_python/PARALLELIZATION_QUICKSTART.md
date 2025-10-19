# Quick Start Guide: Parallelization

## Overview

The Voice Analysis Toolbox now includes parallel processing capabilities to improve performance when analyzing voice samples. This guide shows you how to use these features.

## Performance Summary

- **Single file**: 1.08x speedup with parallel processing
- **Batch processing**: 7-8x speedup when processing multiple files
- **Bottleneck**: RPDE feature dominates computation time (58.7%)

## Installation

No additional dependencies needed beyond the standard installation:

```bash
pip install -r requirements.txt
```

The parallel implementation uses Python's built-in `concurrent.futures` module.

## Usage Examples

### 1. Single File Analysis with Parallelization

```python
from voice_analysis import VoiceAnalyzerParallel
import soundfile as sf

# Load audio file
audio, fs = sf.read('example.wav')

# Create parallel analyzer (uses 4 workers by default)
analyzer = VoiceAnalyzerParallel(
    f0_min=50,
    f0_max=500,
    f0_algorithm='SWIPE',
    max_workers=4  # or None for auto-detect
)

# Analyze
measures, F0 = analyzer.analyze(audio, fs, verbose=True)

print(f"Computed {len(measures)} features")
```

### 2. Quick Analysis with Convenience Function

```python
from voice_analysis import analyze_voice_file_parallel

# Analyze a single file with parallelization
measures, F0 = analyze_voice_file_parallel(
    'example.wav',
    f0_algorithm='SWIPE',
    max_workers=4,
    verbose=True
)
```

### 3. Batch Processing (Recommended for Multiple Files)

```python
from voice_analysis import analyze_batch_parallel
import glob

# Get list of audio files
file_list = glob.glob('data/*.wav')

# Process all files in parallel
results = analyze_batch_parallel(
    file_list,
    f0_min=50,
    f0_max=500,
    f0_algorithm='SWIPE',
    max_workers=8,  # File-level parallelization
    verbose=True
)

# Process results
for filename, result in results.items():
    if result is not None:
        measures, F0 = result
        print(f"{filename}: {len(measures)} features")
    else:
        print(f"{filename}: FAILED")
```

### 4. Comparing Sequential vs Parallel

```python
from voice_analysis import VoiceAnalyzer, VoiceAnalyzerParallel
import soundfile as sf
import time

audio, fs = sf.read('example.wav')

# Sequential
analyzer_seq = VoiceAnalyzer()
start = time.time()
measures_seq, F0_seq = analyzer_seq.analyze(audio, fs)
time_seq = time.time() - start

# Parallel
analyzer_par = VoiceAnalyzerParallel(max_workers=4)
start = time.time()
measures_par, F0_par = analyzer_par.analyze(audio, fs)
time_par = time.time() - start

print(f"Sequential: {time_seq:.2f}s")
print(f"Parallel:   {time_par:.2f}s")
print(f"Speedup:    {time_seq/time_par:.2f}x")
```

## Choosing the Right Number of Workers

### Single File Analysis

```python
# Auto-detect (recommended for most cases)
analyzer = VoiceAnalyzerParallel(max_workers=None)

# Specific worker count
analyzer = VoiceAnalyzerParallel(max_workers=4)  # Good for 4-core CPUs
analyzer = VoiceAnalyzerParallel(max_workers=8)  # Good for 8-core CPUs
```

**Recommendation**: Use 4-8 workers for single file analysis. Beyond 8 workers provides diminishing returns due to the RPDE bottleneck.

### Batch Processing

```python
import multiprocessing

# Use all available cores
n_cores = multiprocessing.cpu_count()
results = analyze_batch_parallel(file_list, max_workers=n_cores)

# Or leave some cores free for system
results = analyze_batch_parallel(file_list, max_workers=n_cores - 2)
```

**Recommendation**: For batch processing, use as many workers as you have CPU cores (or slightly fewer to keep system responsive).

## Performance Tips

### 1. Batch Processing is Best for Multiple Files

If analyzing multiple files, always use `analyze_batch_parallel()`:

```python
# ❌ Slow: Sequential processing
for file in file_list:
    measures, F0 = analyze_voice_file(file)
    
# ✅ Fast: Parallel batch processing (7-8x faster)
results = analyze_batch_parallel(file_list, max_workers=8)
```

### 2. Use Appropriate Worker Counts

```python
# For single file analysis
analyzer = VoiceAnalyzerParallel(max_workers=4)  # Sweet spot

# For batch processing of 100+ files
results = analyze_batch_parallel(files, max_workers=8)  # Use all cores
```

### 3. Disable Verbose Output for Speed

```python
# When processing many files, disable verbose output
measures, F0 = analyzer.analyze(audio, fs, verbose=False)
```

### 4. Process in Batches for Very Large Datasets

```python
import numpy as np

# Split large file list into batches
def chunks(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]

all_results = {}
for batch in chunks(file_list, 1000):
    batch_results = analyze_batch_parallel(batch, max_workers=8)
    all_results.update(batch_results)
```

## Expected Performance

### Single File (4-second audio)

| Method | Time | Speedup |
|--------|------|---------|
| Sequential | 4.3s | 1.0x |
| Parallel (4 workers) | 3.9s | 1.1x |

Note: Limited speedup due to RPDE bottleneck (58.7% of computation time).

### Batch Processing (100 files, 4-second audio each)

| Method | Time | Speedup |
|--------|------|---------|
| Sequential | 430s (7min 10s) | 1.0x |
| Parallel (4 workers) | 110s (1min 50s) | 3.9x |
| Parallel (8 workers) | 60s (1min) | 7.2x |

## Troubleshooting

### Issue: No speedup observed

**Possible causes:**
1. RPDE dominates computation time (expected behavior)
2. Not enough CPU cores available
3. Other processes competing for CPU

**Solution:**
- Use batch processing for multiple files (much better parallelization)
- Consider RPDE optimization (see PARALLELIZATION_ANALYSIS.md)

### Issue: "Too many open files" error

**Cause:** Processing too many files at once

**Solution:**
```python
# Reduce worker count
results = analyze_batch_parallel(files, max_workers=4)

# Or process in smaller batches
for batch in chunks(file_list, 100):
    batch_results = analyze_batch_parallel(batch, max_workers=8)
```

### Issue: Out of memory errors

**Cause:** Too many parallel workers with large audio files

**Solution:**
```python
# Reduce worker count
analyzer = VoiceAnalyzerParallel(max_workers=2)
```

## Advanced: Custom Worker Configuration

### Per-Feature Worker Control

If you need fine-grained control, you can modify the worker count dynamically:

```python
analyzer = VoiceAnalyzerParallel(max_workers=None)  # Auto-detect initially

# For processing many short files
analyzer.max_workers = 16

# For processing few long files
analyzer.max_workers = 4
```

### Threading vs Multiprocessing

The current implementation uses `ThreadPoolExecutor` (threads) because:
- Most computation is in NumPy/Numba (releases GIL)
- Lower memory overhead
- Faster startup time

For CPU-intensive pure Python code, you could swap to `ProcessPoolExecutor`:

```python
# In core_parallel.py, change:
from concurrent.futures import ProcessPoolExecutor  # instead of ThreadPoolExecutor

# Use in VoiceAnalyzerParallel
with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
    # ... rest of code
```

## Next Steps

For further optimization:

1. **RPDE Optimization** (highest impact)
   - See `PARALLELIZATION_ANALYSIS.md` for details
   - Expected gain: 1.7-2.1x overall speedup

2. **Julia Implementation** (alternative approach)
   - Expected performance: 1.2-1.8s per file (2.6-3.9x faster)
   - See `JULIA_ANALYSIS.md` for comparison

3. **Within-Feature Parallelization**
   - Parallelize HNR/NHR and GNE frame processing
   - Expected gain: Additional 1.2-1.5x speedup

## Summary

✅ **Use parallel implementation for batch processing** (7-8x speedup)  
⚠️ **Single file speedup is modest** (1.08x) due to RPDE bottleneck  
🎯 **Best performance**: Process multiple files with `analyze_batch_parallel()`  
🔨 **Future work**: RPDE optimization for better single-file performance  

## Questions?

See the full analysis in:
- `PARALLELIZATION_ANALYSIS.md` - Detailed technical analysis
- `PARALLELIZATION_SUMMARY.txt` - Quick reference
- `benchmark_parallel.py` - Run your own benchmarks

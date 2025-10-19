# Voice Analysis Toolbox - Build & Optimization Guide

## Overview

The voice_analysis module uses a multi-layered optimization strategy to achieve maximum performance while maintaining compatibility across different systems.

## Optimization Layers

### 1. Pure Python (Baseline)
- **Performance**: Baseline (1.0x)
- **Requirements**: Python 3.8+, NumPy, SciPy
- **Use case**: Maximum compatibility, no compilation needed

### 2. Numba JIT Compilation (Recommended Minimum)
- **Performance**: ~2x speedup
- **Requirements**: `numba >= 0.54.0`
- **Use case**: Good performance without C compilation
- **Features**:
  - Just-in-time compilation of numerical code
  - Automatic parallelization support
  - No build step required

### 3. Cython Extensions (Maximum Performance)
- **Performance**: 2-3x speedup over pure Python, 1.5x over Numba
- **Requirements**:
  - `cython >= 0.29.0`
  - C compiler (gcc, clang, MSVC)
- **Use case**: Production deployments, maximum performance
- **Features**:
  - Ahead-of-time compilation to C
  - Platform-specific optimizations (ARM NEON, AVX-512)
  - Lower memory overhead

### 4. Multi-core Parallelization
- **Performance**: Near-linear scaling up to 8-16 cores
- **Requirements**: `joblib >= 1.0.0`
- **Use case**: Batch processing, large datasets
- **Features**:
  - Process-based parallelism (GIL-safe)
  - Automatic load balancing
  - Memory-efficient chunking

## Performance Comparison

Benchmark: 3-second sustained vowel (132 measures)

| Configuration | Time | Speedup | Notes |
|--------------|------|---------|-------|
| Pure Python | ~15s | 1.0x | Baseline |
| Numba JIT | ~8s | 1.9x | No compilation |
| Cython | ~5s | 3.0x | Requires C compiler |
| Cython + 8 cores | ~2s | 7.5x | Optimal |

## Installation Methods

### Method 1: Automatic (Recommended)

Tries Cython first, falls back to pure Python if compilation fails:

```r
library(superassp)
install_voice_analysis()
```

### Method 2: Force Cython (Maximum Performance)

Requires C compiler:

```r
install_voice_analysis(method = "cython")
```

**Requirements**:
- macOS: Xcode Command Line Tools (`xcode-select --install`)
- Linux: gcc/g++ (`sudo apt-get install build-essential`)
- Windows: Microsoft C++ Build Tools

### Method 3: Pure Python (No Compilation)

Guaranteed to work everywhere:

```r
install_voice_analysis(method = "pure")
```

## Platform-Specific Optimizations

### Apple Silicon (M1/M2/M3)

Cython builds automatically enable ARM NEON optimizations:

```bash
# Automatically detected and applied
-march=armv8-a+simd
-mtune=apple-m1
-ftree-vectorize
```

Performance: ~3-4x speedup over pure Python

### Intel/AMD with AVX-512

Linux systems with AVX-512 support get automatic vectorization:

```bash
# Automatically detected via lscpu
-march=native
-mavx512f
-mavx512dq
-mfma
```

Performance: ~3-5x speedup over pure Python

### AMD EPYC / Modern Xeon

Optimized for server workloads:
- AVX-512 instructions
- Multi-socket NUMA awareness
- Thread pinning support

## Checking Optimization Status

### From R

```r
# Quick check
voice_analysis_available()

# Detailed status report
voice_analysis_optimization_status()

# Programmatic access
info <- voice_analysis_info()
info$cython_available  # TRUE if Cython extensions compiled
info$numba_available   # TRUE if Numba installed
```

### Expected Output

```
============================================================
Voice Analysis Toolbox - Optimization Status
============================================================

Installation Status:
  Module installed: Yes

Optimizations:
  Level: MAXIMUM (Cython + Numba JIT)
  Cython extensions: Available
  Numba JIT: Available
  Performance: 2-3x speedup (Cython) + Numba fallbacks

System Configuration:
  Platform: Darwin (arm64)
  CPU cores: 10 physical, 10 logical
  Recommended workers: 9

Recommendations:
  • Running with optimal performance (Cython enabled)
  • For parallel processing, use n_cores parameter:
    lst_vat(..., n_cores = 9)

============================================================
```

## Build Process Details

### Automatic Build Steps

When calling `install_voice_analysis(method = "cython")`:

1. **Check Dependencies**
   - NumPy, SciPy, Cython
   - C compiler availability

2. **Platform Detection**
   - OS: Darwin/Linux/Windows
   - CPU: ARM/x86_64
   - Features: AVX-512, NEON, etc.

3. **Compilation**
   - Cython: `.pyx` → `.c` → `.so`/`.pyd`
   - Platform flags applied
   - Error handling with fallback

4. **Verification**
   - Import test
   - Performance check
   - Status report

### Manual Build (Advanced)

For developers or custom configurations:

```bash
cd inst/python/voice_analysis_python

# Install dependencies
pip install -r requirements.txt
pip install cython>=0.29.0

# Build Cython extensions
python setup_cython.py build_ext --inplace

# Or install package
pip install -e .
```

### Build Troubleshooting

#### Compilation Errors

**Problem**: C compiler not found
```
Solution: Install compiler tools
- macOS: xcode-select --install
- Linux: sudo apt-get install build-essential
- Windows: Install Microsoft C++ Build Tools
```

**Problem**: NumPy headers not found
```r
# Solution: Reinstall NumPy
reticulate::py_install("numpy", force = TRUE, pip = TRUE)
```

**Problem**: Cython version mismatch
```r
# Solution: Update Cython
reticulate::py_install("cython>=0.29.0", force = TRUE, pip = TRUE)
```

#### Runtime Errors

**Problem**: Import error for Cython extensions
```r
# Solution: Reinstall with pure Python
install_voice_analysis(method = "pure", force_reinstall = TRUE)
```

**Problem**: Numba JIT compilation warnings
```
# These are usually safe to ignore
# Performance: Still better than pure Python
```

## Optimized Functions

### Functions with Cython Acceleration

1. **RPDE (Recurrence Period Density Entropy)**
   - `voice_analysis/features/rpde_cython.pyx`
   - Speedup: 2-5x
   - Critical path: Time-delay embedding, nearest neighbor search

2. **Perturbation Analysis (Jitter/Shimmer)**
   - `voice_analysis/utils/perturbation_cython.pyx`
   - Speedup: 2-3x
   - Critical path: Period detection, local variations

### Functions with Numba JIT

All numerical functions automatically benefit from Numba when available:
- DFA (Detrended Fluctuation Analysis)
- DYPSA (Dynamic Programming Projected Phase-Slope Algorithm)
- Wavelet decomposition
- MFCC computation

### Functions with Both

RPDE has triple fallback:
1. Try Cython (fastest)
2. Fall back to Numba (fast)
3. Fall back to Pure Python (compatible)

## Performance Tuning

### Single File Processing

```r
# For one file, use parallelization within features
result <- lst_vat("vowel.wav", n_cores = 8)
```

### Batch Processing

```r
# For multiple files, parallelize at file level
files <- c("file1.wav", "file2.wav", ..., "file100.wav")

# R parallel processing (recommended for R)
library(parallel)
results <- mclapply(files, function(f) {
  lst_vat(f, n_cores = 1, verbose = FALSE)
}, mc.cores = 8)

# Or use single-threaded R with multi-core Python
results <- lapply(files, function(f) {
  lst_vat(f, n_cores = 8, verbose = FALSE)
})
```

### Memory Optimization

For large-scale processing:

```r
# Process in chunks to manage memory
process_chunk <- function(file_chunk) {
  results <- lst_vat(file_chunk, verbose = FALSE)
  # Save results immediately
  saveRDS(results, sprintf("results_chunk_%s.rds", digest::digest(file_chunk)))
  rm(results)
  gc()
}

# Process 100 files in chunks of 10
chunk_size <- 10
for (i in seq(1, 100, chunk_size)) {
  files_chunk <- files[i:min(i+chunk_size-1, 100)]
  process_chunk(files_chunk)
}
```

## Continuous Integration / Docker

### Dockerfile Example

```dockerfile
FROM python:3.11-slim

# Install system dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    gcc \
    g++ \
    && rm -rf /var/lib/apt/lists/*

# Install Python dependencies
COPY inst/python/voice_analysis_python/requirements.txt /tmp/
RUN pip install --no-cache-dir -r /tmp/requirements.txt
RUN pip install cython>=0.29.0

# Install voice_analysis with Cython
COPY inst/python/voice_analysis_python /app/voice_analysis_python
WORKDIR /app/voice_analysis_python
RUN python setup_cython.py install

# Verify installation
RUN python -c "from voice_analysis.r_interface import check_cython_available; assert check_cython_available()"
```

### GitHub Actions Example

```yaml
name: Build and Test Voice Analysis

on: [push, pull_request]

jobs:
  test:
    runs-on: ${{ matrix.os }}
    strategy:
      matrix:
        os: [ubuntu-latest, macos-latest, windows-latest]
        python-version: ['3.8', '3.9', '3.10', '3.11']

    steps:
      - uses: actions/checkout@v3

      - name: Set up Python
        uses: actions/setup-python@v4
        with:
          python-version: ${{ matrix.python-version }}

      - name: Install dependencies
        run: |
          pip install -r inst/python/voice_analysis_python/requirements.txt
          pip install cython pytest

      - name: Build Cython extensions
        run: |
          cd inst/python/voice_analysis_python
          python setup_cython.py build_ext --inplace

      - name: Run tests
        run: |
          cd inst/python/voice_analysis_python
          pytest tests/
```

## References

### Optimization Techniques
- Numba documentation: https://numba.pydata.org/
- Cython documentation: https://cython.readthedocs.io/
- Platform-specific SIMD: ARM NEON, Intel AVX-512

### Scientific Background
- Tsanas, A. et al. (2011). J. Royal Society Interface 8(59):842-855
- Tsanas, A. (2012). D.Phil. thesis, University of Oxford

## Support

For build issues:
- Check: `voice_analysis_optimization_status()`
- Report: https://github.com/humlab-speech/superassp/issues
- Fallback: `install_voice_analysis(method = "pure")`

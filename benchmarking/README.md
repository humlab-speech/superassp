# Benchmarking Suite for superassp

This directory contains comprehensive benchmarking scripts for comparing performance across different DSP implementations in the superassp package.

## Overview

The superassp package provides multiple implementations of signal processing functions:

1. **Original wrassp** - C-based ASSP library implementations
2. **Parselmouth-optimized** - Memory-based Python/Parselmouth implementations (no file I/O)
3. **Python DSP** - Pure Python implementations (openSMILE, WORLD, etc.)

These benchmarks help measure:
- Performance improvements from optimization
- Speed vs accuracy tradeoffs
- Memory usage patterns
- Cross-implementation equivalence

## Benchmark Scripts

### `benchmark_suite.R`
**Comprehensive benchmark covering all major optimizations**

Compares:
- Parselmouth-optimized vs original Praat-calling functions
- SuperASP vs wrassp implementations

Functions benchmarked:
- Formant analysis (Burg method)
- FormantPath analysis
- Pitch tracking
- Intensity analysis
- Spectral moments
- ACF, RFC, RMS, ZCR analysis

Output:
- `benchmark_results.rds` - Full benchmark object with all metrics
- `benchmark_summary.csv` - Summary table with median times and memory

### `benchmark_suite_simple.R`
**Quick benchmark for rapid testing**

Simplified version focusing on:
- Key Parselmouth optimizations only
- Single test file
- Faster execution for CI/CD

### `benchmark_python_ssff.R`
**Python SSFF function benchmarks**

Specifically tests:
- Memory-based vs file-based Python DSP
- SSFF output file generation
- Audio loading strategies (av vs direct file I/O)

### `benchmark_python_memory_improvements.R`
**Memory optimization analysis**

Measures impact of:
- Direct memory processing (numpy arrays)
- Eliminating WAV conversion
- Bypassing temporary files
- av package for audio loading

### `benchmark_opensmile_slice_functions.R`
**openSMILE function benchmarks**

Tests:
- ComParE_2016 features
- GeMAPS features
- Memory-based vs file-based processing

## Running Benchmarks

### Prerequisites

```r
# Required packages
install.packages(c("bench", "dplyr", "tidyr"))

# Optional for full functionality
reticulate::py_install("praat-parselmouth")
reticulate::py_install("opensmile")
```

### Basic Usage

```r
# Run comprehensive benchmark suite
source("benchmarking/benchmark_suite.R")

# Quick test
source("benchmarking/benchmark_suite_simple.R")
```

### From Command Line

```bash
cd /path/to/superassp

# Full benchmark
Rscript benchmarking/benchmark_suite.R

# Quick benchmark
Rscript benchmarking/benchmark_suite_simple.R
```

## Results Interpretation

### Speedup Metrics

Speedup is calculated as:
```
speedup = median_time_original / median_time_optimized
```

Typical speedups:
- **Parselmouth optimizations**: 10-20x faster (no file I/O)
- **Memory-based Python**: 5-15x faster (no temporary files)
- **SuperASP vs wrassp**: Comparable (minor differences)

### Example Output

```
Benchmarking: Formant Analysis (Burg)
  File: sv1.wav
    Original: 245.3 ms, Optimized: 18.7 ms, Speedup: 13.12x
```

### Memory Usage

Memory is reported in MB. Optimized versions typically use:
- Less peak memory (no file copies)
- More consistent allocation (predictable numpy arrays)

## Benchmark Comparison Categories

### 1. Formant Estimation Procedures

Multiple methods for formant frequency estimation:

| Method | Implementation | Speed | Accuracy | Notes |
|--------|---------------|-------|----------|-------|
| Burg (Praat) | `praat_formant_burg()` | Medium | High | Traditional method |
| Burg (optimized) | `praat_formant_burg_opt()` | Fast | High | Memory-based |
| FormantPath | `praat_formantpath_burg_opt()` | Fast | Highest | Adaptive ceiling |
| ASSP (wrassp) | `forest()` | Fast | High | Alternative algorithm |

**Comparison criteria**:
- Processing speed
- Formant tracking accuracy
- Robustness to noise
- Parameter sensitivity

### 2. Pitch / F0 Estimation Procedures

Multiple fundamental frequency estimation algorithms:

| Method | Function | Algorithm Type | Strengths |
|--------|----------|---------------|-----------|
| Praat AC | `praat_pitch_opt()` | Autocorrelation | General purpose |
| Praat CC | `praat_pitch_opt()` | Cross-correlation | Robust |
| RAPT | `rapt()` | Autocorrelation | Low-frequency voices |
| REAPER | `reaper()` | Epoch detection | High quality |
| SWIPE | `swipe()` | Sawtooth waveform | No voicing decision |
| Dio (WORLD) | `dio()` | Instantaneous freq | Fast, real-time |
| Harvest (WORLD) | `harvest()` | Autocorrelation | High quality |
| SPICE | `spice()` | Neural network | Modern ML-based |

**Comparison criteria**:
- Gross pitch error rate (GPE)
- Fine pitch error (FPE)
- Voicing decision accuracy
- Processing speed
- Robustness to noise

### 3. Voice Quality Analysis

Clinical voice assessment functions:

| Measure | Function | Speed | Clinical Use |
|---------|----------|-------|--------------|
| Voice Report | `praat_voice_report_opt()` | Fast | Standard assessment |
| AVQI | `praat_avqi()` | Medium | Dysphonia severity |
| DSI | `praat_dsi()` | Medium | Dysphonia index |
| Tremor | `praat_voice_tremor()` | Slow | Tremor analysis |

## Benchmarking Best Practices

### 1. Test File Selection

Use diverse test files:
- Short files (~3s) for rapid iteration
- Medium files (~10s) for realistic workloads
- Long files (~30s+) for stress testing
- Various sample rates (8kHz, 16kHz, 44.1kHz, 48kHz)
- Different voice types (male, female, children)
- Clean and noisy recordings

### 2. Statistical Rigor

- Run multiple iterations (≥10) per benchmark
- Use median times (robust to outliers)
- Check for memory leaks (repeated runs)
- Monitor CPU and RAM usage
- Test on representative hardware

### 3. Fair Comparisons

When comparing implementations:
- Use identical input parameters
- Verify output equivalence first
- Account for algorithmic differences
- Document accuracy vs speed tradeoffs
- Consider real-world usage patterns

### 4. Interpretation Guidelines

**Speed improvements**:
- < 1.5x: Marginal, may not be perceptible
- 1.5-5x: Noticeable improvement
- 5-20x: Major optimization success
- > 20x: Dramatic (usually from eliminating I/O)

**When to optimize**:
- Batch processing large datasets
- Real-time/interactive applications
- Repeated analysis in workflows
- Memory-constrained environments

**When not to optimize**:
- One-time analysis
- Small datasets
- Accuracy is paramount
- Existing code works well

## Output Files

Benchmark scripts generate:

| File | Description | Use Case |
|------|-------------|----------|
| `benchmark_results.rds` | Full bench::mark objects | Detailed analysis, plotting |
| `benchmark_summary.csv` | Summary table | Quick review, reports |
| `*_comparison.csv` | Function-specific results | Targeted analysis |

## Integration with Vignettes

Results from these benchmarks are used in:

1. **`vignettes/benchmark_report.qmd`** - Main benchmarking report
2. **Performance optimization guides** - Best practices documentation
3. **Function documentation** - Performance notes in `?function` help

To regenerate vignette benchmarks:

```r
# Run benchmarks
source("benchmarking/benchmark_suite.R")

# Build vignette
quarto::quarto_render("vignettes/benchmark_report.qmd")
```

## Continuous Benchmarking

For package development:

1. Run `benchmark_suite_simple.R` before major commits
2. Full `benchmark_suite.R` for release candidates
3. Compare against baseline results
4. Document any performance regressions
5. Update vignettes if significant changes

##Contributing New Benchmarks

When adding new DSP functions:

1. Add comparison to appropriate benchmark script
2. Use consistent naming (`function_name` vs `function_name_opt`)
3. Document expected speedup range
4. Verify output equivalence
5. Update this README

### Template for New Benchmark

```r
list(
  name = "My New Function",
  orig = "original_function",
  opt = "original_function_opt",
  args = list(param1 = value1, param2 = value2)
)
```

## Troubleshooting

**Benchmark fails with "function not found"**:
- Check function is exported from package
- Verify Python modules installed for optimized versions
- Ensure test files exist

**Inconsistent results**:
- CPU throttling (laptop on battery)
- Background processes consuming resources
- Insufficient iterations

**Memory errors**:
- Large test files with limited RAM
- Memory leak in function being tested
- Try smaller test files or fewer iterations

## References

- [`bench` package documentation](https://bench.r-lib.org/)
- [Praat algorithms](https://www.fon.hum.uva.nl/praat/manual/Algorithms.html)
- [WORLD vocoder](https://github.com/mmorise/World)
- [openSMILE feature extraction](https://audeering.github.io/opensmile-python/)

## See Also

- Main benchmarking vignette: `vignette("benchmark_report", package = "superassp")`
- Performance tips: `vignette("performance", package = "superassp")`
- Function comparison tables: Package documentation

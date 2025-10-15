# SuperASP Testing and Benchmarking Suite

This directory contains comprehensive test suites and benchmarking tools for the SuperASP package.

## Overview

The testing infrastructure includes:

1. **Equivalence Tests** - Verify that different implementations produce equivalent results
2. **Benchmark Suite** - Measure and compare performance across implementations
3. **Quarto Report** - Generate detailed performance analysis reports

## Test Files

### Equivalence Tests

#### `testthat/test-equivalence-comprehensive.R`
Comprehensive equivalence testing that:
- Compares Parselmouth-optimized functions (`*_opt`) with original Praat-calling versions
- Compares SuperASP DSP functions with wrassp implementations
- Tests across all audio files in `signalfiles/`
- Validates:
  - Track names match
  - Number of frames are equivalent
  - Values are highly correlated (r > 0.95 for Parselmouth, r > 0.99 for SuperASP/wrassp)

#### `testthat/test-equivalence-simple.R`
Simplified equivalence tests for CI/CD:
- Tests core functionality with single representative file
- Faster execution for continuous integration
- Validates basic compatibility

### Running Equivalence Tests

```r
# From R console
library(testthat)
library(superassp)

# Run all tests
test_dir("tests/testthat")

# Run specific test
test_file("tests/testthat/test-equivalence-comprehensive.R")
```

```bash
# From command line
Rscript -e "testthat::test_dir('tests/testthat')"
```

## Benchmark Suite

### `benchmark_suite.R` (Full Version)
Comprehensive benchmarking including:
- Parselmouth-optimized vs original Praat-calling functions
- SuperASP vs wrassp implementations
- Multiple audio files with different characteristics
- Detailed timing and memory usage statistics

### `benchmark_suite_simple.R` (Simplified Version)
Focused benchmarking for:
- SuperASP vs wrassp core functions (acfana, rfcana, rmsana, zcrana)
- Quick performance comparison
- Reduced execution time

### Running Benchmarks

```bash
# From tests directory
cd tests
Rscript benchmark_suite_simple.R
```

### Benchmark Output

Benchmarks generate:
- `benchmark_results.rds` - Detailed results for analysis
- `benchmark_summary.csv` - CSV export for external tools
- Console output with real-time results

## Performance Report

### `../vignettes/benchmark_report.qmd`

Quarto document that generates:
- Visual performance comparisons
- Statistical analysis of speedups
- Memory usage analysis
- Detailed conclusions and recommendations

### Rendering the Report

```r
# From R
quarto::quarto_render("vignettes/benchmark_report.qmd")
```

```bash
# From command line
quarto render vignettes/benchmark_report.qmd
```

The report includes:
- Bar charts showing speedup factors
- Box plots comparing execution times
- Memory allocation comparisons
- Summary tables with statistics
- Conclusions and recommendations

## Test Data

Audio files are organized in `signalfiles/` subdirectories:
- `AVQI/input/` - Voice quality assessment files
- `DSI/input/` - Dysphonia severity index files
- `generated/` - Synthetic test signals

Different files test various scenarios:
- Short files (< 3 seconds) - sv1.wav, cs1.wav
- Medium files (3-10 seconds) - fh1.wav, im1.wav
- Stereo files - vowel14s_stereo.wav
- Synthetic signals - sine10m.wav

## Expected Results

### Parselmouth Optimization
- **Speedup**: 5-20x faster than original Praat-calling functions
- **Accuracy**: High correlation (r > 0.95) with original implementations
- **Memory**: Similar or slightly higher memory usage

### SuperASP vs wrassp
- **Performance**: Comparable, within 10-30% of each other
- **Accuracy**: Nearly identical results (r > 0.99)
- **Compatibility**: Same AsspDataObj structure

## Continuous Integration

For CI/CD pipelines, use:

```bash
# Quick equivalence check
Rscript -e "testthat::test_file('tests/testthat/test-equivalence-simple.R')"

# Optional: Quick benchmark
Rscript tests/benchmark_suite_simple.R
```

## Troubleshooting

### Parselmouth Tests Skipped
If you see "Parselmouth not available" messages:
```r
# Install parselmouth via pip
reticulate::py_install("praat-parselmouth")
```

### Memory Issues
For large-scale benchmarking:
- Reduce number of test files
- Decrease benchmark iterations
- Run benchmarks sequentially

### Test Failures
Common causes:
1. **Path issues**: Ensure working directory is `tests/`
2. **Missing dependencies**: Install `bench`, `testthat`, `wrassp`
3. **Python setup**: Configure reticulate properly for Parselmouth

## Adding New Tests

### Adding Equivalence Tests

```r
test_that("New function produces valid output", {
  result <- new_function("test_file.wav", toFile = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_true(length(names(result)) > 0)
})
```

### Adding Benchmark Tests

```r
bm <- bench::mark(
  implementation_a = function_a(test_file, toFile = FALSE),
  implementation_b = function_b(test_file, toFile = FALSE),
  iterations = 10,
  check = FALSE
)
```

## Contact

For issues or questions about testing:
- Check package documentation: `?superassp`
- Report issues: https://github.com/anthropics/superassp/issues
- Review test output carefully for clues

## References

- testthat documentation: https://testthat.r-lib.org/
- bench package: https://bench.r-lib.org/
- Quarto documentation: https://quarto.org/

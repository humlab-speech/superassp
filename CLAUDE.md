# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**superassp** is an R package that extends the wrassp package with speech signal processing capabilities from multiple frameworks (Praat, Python/SPTK, C++/ASSP, ESTK). All functions provide a unified wrassp-like interface, outputting SSFF files or AsspDataObj objects compatible with the emuR framework.

**Key Design Philosophy**: The package uses *Depends* (not *Imports*) for wrassp, so all wrassp functions are automatically available to users.

## Building and Testing

### Package Development Commands

```r
# Install package from source
devtools::install_github("humlab-speech/superassp", dependencies = "Imports")

# Or install locally
devtools::install()

# Load for development
devtools::load_all()

# Run all tests
devtools::test()

# Run specific test file
testthat::test_file("tests/testthat/test-sptk-pitch.R")

# Generate documentation
devtools::document()

# Check package
devtools::check()
```

### C/C++ Compilation

```bash
# Package compiles automatically during R CMD INSTALL
# Makevars in src/ controls compilation flags

# Clean and rebuild
R CMD INSTALL --preclean --no-multiarch --with-keep.source .
```

### Running Benchmarks

```r
# Comprehensive benchmark suite
source("benchmarking/benchmark_suite.R")

# Or from installed package
source(system.file("benchmarks", "run_benchmarks.R", package = "superassp"))
```

## Architecture

### Three-Layer DSP Processing Architecture

The package implements DSP functions through a three-layer architecture:

**Layer 1: Core DSP Implementations**
- **C/C++ Native**: ASSP library (`src/assp/`), ESTK (`src/ESTK/`), SPTK (`src/SPTK/`)
- **Praat Scripts**: Located in `inst/praat/`, called via Parselmouth
- **Python**: Via reticulate for pysptk, opensmile, and other Python libraries

**Layer 2: Low-Level Functions** (C++/Rcpp)
- Direct bindings to C++ implementations: `rapt_cpp()`, `swipe_cpp()`, `reaper_cpp()`, `dio_cpp()`, `estk_pda_cpp()`, `estk_pitchmark_cpp()`
- Require pre-loaded `AsspDataObj` as input
- Exposed in `R/RcppExports.R` (auto-generated)
- Fast but less user-friendly

**Layer 3: High-Level R Wrappers** (Recommended for users)
- Full-featured DSP functions: `rapt()`, `swipe()`, `reaper()`, `dio()`, `ksvfo()`, `mhspitch()`, etc.
- Standard interface with parameters like `listOfFiles`, `toFile`, `beginTime`, `endTime`, `minF`, `maxF`, `windowShift`
- Handle any media format via av package (WAV, MP3, MP4, video files)
- Automatic parallel processing for batch operations
- Support both in-memory (`toFile=FALSE`) and file output (`toFile=TRUE`)

### Media Processing Pipeline

**Modern Load-and-Process Pattern** (Preferred):
1. `av_to_asspDataObj()` - Load any media format into memory
2. `processMediaFiles_LoadAndProcess()` - Unified processing with parallel support
3. `.External("performAsspMemory", ...)` - In-memory DSP via C interface

**Key Functions** (`R/av_helpers.R`):
- `av_to_asspDataObj()`: Convert any media file to AsspDataObj using av package
- `processMediaFiles_LoadAndProcess()`: Batch processing with automatic parallelization
- Handles time windowing, format conversion, and parallel processing internally

### AsspDataObj Data Structure

The central data structure throughout the package:

```r
# Structure
obj <- list(
  audio = matrix(samples),  # INT16 sample data
  # or other tracks like f0, formants, etc.
)
attr(obj, "sampleRate") <- numeric(sample_rate)
attr(obj, "startTime") <- numeric(start_time)
attr(obj, "startRecord") <- integer(start_frame)
attr(obj, "endRecord") <- integer(end_frame)
attr(obj, "trackFormats") <- character(format_vector)
class(obj) <- "AsspDataObj"
```

### Parallel Processing

All DSP functions support automatic parallelization:
- **Auto-enabled**: For 2+ files
- **Platform-aware**: Fork-based (Unix/Mac) vs socket cluster (Windows)
- **Thread-safe**: Each worker processes independently via `av_to_asspDataObj()` → `performAsspMemory()`
- Controlled by `parallel` and `n_cores` parameters (passed through `...`)

## Adding New DSP Functions

### For C++ Native Implementations (SPTK/ESTK style)

1. **Implement C++ function** in `src/`:
   - Accept `AsspDataObj` via Rcpp
   - Return list with results (f0, times, sample_rate, n_frames, etc.)
   - Use `// [[Rcpp::export]]` to expose to R

2. **Create R wrapper** in `R/superassp_*.R`:
   - Follow pattern from `R/superassp_sptk_pitch.R` (rapt, swipe, reaper, dio)
   - Load audio via `av_to_asspDataObj()`
   - Call C++ function
   - Convert results to AsspDataObj with proper attributes
   - Handle file output with `write.AsspDataObj()` if `toFile=TRUE`
   - Support batch processing with progress bars

3. **Add tests** in `tests/testthat/test-*.R`:
   - Basic functionality with single file
   - Custom parameters (F0 range, time windowing)
   - Batch processing
   - File I/O modes (`toFile=TRUE` and `FALSE`)
   - Non-WAV media formats
   - Error handling

### For ASSP Library Functions

Use the unified `processMediaFiles_LoadAndProcess()` helper:

```r
my_dsp_function <- function(listOfFiles, beginTime = 0.0, endTime = 0.0,
                           param1 = default1, toFile = TRUE,
                           outputDirectory = NULL, verbose = TRUE) {

  # Setup
  explicitExt <- "ext"
  newTracknames <- c("track1", "track2")
  nativeFiletypes <- c("wav", "au")

  # Time parameter normalization
  beginTime <- if(is.null(beginTime)) 0.0 else beginTime
  endTime <- if(is.null(endTime)) 0.0 else endTime
  n_files <- length(listOfFiles)
  beginTime <- fast_recycle_times(beginTime, n_files)
  endTime <- fast_recycle_times(endTime, n_files)

  # Use unified processor
  result <- processMediaFiles_LoadAndProcess(
    listOfFiles = listOfFiles,
    beginTime = beginTime,
    endTime = endTime,
    nativeFiletypes = nativeFiletypes,
    fname = "assp_function_name",  # ASSP C function name
    toFile = toFile,
    verbose = verbose,
    param1 = param1,  # Pass DSP parameters
    explicitExt = explicitExt,
    outputDirectory = outputDirectory
  )

  externalRes <- result$externalRes

  # Rename tracks if needed
  if(!toFile && !is.null(newTracknames)) {
    externalRes <- fast_rename_tracks(externalRes, newTracknames)
  }

  # Simplify single file output
  if(n_files == 1) externalRes <- externalRes[[1]]

  return(externalRes)
}

# Set function attributes
attr(my_dsp_function, "ext") <- "ext"
attr(my_dsp_function, "tracks") <- c("track1", "track2")
attr(my_dsp_function, "outputType") <- "SSFF"
attr(my_dsp_function, "nativeFiletypes") <- c("wav", "au")
```

### For Praat Functions

1. **Create Praat script** in `inst/praat/`:
   - Accept input file path and output CSV path as arguments
   - Compute DSP analysis
   - Write results to CSV table
   - Return table file path

2. **Create R wrapper**:
   - Use `praat_formant_burg` as template for signal tracks
   - Call Praat script via reticulate/Parselmouth
   - Read CSV output
   - Convert to AsspDataObj with proper track names
   - Handle SSFF file writing

3. **Keep praat_ prefix** for consistency

## Common Patterns

### Testing Patterns

```r
test_that("function works with single file", {
  skip_if_not_installed("superassp")

  test_wav <- system.file("samples", "sustained", "a1.wav",
                         package = "superassp")
  skip_if(test_wav == "", "Test file not found")

  result <- my_function(test_wav, toFile = FALSE, verbose = FALSE)

  expect_s3_class(result, "AsspDataObj")
  expect_true("track_name" %in% names(result))
  # Additional assertions...
})
```

### Rcpp Helper Functions

Fast utility functions in `src/dsp_helpers.cpp`:
- `fast_file_ext()`: Extract file extensions
- `fast_recycle_times()`: Recycle time parameters efficiently
- `fast_rename_tracks()`: Batch rename tracks in AsspDataObj list
- `fast_is_native()`, `fast_is_lossless()`: Format checks

## Important Conventions

### Naming
- **DSP functions**: Follow wrassp conventions (e.g., `mhspitch`, `ksvfo`, `forest`)
- **SPTK wrappers**: Use lowercase names (`rapt`, `swipe`, `reaper`, `dio`)
- **Low-level C++**: Add `_cpp` suffix (`rapt_cpp`, `swipe_cpp`)
- **Praat functions**: Use `praat_` prefix
- **Helper functions**: Clear descriptive names

### Parameters
Standard DSP function parameters:
- `listOfFiles`: Input file path(s)
- `beginTime`, `endTime`: Time windowing (seconds)
- `windowShift`: Frame shift (milliseconds)
- `minF`, `maxF`: F0 range for pitch tracking
- `toFile`: Write to file (TRUE) or return object (FALSE)
- `explicitExt`: Output file extension
- `outputDirectory`: Where to write files
- `verbose`: Progress messages

### File Locations
- `R/superassp_*.R`: Main DSP function implementations
- `R/av_helpers.R`: Media loading and processing helpers
- `src/*.cpp`: C++ implementations and Rcpp bindings
- `src/assp/`: ASSP C library
- `src/SPTK/`: SPTK submodule
- `src/ESTK/`: Edinburgh Speech Tools submodule
- `tests/testthat/test-*.R`: Test files
- `benchmarking/`: Benchmark scripts

### SSFF Format
All track-based outputs use SSFF (Simple Signal File Format):
- Binary format for time-series data
- Written via `write.AsspDataObj()`
- Track names should be descriptive (e.g., "pitch[Hz]", "fm", "bw")

## Performance Considerations

- Prefer C++ implementations over Python (2-3x faster)
- Use `processMediaFiles_LoadAndProcess()` for automatic parallelization
- Avoid file I/O when possible - process in memory
- SPTK wrappers (`rapt`, `swipe`, etc.) recommended for users
- Low-level `_cpp` functions for advanced use cases only

## Dependencies

- **R packages**: wrassp (Depends), av, reticulate, Rcpp, parallel, cli, rlang
- **System**: C++11 compiler
- **Optional**: Praat (for Praat functions), Python with pysptk/parselmouth
- **Submodules**: SPTK (src/SPTK), ESTK (src/ESTK)

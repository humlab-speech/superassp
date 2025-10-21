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
  - **Automatic fallback**: Tries av (FFmpeg) first, falls back to wrassp for niche formats
  - Supports: wav, mp3, mp4, flac, ogg, aac, opus (via av) + au, kay, nist, nsp, csre, ssff (via wrassp)
  - Time windowing and resampling supported (av only; wrassp warns if resampling requested)
- `processMediaFiles_LoadAndProcess()`: Batch processing with automatic parallelization
- Handles time windowing, format conversion, and parallel processing internally

### S7 AVAudio Class (v0.6.0)

**In-Memory Audio Processing**: The AVAudio S7 class enables efficient in-memory audio workflows:

```r
# Load audio into memory
audio <- read_avaudio("speech.wav", sample_rate = 16000, channels = 1)

# S7 automatic dispatch - works with any DSP function
f0 <- trk_rapt(audio, toFile = FALSE)
formants <- trk_forest(audio, toFile = FALSE)
```

**Key Features**:
- Preprocessing: Automatic resampling, channel mixing, time windowing
- Zero file I/O: Process entirely in memory
- S7 dispatch: All `trk_*` and `lst_*` functions accept AVAudio objects
- Implementation: `R/s7_avaudio.R`, `R/s7_methods.R`

**S7 Method Registration Pattern**:
```r
# Define generic (if not exists)
trk_function <- new_generic("trk_function", "x")

# Character method (file path)
method(trk_function, class_character) <- function(x, ...) {
  # Original implementation
}

# AVAudio method (in-memory)
method(trk_function, AVAudio) <- function(x, ...) {
  # Convert AVAudio to temp file or process directly
}
```

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

### Function Naming Conventions

All DSP functions must follow these prefixes:

- **`trk_*`**: Time-series tracks (e.g., `trk_rapt`, `trk_forest`, `trk_swiftf0`)
  - Return signal tracks that follow the audio waveform
  - Output: AsspDataObj with time-aligned tracks (F0, formants, energy, etc.)
  - Examples: pitch tracking, formant tracking, energy analysis

- **`lst_*`**: Summary statistics (e.g., `lst_voice_sauce`, `lst_vat`, `lst_covarep_vq`)
  - Return aggregate measures that summarize audio properties
  - Output: Data frame or list with scalar/vector values
  - Examples: jitter, shimmer, HNR, spectral features

### For Python Module Integrations

**Pattern for Python-based DSP functions**:

1. **Installation helpers** (`R/install_voice_analysis.R` or similar):
```r
install_module <- function(envname = NULL, method = "auto", ...) {
  reticulate::py_install("module-name", envname = envname, ...)
}

module_available <- function() {
  reticulate::py_module_available("module_name")
}

module_info <- function() {
  # Return module specifications
}
```

2. **Main function** (follow existing patterns):
   - For tracks: `trk_swiftf0()`, `trk_crepe()`, `trk_pyin()`
   - For summaries: `lst_voice_sauce()`, `lst_vat()`, `lst_covarep_vq()`

3. **Audio loading**:
```r
# Load via av package
audio_data <- av::read_audio_bin(file, channels = 1)
sample_rate <- attr(audio_data, "sample_rate")

# Convert to numpy array for Python
audio_float <- as.numeric(audio_data) / 2147483647.0  # INT32_MAX
np <- reticulate::import("numpy")
audio_array <- np$array(audio_float, dtype = "float32")
```

4. **Return format**:
   - Tracks: Convert Python output to AsspDataObj with proper attributes
   - Summaries: Return data frame with descriptive column names

**Examples**:
- `trk_swiftf0()`: Deep learning pitch tracker (R/ssff_python_swiftf0.R:112)
- `lst_voice_sauce()`: VoiceSauce voice quality measures
- `lst_vat()`: Voice Analysis Toolbox (132 dysphonia measures)

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
   - S7 AVAudio dispatch (if applicable)
   - Error handling

4. **Register S7 methods** in `R/s7_methods.R` (for automatic AVAudio dispatch):
```r
# Import S7
library(S7)

# Create or extend generic
trk_myfunction <- new_generic("trk_myfunction", "x")

# Character method (file paths)
method(trk_myfunction, class_character) <- function(x, ...) {
  # Your existing implementation
}

# AVAudio method (in-memory audio)
method(trk_myfunction, AVAudio) <- function(x, toFile = FALSE,
                                            explicitExt = "ext", ...) {
  if (toFile) {
    # Create temp file from AVAudio
    temp_file <- write_avaudio_temp(x, ext = "wav")
    on.exit(unlink(temp_file), add = TRUE)
    result <- trk_myfunction(temp_file, toFile = toFile,
                             explicitExt = explicitExt, ...)
    return(result)
  } else {
    # Process directly from AVAudio
    temp_file <- write_avaudio_temp(x, ext = "wav")
    on.exit(unlink(temp_file), add = TRUE)
    result <- trk_myfunction(temp_file, toFile = FALSE, ...)
    return(result)
  }
}
```

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
- **Track functions**: Use `trk_` prefix (e.g., `trk_rapt`, `trk_swiftf0`, `trk_forest`, `trk_mhspitch`)
- **Summary functions**: Use `lst_` prefix (e.g., `lst_voice_sauce`, `lst_vat`, `lst_covarep_vq`)
- **Low-level C++**: Add `_cpp` suffix (`rapt_cpp`, `swipe_cpp`, `estk_pda_cpp`)
- **Praat functions**: Use `praat_` prefix (legacy, but keep for existing functions)
- **Installation helpers**: `install_*`, `*_available`, `*_info` patterns
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
- `R/ssff_*.R`: DSP function implementations (track-based)
  - `R/ssff_c_assp_*.R`: ASSP C library wrappers
  - `R/ssff_python_*.R`: Python-based implementations
- `R/list_*.R`: Summary statistic functions
- `R/superassp_*.R`: Legacy DSP implementations (being migrated to ssff_*)
- `R/av_helpers.R`: Media loading and processing helpers
- `R/s7_avaudio.R`: S7 AVAudio class definition
- `R/s7_methods.R`: S7 method registrations for DSP functions
- `R/install_*.R`: Python module installation helpers
- `src/*.cpp`: C++ implementations and Rcpp bindings
- `src/assp/`: ASSP C library
- `src/SPTK/`: SPTK submodule (pitch tracking, MFCC, etc.)
- `src/ESTK/`: Edinburgh Speech Tools submodule
- `inst/python/`: Python modules (voice_analysis_python, covarep_python, etc.)
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

- **R packages**: wrassp (Depends), av, reticulate, Rcpp, S7, parallel, cli, rlang
- **System**: C++11 compiler
- **Optional Python modules**:
  - `swift-f0`: Deep learning pitch tracker (install via `install_swiftf0()`)
  - `pysptk`, `parselmouth`: Alternative pitch/formant implementations
  - `voice_analysis_python`: 132 dysphonia measures (install via `install_voice_analysis()`)
  - Others as needed for specific functions
- **Submodules**: SPTK (src/SPTK), ESTK (src/ESTK)

## Key Recent Additions (v0.6.0+)

### S7 AVAudio Class (v0.6.0)
- In-memory audio processing with automatic dispatch
- All DSP functions accept both file paths and AVAudio objects
- See: `R/s7_avaudio.R:109`, `R/s7_methods.R:1`

### Swift-F0 Deep Learning Pitch Tracker
- Fast CNN-based F0 detection (~90-130ms for 3s audio)
- Installation: `install_swiftf0()`
- Function: `trk_swiftf0()` (R/ssff_python_swiftf0.R:112)
- Frequency range: 46.875-2093.75 Hz (G1 to C7)

### Automatic Format Fallback
- `av_to_asspDataObj()` now tries av first, falls back to wrassp
- Handles niche formats: au, kay, nist, nsp, csre, ssff
- Graceful degradation with warnings (R/av_helpers.R:109)

### Python Module Integrations
- VoiceSauce: `lst_voice_sauce()` - 34 voice quality measures
- Voice Analysis Toolbox: `lst_vat()` - 132 dysphonia measures
- COVAREP: `lst_covarep_vq()` - Voice quality parameters
- Each has `install_*()`, `*_available()`, `*_info()` helpers

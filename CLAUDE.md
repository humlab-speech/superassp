# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

**superassp** is an R package providing comprehensive speech signal processing capabilities from multiple frameworks (Praat, Python/SPTK, C++/ASSP, ESTK, OpenSMILE). All functions provide a unified interface, outputting SSFF files or AsspDataObj objects compatible with the emuR framework.

**Key Design Philosophy**:
- **Self-contained**: Complete DSP implementations from ASSP, SPTK, ESTK, OpenSMILE (C++), and Python frameworks
- **Independent**: No required external dependencies beyond R and C++ compiler
- **wrassp is NOT required**: superassp includes its own ASSP C library and can operate completely independently
- **Universal media support**: All functions accept any media format (WAV, MP3, MP4, video) via av package
- **Performance-optimized**: Native C++ implementations preferred for speed, Python available for specialized algorithms

**Note on wrassp**: The wrassp package is a separate R package that also provides ASSP-based functions. While superassp and wrassp can be used together, **superassp does not require wrassp** and includes its own complete ASSP library implementation. All superassp DSP functions work independently.

## Quick Reference - Most Common Tasks

```r
# ============================================
# DEVELOPMENT WORKFLOW
# ============================================

# Load package for development
devtools::load_all()

# After changing C++ code
Rcpp::compileAttributes()  # Update R exports
devtools::document()       # Regenerate docs
devtools::load_all()       # Reload

# After changing R code or docs
devtools::document()       # Regenerate docs
devtools::load_all()       # Reload

# Test and check
devtools::test()           # Run all tests
testthat::test_file("tests/testthat/test-filename.R")  # Single test
devtools::check()          # Full package check

# Before committing
devtools::document()       # ALWAYS regenerate docs
devtools::test()          # ALWAYS run tests

# ============================================
# GIT WORKFLOWS
# ============================================

# Initialize submodules (first time)
git submodule update --init --recursive

# Update submodules to latest
git submodule update --remote --recursive

# Commit changes
devtools::document()  # In R
devtools::test()      # In R
git add R/file.R man/file.Rd
git commit -m "feat: Add new DSP function"

# ============================================
# ADDING NEW DSP FUNCTIONS
# ============================================

# C++ function (recommended for performance)
# 1. Create src/myfunction.cpp with // [[Rcpp::export]]
# 2. Rscript -e "Rcpp::compileAttributes()"
# 3. Create R/ssff_cpp_myfunction.R wrapper
# 4. devtools::document()
# 5. Add tests in tests/testthat/test-myfunction.R

# Python function (for specialized algorithms)
# 1. Create inst/python/mymodule/ or inst/python/myscript.py
# 2. Create R/install_mymodule.R helpers
# 3. Create R/ssff_python_myfunction.R wrapper
# 4. Use av::read_audio_bin() for audio (NOT librosa)
# 5. devtools::document()
# 6. Add tests with module availability checks
```

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

# Generate documentation from roxygen2 comments
devtools::document()

# Check package (runs tests, checks documentation, etc.)
devtools::check()

# Build package
devtools::build()
```

**Documentation Notes**:
- All functions use roxygen2 comments (`#'`) for documentation
- Run `devtools::document()` after modifying function signatures or documentation
- Documentation is exported to `man/*.Rd` files automatically
- NAMESPACE is auto-generated - never edit manually
- **When to regenerate docs**:
  - After adding/removing function parameters
  - After changing function signatures in `// [[Rcpp::export]]`
  - After modifying roxygen2 `@param`, `@return`, `@examples` comments
  - After running `Rcpp::compileAttributes()` (updates `R/RcppExports.R`)
  - Always regenerate docs before committing changes

### C/C++ Compilation

```bash
# Package compiles automatically during R CMD INSTALL
# Makevars in src/ controls compilation flags

# Clean and rebuild
R CMD INSTALL --preclean --no-multiarch --with-keep.source .

# Rebuild only C++ code without reinstalling
Rscript -e "Rcpp::compileAttributes()"
R CMD SHLIB src/*.cpp
```

**Important Compilation Notes**:
- `src/Makevars` defines all compilation settings (includes, source lists, flags)
- SPTK, ESTK, and ASSP libraries compile as part of package build
- After modifying C++ code, run `Rcpp::compileAttributes()` to update `R/RcppExports.R`
- Submodules (SPTK, ESTK) are included as git submodules - update with `git submodule update --init --recursive`

**Git Submodule Management**:
```bash
# Initialize all submodules (required for first build)
git submodule update --init --recursive

# Update all submodules to latest commits
git submodule update --remote --recursive

# Check submodule status
git submodule status

# Current submodules:
# - src/SPTK: Speech Signal Processing Toolkit (pitch, MFCC, spectral analysis)
# - src/ESTK: Edinburgh Speech Tools (pitch detection, pitchmarking)
# - src/tcl-snack: Snack Sound Toolkit (reference implementations)
# - inst/onnx/swift-f0: Swift-F0 deep learning pitch tracker
# - inst/python/DeepFormants: Deep learning formant tracking
```

**⚠️ CRITICAL: Never modify submodule code directly**
- Submodules point to external repositories (SPTK, ESTK)
- Changes must be made in the upstream repository
- Update submodule commit references after upstream changes
- Test thoroughly after submodule updates (may affect DSP behavior)
- **Exception**: OpenSMILE is NOT a submodule - it's bundled in `src/opensmile/` and can be modified directly

### Running Benchmarks

```r
# Comprehensive benchmark suite
source("inst/benchmarking/r/benchmark_suite.R")

# Or from installed package
source(system.file("benchmarking", "r", "run_benchmarks.R", package = "superassp"))
```

## Common Development Workflows

### Adding a new C++ DSP function

1. **Implement C++ function** in `src/yourfunction.cpp`:
```cpp
#include <Rcpp.h>
// [[Rcpp::export]]
Rcpp::List yourfunction_cpp(Rcpp::List audio_obj, double param1) {
  // Implementation
  return Rcpp::List::create(
    Rcpp::Named("result") = result,
    Rcpp::Named("sample_rate") = sample_rate
  );
}
```

2. **Update C++ exports**: `Rscript -e "Rcpp::compileAttributes()"`
3. **Create R wrapper** in `R/ssff_cpp_yourfunction.R`
4. **Add roxygen2 documentation** with `#'` comments
5. **Regenerate docs**: `devtools::document()`
6. **Update `src/Makevars`** if adding new source files
7. **Add tests** in `tests/testthat/test-yourfunction.R`
8. **Test locally**: `devtools::test()` and `devtools::check()`

### Adding a new Python DSP function

1. **Create Python script** (if needed) in `inst/python/yourscript.py`
2. **Create installation helper** in `R/install_yourmodule.R`
3. **Create R wrapper** in `R/ssff_python_yourfunction.R` or `R/list_python_yourfunction.R`
4. **Use `av::read_audio_bin()`** for audio loading (NOT librosa)
5. **Add roxygen2 documentation**
6. **Regenerate docs**: `devtools::document()`
7. **Add tests** including Python module availability checks
8. **Test with multiple media formats** (WAV, MP3, MP4)

### Modifying existing DSP functions

1. **Read the function** to understand current implementation
2. **Check tests** in `tests/testthat/test-*.R` for expected behavior
3. **Make changes** preserving function signature if possible
4. **Update roxygen2 docs** if parameters/behavior changed
5. **Regenerate docs**: `devtools::document()`
6. **Run tests**: `devtools::test()`
7. **Run benchmarks** (if performance-critical)
8. **Update NEWS.md** with changes

### Git commit workflow

```bash
# After making changes
devtools::document()  # In R - regenerate documentation
devtools::test()      # Run tests

# Stage changes
git add R/modified_file.R man/modified_function.Rd

# Commit with conventional commit message
git commit -m "feat: Add new DSP function for pitch tracking"
# or
git commit -m "fix: Correct time windowing in trk_rapt"
# or
git commit -m "docs: Update CLAUDE.md with development workflows"
```

## Troubleshooting

### C++ compilation errors

**Problem**: `undefined reference to SPTK::...`
```bash
# Solution: Ensure submodules are initialized
git submodule update --init --recursive
# Then rebuild
devtools::clean_dll()
devtools::load_all()
```

**Problem**: `Rcpp function not found`
```r
# Solution: Regenerate Rcpp exports
Rcpp::compileAttributes()
devtools::document()
devtools::load_all()
```

**Problem**: Compilation fails with missing headers
```bash
# Solution: Check src/Makevars includes
# Ensure PKG_CPPFLAGS includes all necessary -I flags
```

### Documentation issues

**Problem**: Function not appearing in NAMESPACE
```r
# Solution: Add @export to roxygen2 comments
#' @export
devtools::document()
```

**Problem**: Changes to C++ function not reflected in docs
```r
# Solution: Full regeneration workflow
Rcpp::compileAttributes()  # Updates R/RcppExports.R
devtools::document()       # Updates man/*.Rd and NAMESPACE
```

### Python integration issues

**Problem**: `Python module not found`
```r
# Solution: Check reticulate configuration
reticulate::py_config()
# Install module
install_yourmodule()  # Use package installation helper
```

**Problem**: `av::read_audio_bin` not working
```r
# Solution: Check av package installation
install.packages("av")
# Check FFmpeg availability
av::av_video_info(test_file)  # Should work if FFmpeg is available
```

### Testing issues

**Problem**: Tests fail with "file not found"
```r
# Solution: Use system.file() for test data
test_file <- system.file("samples", "sustained", "a1.wav", package = "superassp")
skip_if(test_file == "", "Test file not found")
```

**Problem**: Tests timeout on parallel processing
```r
# Solution: Disable parallel for tests
result <- trk_rapt(files, toFile = FALSE, verbose = FALSE, parallel = FALSE)
```

## Architecture

### Function Categories and Discovery

The package provides **75+ DSP functions** organized into 9 main categories:

1. **Pitch/F0 Tracking** (17 functions): `trk_rapt()`, `trk_swipe()`, `trk_dio()`, `trk_harvest()`, `trk_reaper()`, `trk_swiftf0()`, `trk_crepe()`, `trk_sacc()`, `trk_pyin()`, `trk_yin()`, `trk_yaapt()`, etc.
2. **Formant Analysis** (7 functions): `trk_forest()`, `trk_deepformants()`, `trk_formants_tvwlp()`, `trk_formantp()`, etc.
3. **Spectral Analysis** (6 functions): `trk_dftSpectrum()`, `trk_cssSpectrum()`, `trk_lpsSpectrum()`, `trk_cepstrum()`, etc.
4. **Energy & Amplitude** (4 functions): `trk_rmsana()`, `trk_zcrana()`, `trk_acfana()`, `trk_intensityp()`
5. **Voice Quality** (10 functions): `lst_vat()` (132 measures), `lst_voice_sauce()` (40+ params), `trk_brouhaha()`, `trk_creak_union()`, etc.
6. **Prosody & Intonation** (2 functions): `lst_dysprosody()` (193 features), `lst_voxit()` (11 measures)
7. **Source-Filter Decomposition** (3 functions): `trk_gfmiaif()`, `trk_covarep_iaif()`, `trk_excite()`
8. **OpenSMILE Feature Sets** (5 groups): `lst_GeMAPS()` (62 features), `lst_eGeMAPS()` (88 features), `lst_emobase()`, `lst_ComParE_2016()`
9. **Acoustic Features** (3 functions): `trk_mfcc()`, `trk_lp_analysis()`, `trk_npy_import()`

**For detailed function descriptions and usage recommendations, see `PKGDOWN_FUNCTION_GROUPING.md`.**

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

### Parselmouth In-Memory Processing (v0.8.7+)

**CRITICAL: All parselmouth-based functions use in-memory processing - NO temporary files!**

**Helper Functions** (`R/parselmouth_helpers.R`):

```r
# Convert av audio data to parselmouth Sound object
sound <- av_to_parselmouth_sound(audio_data)

# Complete workflow: load file → Sound object
sound <- av_load_for_parselmouth(
  file_path = "audio.mp3",
  start_time = 1.0,
  end_time = 3.0,
  channels = 1,
  target_sample_rate = 16000
)

# Check if parselmouth is available
if (parselmouth_available()) {
  # Process audio
}
```

**Architecture Pattern 1**: Python creates Sound from numpy (most functions)
```r
# R side
audio_data <- av_load_for_python(file_path, start_time, end_time)

# Python side
sound = pm.Sound(audio_data['audio_np'], sampling_frequency=audio_data['sample_rate'])
```

**Architecture Pattern 2**: R creates Sound directly (newer functions)
```r
# R side
sound <- av_load_for_parselmouth(file_path, start_time, end_time)

# Python side
def function_from_sound(sound, params...):
    # sound is already a parselmouth.Sound object!
    result = process(sound, params)
    return result
```

**Functions Using Pattern 1** (6 functions):
- `lst_avqip()` - AVQI voice quality index
- `lst_dsip()` - Dysphonia Severity Index
- `lst_voice_reportp()` - Praat voice report
- `lst_voice_tremorp()` - Voice tremor analysis
- `trk_sacc()` - SAcC pitch tracking

**Functions Using Pattern 2** (2 functions):
- `lst_dysprosody()` - 193 prosodic features
- `trk_pitchp()` - Praat pitch tracking (multiple methods)

**Adding New Parselmouth Functions**:

For Python scripts in `inst/python/`:
```python
def function_from_sound(sound, param1, param2, ...):
    """
    Process audio from parselmouth Sound object.

    Parameters:
    -----------
    sound : parselmouth.Sound
        Sound object (not a file path)
    """
    # sound is already loaded - use directly!
    snd = sound
    # ... processing ...
    return result
```

For R wrappers:
```r
# Load audio and convert to Sound
sound <- av_load_for_parselmouth(
  file_path = file_path,
  start_time = if (bt > 0) bt else NULL,
  end_time = if (et > 0) et else NULL,
  channels = 1
)

# Call Python function
result <- py$function_from_sound(sound, param1, param2)
```

**Migration Status**: See `PYTHON_INMEMORY_MIGRATION_PLAN.md` for remaining functions.

**Performance**: 38% faster than file-based approach (eliminates disk I/O).

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
  - **REQUIRED**: All `trk_*` functions MUST have:
    - `toFile` parameter (default: FALSE for backward compatibility)
    - `explicitExt` parameter specifying output file extension
    - `outputDirectory` parameter (default: NULL = same as input)
    - Function attributes: `ext`, `tracks`, `outputType`, `nativeFiletypes`
    - Parameter and attribute consistency: `explicitExt` default must match `attr(*, "ext")`

- **`lst_*`**: Summary statistics (e.g., `lst_voice_sauce`, `lst_vat`, `lst_covarep_vq`)
  - Return aggregate measures that summarize audio properties
  - Output: Data frame or list with scalar/vector values
  - Examples: jitter, shimmer, HNR, spectral features
  - Extension attributes optional (most return data frames, not SSFF files)

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


### For lst_* Functions with JSON Track Format (JSTF)

**NEW (v0.10.0)**: `lst_*` functions can now write time-sliced results to JSON Track Format files for efficient multi-slice storage.

**Pattern**: JSON-based file output for list-producing DSP functions

1. **Add toFile parameters**:
```r
lst_function <- function(listOfFiles,
                        beginTime = 0.0,
                        endTime = 0.0,
                        toFile = FALSE,              # NEW
                        explicitExt = "ext",         # NEW
                        outputDirectory = NULL,      # NEW
                        verbose = TRUE,
                        ...) {
  
  # Process audio
  audio_data <- av::read_audio_bin(file, ...)
  results <- your_dsp_processing(audio_data)
  
  if (toFile) {
    # Create JsonTrackObj
    json_obj <- create_json_track_obj(
      results = results,
      function_name = "lst_function",
      file_path = file,
      sample_rate = attr(audio_data, "sample_rate"),
      audio_duration = length(audio_data) / sample_rate,
      beginTime = beginTime,
      endTime = endTime,
      parameters = list(...)
    )
    
    # Write to file
    output_path <- file.path(
      outputDirectory %||% dirname(file),
      paste0(tools::file_path_sans_ext(basename(file)), ".", explicitExt)
    )
    
    write_json_track(json_obj, output_path)
    return(invisible(output_path))
  }
  
  return(results)  # In-memory mode
}

# Set function attributes
attr(lst_function, "ext") <- "ext"
attr(lst_function, "outputType") <- "JSTF"
attr(lst_function, "format") <- "JSON"
```

2. **Register extension** in `inst/extdata/json_extensions.csv`:
```csv
function,extension,description,fields,format
lst_function,ext,Description of output,N,JSTF
```

3. **JSON Track Format benefits**:
   - **Efficient**: Avoids field name duplication across slices (~99% reduction)
   - **Fast reading**: RcppSimdJson provides 3x faster parsing than jsonlite
   - **Human-readable**: JSON format is text-based and debuggable
   - **Flexible**: Supports complex nested structures (lists, matrices, vectors)
   - **Compatible**: Converts to data.frame/tibble like AsspDataObj

4. **Example JSTF file structure**:
```json
{
  "format": "JSTF",
  "version": "1.0",
  "function": "lst_vat",
  "file_path": "audio.wav",
  "sample_rate": 16000,
  "audio_duration": 5.0,
  "field_schema": {
    "jitter": "numeric",
    "shimmer": "numeric",
    "hnr": "numeric"
  },
  "slices": [
    {
      "begin_time": 0.0,
      "end_time": 1.0,
      "values": [85.3, 4.2, 15.7]
    },
    {
      "begin_time": 1.0,
      "end_time": 2.0,
      "values": [88.1, 3.9, 16.2]
    }
  ]
}
```

5. **Usage pattern**:
```r
# Write to file
lst_vat("audio.wav", toFile = TRUE)  # Creates audio.vat

# Read back transparently
track <- read_track("audio.vat")     # Auto-detects JSTF

# Convert to data.frame
df <- as.data.frame(track)
#   begin_time end_time jitter shimmer  hnr
# 1        0.0      1.0   85.3     4.2 15.7
# 2        1.0      2.0   88.1     3.9 16.2

# Or use with tibble
library(dplyr)
as_tibble(track) %>% filter(begin_time > 0.5)
```

6. **Registered JSTF extensions**:
   - `.vat` - Voice Analysis Toolbox (132 measures)
   - `.vsj` - VoiceSauce voice quality (40+ params)
   - `.dyp` - Dysprosody features (193 features)
   - `.vxt` - Voxit measures (11 features)
   - `.gem` - GeMAPS features (62 features)
   - `.egm` - eGeMAPS features (88 features)
   - `.emb` - emobase features
   - `.cmp` - ComParE 2016 features
   - `.cvq` - COVAREP voice quality
   - `.avq` - AVQI index
   - `.dsi` - Dysphonia Severity Index
   - `.vrp` - Praat voice report
   - `.vtr` - Voice tremor analysis
   - `.phn` - Phonological posteriors

7. **Key functions**:
   - `create_json_track_obj()` - Create JsonTrackObj from results
   - `write_json_track()` - Write to JSON file (jsonlite)
   - `read_json_track()` - Read from JSON file (RcppSimdJson)
   - `read_track()` - Unified reader (SSFF or JSTF)
   - `as.data.frame.JsonTrackObj` - Convert to data.frame
   - `as_tibble.JsonTrackObj` - Convert to tibble
   - `append_json_track_slice()` - Add more slices
   - `merge_json_tracks()` - Combine multiple files
   - `subset_json_track()` - Filter by time range
   - `get_jstf_extension()` - Get extension for function

8. **Full specification**: See `JSON_TRACK_FORMAT_SPECIFICATION.md`

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
- **pladdrr functions**: Use `*p()` suffix for migrated functions (e.g., `trk_pitchp()`, `lst_voice_reportp()`)
  - Newer pladdrr functions use standard names (e.g., `trk_cpps()`, `lst_vq()`)
  - All use pladdrr R/C++ implementation (no Python)
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
- `R/ssff_*.R`: DSP function implementations (track-based, outputs time-series)
  - `R/ssff_c_assp_*.R`: ASSP C library wrappers (e.g., trk_forest, trk_ksvfo, trk_mhspitch)
  - `R/ssff_cpp_*.R`: Native C++ implementations (e.g., SPTK, ESTK wrappers)
    - `R/ssff_cpp_sptk_*.R`: SPTK-based functions (trk_rapt, trk_swipe, trk_dio, etc.)
    - `R/ssff_cpp_estk_*.R`: Edinburgh Speech Tools functions
  - `R/ssff_python_*.R`: Python-based implementations (e.g., trk_swiftf0, trk_crepe)
    - `R/ssff_python_pm_*.R`: Parselmouth/Praat-based functions (legacy naming, now use pladdrr)
  - `R/ssff_pladdrr_*.R`: pladdrr-based track functions (e.g., trk_cpps, trk_vuv)
- `R/list_*.R`: Summary statistic functions (outputs data frames/scalars)
  - `R/list_python_*.R`: Python-based summary functions (e.g., lst_vat, lst_voice_sauce)
    - `R/list_python_pm_*.R`: Parselmouth/Praat-based functions (legacy naming, now use pladdrr)
  - `R/list_pladdrr_*.R`: pladdrr-based summary functions (e.g., lst_vq, lst_pharyngeal)
  - `R/list_vat.R`: Voice Analysis Toolbox (132 dysphonia measures)
- `R/superassp_*.R`: Legacy DSP implementations (being migrated to ssff_* naming)
- `R/av_helpers.R`: Media loading and processing helpers (av_to_asspDataObj, processMediaFiles_LoadAndProcess)
- `R/pladdrr_helpers.R`: pladdrr-specific helpers (av_load_for_pladdrr, extract_sound_ptr, ptr_to_praatobj)
- `R/jstf_helpers.R`: JSON Track Format helpers (create_json_track_obj, write_json_track, read_json_track)
- `R/s7_avaudio.R`: S7 AVAudio class definition for in-memory audio
- `R/s7_methods.R`: S7 method registrations for automatic AVAudio dispatch
- `R/install_*.R`: Python module installation helpers (install_swiftf0, install_voice_analysis, etc.)
- `src/*.cpp`: C++ implementations and Rcpp bindings
  - `src/dsp_helpers.cpp`: Fast utility functions (fast_file_ext, fast_recycle_times, etc.)
  - `src/sptk_pitch.cpp`, `src/sptk_mfcc.cpp`: SPTK wrappers
  - `src/estk_pda.cpp`, `src/estk_pitchmark.cpp`: ESTK wrappers
- `src/assp/`: ASSP C library (legacy audio processing)
- `src/SPTK/`: SPTK submodule (pitch tracking, MFCC, spectral analysis)
- `src/ESTK/`: Edinburgh Speech Tools submodule (pitch detection, pitchmarking)
- `src/Makevars`: Compilation configuration (includes, flags, source lists)
- `inst/python/`: Python modules (voice_analysis_python, covarep_python, etc.)
- `inst/praat/`: Praat scripts called via Parselmouth (legacy, now use pladdrr)
- `tests/testthat/test-*.R`: Test files
- `inst/benchmarking/`: Benchmark scripts (R, Python, results, reports)

### SSFF Format
All track-based outputs use SSFF (Simple Signal File Format):
- Binary format for time-series data
- Written via `write.AsspDataObj()`
- Track names should be descriptive (e.g., "pitch[Hz]", "fm", "bw")

### DSP Function Extension Requirements

**CRITICAL**: All `trk_*` functions must implement file output capability with proper extension handling.

**Required Components**:

1. **Function Parameters**:
```r
trk_myfunction <- function(listOfFiles,
                           # ... DSP parameters ...
                           toFile = FALSE,            # REQUIRED
                           explicitExt = "ext",       # REQUIRED - default must match attr
                           outputDirectory = NULL,    # REQUIRED
                           verbose = TRUE) {
  # Implementation
}
```

2. **Function Attributes** (set at end of file):
```r
attr(trk_myfunction, "ext") <- "ext"  # MUST match explicitExt default
attr(trk_myfunction, "tracks") <- c("track1", "track2")
attr(trk_myfunction, "outputType") <- "SSFF"
attr(trk_myfunction, "nativeFiletypes") <- c("wav")
```

3. **File Writing Logic**:
```r
if (toFile) {
  base_name <- tools::file_path_sans_ext(basename(audio_path))
  out_dir <- if (is.null(outputDirectory)) dirname(audio_path) else outputDirectory
  output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))
  write.AsspDataObj(assp_obj, output_path)
  return(invisible(output_path))
}
```

**Extension Naming Conventions**:

Common extension patterns (40+ extensions in use):
- **Pitch tracking**: f0, sf0, yf0, dvf (pitch/F0 tracks)
- **Formants**: fms, pfm, dfm, dvfm (formant frequency tracks)
- **Spectral**: css, dft, lps, cep (spectral analysis)
- **Energy**: rms, zcr, acf, int (energy/amplitude measures)
- **Voice quality**: crk, vad, snr, c50 (voice quality indicators)
- **OpenSMILE**: ogs (eGeMAPS), emo (emobase), cmp (ComParE)

**Verification Checklist**:
- ✅ `explicitExt` parameter default matches `attr(*, "ext")`
- ✅ `toFile` parameter with default FALSE (backward compatibility)
- ✅ `outputDirectory` parameter for flexible output location
- ✅ File writing logic using `write.AsspDataObj()`
- ✅ Function attributes set at end of file
- ✅ Documentation includes `@param toFile`, `@param explicitExt`, `@param outputDirectory`
- ✅ `@return` documents different behavior when `toFile=TRUE` vs `FALSE`

**Reference Documentation**:
- See `DSP_EXTENSION_AUDIT_REPORT.md` for complete extension catalog
- See `FIXES_IMPLEMENTATION_SUMMARY.md` for implementation examples
- All functions audited as of 2025-10-29 - 100% compliance achieved

## Performance Considerations

- Prefer C++ implementations over Python (2-3x faster)
- Use `processMediaFiles_LoadAndProcess()` for automatic parallelization
- Avoid file I/O when possible - process in memory
- SPTK wrappers (`rapt`, `swipe`, etc.) recommended for users
- Low-level `_cpp` functions for advanced use cases only

## Dependencies

- **Required R packages** (Imports):
  - av, Rcpp, S7, parallel, cli, rlang - Core functionality
  - reticulate - Python integration (only when using Python-based functions)
  - tidyr, assertthat, readr, stringr, tools, digest, logger, uuid, R.matlab, dplyr, purrr - Utilities
- **Required System**:
  - C++11 compiler (gcc, clang, or MSVC)
  - No other system dependencies required
- **Optional Python modules** (only for Python-based functions):
  - `swift-f0`: Deep learning pitch tracker (install via `install_swiftf0()`)
  - `brouhaha`: VAD + SNR + C50 estimation (install via `install_brouhaha()`)
  - `deepformants`: Deep learning formant tracking (install via `install_deepformants()`)
  - `sacc`: SAcC pitch tracker (install via `install_sacc()`)
  - `dysprosody`: Prosodic assessment (install via `install_dysprosody()`)
  - `voice_analysis_python`: 132 dysphonia measures (install via `install_voice_analysis()`)
  - `pysptk`, `parselmouth`: Alternative pitch/formant implementations
  - **Note**: OpenSMILE no longer requires Python (native C++ implementation available)
- **Bundled Libraries** (compiled automatically):
  - **ASSP C library**: `src/assp/` - Core signal processing (formants, pitch, spectral analysis)
  - **SPTK C++ library**: `src/SPTK/` - Git submodule, pitch tracking (RAPT, SWIPE, REAPER, etc.)
  - **ESTK C++ library**: `src/ESTK/` - Git submodule, Edinburgh Speech Tools (pitch detection, pitchmarking)
  - **OpenSMILE C++ library**: `src/opensmile/` - Feature extraction (GeMAPS, eGeMAPS, ComParE, emobase)
  - All configured in `src/Makevars`

## Function Modernization Status (v0.6.0+)

### Overview

As of v0.6.0, superassp follows a unified architecture where **all DSP functions**:
1. ✅ Accept any media format via `av` package (WAV, MP3, MP4, video, etc.)
2. ✅ Process in memory using `av_to_asspDataObj()` or `processMediaFiles_LoadAndProcess()`
3. ✅ Support AVAudio S7 class with automatic dispatch (via `.setup_s7_methods()`)
4. ✅ Follow `trk_*` (tracks) or `lst_*` (summaries) naming conventions

**Current Status: 54% modernized** (29 of 54 functions fully compliant)

### Compliant Functions (29)

**C++ SPTK/ESTK Functions (8):** ✅ All modern
- trk_rapt, trk_swipe, trk_dio, trk_harvest, trk_reaper, trk_mfcc, trk_d4c, trk_estk_pitchmark

**C ASSP Functions (11):** ✅ All modern
- trk_forest, trk_mhspitch, trk_ksvfo, trk_acfana, trk_zcrana, trk_rmsana, trk_cepstrum, trk_lp_analysis, trk_cssSpectrum, trk_dftSpectrum, trk_lpsSpectrum

**Python Functions (10):** ✅ Already using av
- trk_swiftf0 (uses av::read_audio_bin)
- lst_vat, lst_voice_sauce (use av_load_for_python helper)
- 4× OpenSmile functions (GeMAPS, eGeMAPS, emobase, ComParE)
- 3× COVAREP functions (iaif, srh, vq)

### Functions Needing Migration (22)

**Python Functions Using librosa.load (11):** ⚠️ Need av migration
- trk_pyin, trk_yin, trk_crepe, trk_yaapt, trk_kaldi_pitch (High priority)
- trk_snackp, trk_snackf, trk_seenc, trk_excite, trk_aperiodicities, reaper_pm (Medium priority)

**Parselmouth Functions (10):** ⚠️ Need av integration
- ssff_python_pm_*.R (6 track functions)
- list_python_pm_*.R (4 summary functions)

**PyTorch Functions (2):** ⚠️ Need av loading
- trk_torch_pitch, trk_torch_mfcc

**Migration Guide:** See `MIGRATION_LIBROSA_TO_AV.md` for step-by-step instructions and reference implementations.

### Deleted Functions (v0.6.1)

The following Python implementations were superseded by faster C++ versions and removed:
- ~~nonopt_rapt~~ → use trk_rapt() (C++ SPTK)
- ~~nonopt_swipe~~ → use trk_swipe() (C++ SPTK)
- ~~dio_python~~ → use trk_dio() (C++ WORLD)
- ~~harvest_python~~ → use trk_harvest() (C++ WORLD)
- ~~Python REAPER~~ → use trk_reaper() (C++ SPTK)

**Benefit:** C++ implementations are 2-3x faster and require no Python dependencies.

## Key Recent Additions (v0.6.0+)

### OpenSMILE C++ Integration (v0.8.0)
- **Direct C++ integration** replacing Python bindings (5.5x faster)
- All OpenSMILE functions now default to C++ mode with Python fallback
- GeMAPS, eGeMAPS, ComParE, emobase all available via native C++
- Zero Python dependency for OpenSMILE features when using `use_cpp = TRUE` (default)
- Located: `src/opensmile_wrapper.cpp`, `inst/opensmile/`
- See v0.8.0 section below for detailed performance benchmarks

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
- Dysprosody: `lst_dysprosody()` - 193 prosodic features (v0.7.1+)
- Voxit: `lst_voxit()` - 11 voice/articulation complexity measures (v0.8.8+)
- Each has `install_*()`, `*_available()`, `*_info()` helpers

### Brouhaha-VAD (v0.8.0+)
- **Voice Activity Detection** with deep learning (50-100x faster than original)
- **SNR Estimation**: Signal-to-Noise Ratio tracking
- **C50 Estimation**: Room clarity/reverberation measure
- Installation: `install_brouhaha()` (R/install_brouhaha.R)
- Function: `trk_brouhaha()` (R/ssff_python_brouhaha.R)
- Location: `inst/python/brouhaha-vad/`
- **Key Features**:
  - Joint VAD + SNR + C50 prediction from single model
  - 50-100x performance improvement through Cython/Numba optimizations
  - 100% faithful to original results (verified)
  - GPU support with automatic device detection
  - Parallel batch processing
  - Custom model support
- **Performance Tiers**:
  - Basic (3-10x): Python vectorization (always active)
  - High (10-30x): + Numba JIT (`install_brouhaha(install_numba=TRUE)`)
  - Maximum (50-100x): + Cython compilation (`install_brouhaha(compile_cython=TRUE)`)
- **Output**: AsspDataObj with 3 tracks: `vad` (binary), `snr` (dB), `c50` (dB)
- **Integration**: Full superassp compliance, emuR compatible, supports any media format
- **Documentation**: See `inst/python/brouhaha-vad/README.md` and `BROUHAHA_INTEGRATION_SUMMARY.md`

### Voxit (v0.8.8+)
- **Voice and Articulation Complexity** measures for prosodic analysis
- **11 features**: Speaking rate, pause statistics, rhythmic complexity, pitch dynamics
- Installation: `install_voxit()` (R/install_voxit.R)
- Function: `lst_voxit()` (R/list_voxit.R)
- Location: `inst/python/voxit/`
- **Key Features**:
  - Temporal features: WPM, pause counts/durations, rhythmic complexity (Lempel-Ziv)
  - Pitch features: Range, entropy, velocity, acceleration (with Savitzky-Golay smoothing)
  - Requires word alignments (CSV format)
  - Uses SAcC for pitch tracking (`install_sacc()`)
  - 2-3x speedup with Numba JIT optimization
  - 3-5x speedup with Cython compilation
- **Performance Tiers**:
  - Standard: Pure Python with NumPy (~200ms/5s audio)
  - Numba: JIT compilation (~80ms/5s audio)
  - Cython: Compiled C extensions (~60ms/5s audio)
- **Output**: Named list with 11 prosodic/rhythmic features
- **Integration**: Full superassp compliance, av package audio loading, parallel processing
- **Documentation**: See `inst/python/voxit/README.md` and `VOXIT_INTEGRATION_SUMMARY.md`


### Legacy STRAIGHT Vocoder (v0.9.0+)
- **High-quality vocoder** for speech analysis and synthesis
- **Components**: F0 extraction, spectral analysis, aperiodicity, synthesis
- Installation: `install_legacy_straight()` (R/install_legacy_straight.R)
- Functions: `trk_straight_f0()`, `trk_straight_spec()` (R/ssff_python_straight_f0.R, R/ssff_python_straight_spec.R)
- Location: `inst/python/legacy_STRAIGHT/`
- **Key Features**:
  - Multi-cue F0 extraction (~91% frame accuracy, ~96.5% mean F0 accuracy)
  - Pitch-adaptive spectral analysis (99.996% MATLAB correlation)
  - High-quality synthesis (99.99% MATLAB correlation)
  - Numba JIT optimization (~20% speedup)
  - No external dependencies beyond NumPy/SciPy
- **Performance**:
  - F0 extraction: ~0.68s for 0.79s audio (0.86x RT with Numba)
  - Spectral analysis: ~0.15s
  - Synthesis: ~0.05s
- **Accuracy vs MATLAB**:
  - F0: ~91% frame accuracy (known limitation: octave errors < 100 Hz)
  - Spectral: 99.996% correlation
  - Synthesis: 99.99% correlation
- **Status**: Integrated, documented, tested. Segfault issue being resolved.
- **Documentation**: See `inst/python/legacy_STRAIGHT/README.md`, `STRAIGHT_95_PERCENT_PLAN.md`, `STRAIGHT_INTEGRATION_STATUS_NOV1_2025.md`
- **Roadmap**: Systematic improvement to 95% accuracy planned (2-3 weeks)

### Pladdrr Integration (v0.11.2) - COMPLETE ✅

**Achievement**: Complete migration from Python's parselmouth to R's pladdrr

- **14 core functions** migrated or created (10 track, 4 summary)
- **Pure R/C++ implementation** - No Python for Praat-based functions
- **100% plabench coverage** - All 16 reference implementations ported
- **Performance**: 2-15x faster than parselmouth equivalents
- **Located**: `R/ssff_pladdrr_*.R`, `R/list_pladdrr_*.R`

**Migration Batches**:
- **Batch 1** (Sessions 3-4): `trk_intensityp()`, `trk_pitchp()`, `trk_formantp()`
- **Batch 2** (Session 5): `lst_voice_reportp()`, `lst_dsip()`, `lst_voice_tremorp()`, `lst_avqip()`
- **Batch 3** (Session 6): `trk_spectral_momentsp()`, `trk_praatsaucep()` (36 measures!)
- **Phase 4** (Session 7): `trk_cpps()`, `trk_vuv()`, `lst_vq()`, `lst_pharyngeal()` (68 measures!)

**New Functions (Phase 4)**:
- `trk_cpps()` - Cepstral Peak Prominence Smoothed (voice quality)
- `trk_vuv()` - Voice/Unvoiced Detection (dual output: TextGrid/SSFF)
- `lst_vq()` - Voice quality summary (36 measures)
- `lst_pharyngeal()` - Pharyngeal voice quality (68 measures - most comprehensive!)

**Integrated Functions**:
- `trk_formantpathp()` - MERGED into `trk_formantp()` (HMM tracking)
- MOMEL - INTEGRATED in `lst_dysprosody()`
- INTSINT - INTEGRATED in `lst_dysprosody()`

**Helper Infrastructure**:
- `pladdrr_helpers.R` - Audio loading, pointer extraction, format conversion
- `jstf_helpers.R` - JSON Track Format I/O for lst_* functions
- `av_load_for_pladdrr()` - Flexible audio loading with resampling
- `extract_sound_ptr()` - Extract C pointers for direct API access
- `ptr_to_praatobj()` - Convert pointers to R6 PraatObj

**pladdrr API Patterns**:
- **Direct API (Tier 2)**: `to_pitch_cc_direct()`, `to_formant_direct()`, etc.
- **Ultra API (Tier 4)**: Batch operations for 5-10x speedups
- **R6 Methods**: `sound$extract_part()`, `spectrum$to_ltas_1to1()`, etc.
- **Internal API**: Namespace access for specialized functions

**Key Features**:
- JSTF output for all lst_* functions (JSON Track Format)
- Dual output format support (trk_vuv: TextGrid + SSFF)
- Two-pass adaptive pitch (speaker-specific F0 range)
- Ultra API batch operations (jitter/shimmer, multi-band HNR)
- Comprehensive voice quality (68 pharyngeal measures!)

**Performance Wins**:
- `lst_vq()`: 5-10x faster jitter/shimmer (Ultra API batch ops)
- `lst_vq()`: 2-2.5x faster multi-band HNR
- `lst_pharyngeal()`: 15.7x faster vs pladdrr v4.8.14
- Overall: 2-15x faster than parselmouth

**Requirements**:
- pladdrr >= 4.8.16 (formant bug fixes)
- Formant+intensity integration reported fixed (testing pending)
- Formant window extraction bug reportedly fixed

**Known Issues (Testing Pending)**:
- `trk_formantp()` intensity disabled (segfault workaround, reported fixed)
- `lst_pharyngeal()` uses full-sound formants (window bug workaround, reported fixed)
- Both workarounds can be removed after testing with latest pladdrr

**Documentation**:
- See `PLADDRR_MIGRATION_STATUS.md` for complete status
- See `PLADDRR_FINAL_STATUS.md` for project analysis
- See `SESSION_7_SUMMARY.md` for Phase 4 details
- See `NEWS.md` v0.11.2 for release notes

**Timeline**:
- Started: 2026-02-03 (Session 3)
- Completed: 2026-02-06 (Session 7)
- Duration: 4 days (7 sessions)
- **20 days ahead of schedule!** 🚀

## Package Version History

### v0.8.7 (Current - In Development)
- **Parselmouth In-Memory Processing**: All parselmouth functions now use in-memory processing
  - Added `av_load_for_parselmouth()` and `av_to_parselmouth_sound()` helpers
  - Updated `lst_dysprosody()` to eliminate temporary files (38% faster)
  - Updated `trk_pitchp()` for in-memory Sound object processing
  - Pattern established for remaining functions (5 pending migration)
  - Full documentation in AV_TO_PARSELMOUTH_STRATEGY.md
  - Comprehensive audit in PARSELMOUTH_FUNCTIONS_AUDIT.md
  - Migration plan in PYTHON_INMEMORY_MIGRATION_PLAN.md
- **Architecture**: Two patterns for parselmouth integration documented
  - Pattern 1: Python creates Sound from numpy (6 functions already compliant)
  - Pattern 2: R creates Sound directly (2 functions migrated, 5 pending)

### v0.8.6
- **Brouhaha-VAD Integration**: Voice Activity Detection + SNR + C50 estimation (v0.8.0-0.8.3)
  - 50-100x performance improvement through optimizations
  - Cython compilation support for maximum speed
  - Numba JIT for instant speedup without compilation
  - GPU support with automatic device detection
  - Parallel batch processing
  - Full superassp interface compliance
  - emuR database integration
- New functions: `trk_brouhaha()`, `install_brouhaha()`, `brouhaha_available()`, `brouhaha_info()`
- Complete Python module in `inst/python/brouhaha-vad/`
- Comprehensive documentation and integration guide

### v0.7.1
- Added dysprosody prosodic assessment module (193 features)
- MOMEL-INTSINT pitch target extraction
- Spectral tilt with Iseli-Alwan harmonic correction
- Full av package integration for universal media support
- Performance: ~0.16-0.44s per file (14x realtime)

### v0.7.0
- **Complete librosa migration** - All functions now use av package
- **Universal media format support** - WAV, MP3, MP4, video files, etc.
- **Deleted redundant PyTorch functions** (use faster C++ alternatives)
- **6 functions migrated**: trk_pyin, trk_yin, trk_crepe, trk_yaapt, trk_seenc, trk_excite
- **Breaking changes**: Removed trk_kaldi_pitch, trk_torch_pitch, trk_torch_mfcc

### v0.6.0
- Introduced S7 AVAudio class for in-memory processing
- Swift-F0 deep learning pitch tracker integration
- Automatic format fallback (av → wrassp for niche formats)
- Modernized 29 of 54 functions to unified architecture

### Earlier Versions
See git history and NEWS.md for complete version history.

## Key References

### Documentation Files
- **MIGRATION_LIBROSA_TO_AV.md**: Guide for migrating Python functions from librosa to av
- **MIGRATION_EXAMPLE.md**: Step-by-step migration examples
- **README.md**: User-facing documentation with benchmarks
- **NEWS.md**: Comprehensive version history and changelog
- **INDEX.md**: Function inventory and quick reference

### Technical Notes
- **EMUR_COMPATIBILITY_ANALYSIS.md**: Integration with emuR framework
- **OPTIMIZATION_PROPOSAL_MEMORY_DSP.md**: Memory optimization strategies
- **PARSELMOUTH_MEMORY_OPTIMIZATION.md**: Praat/Parselmouth optimizations
- **UNIFORM_PLACEHOLDER_STRATEGY.md**: Handling missing values in DSP output

### Parselmouth In-Memory Processing Documentation
- **AV_TO_PARSELMOUTH_STRATEGY.md**: Complete strategy for av → parselmouth Sound conversion
- **PARSELMOUTH_FUNCTIONS_AUDIT.md**: Comprehensive audit of all parselmouth functions
- **PYTHON_INMEMORY_MIGRATION_PLAN.md**: Systematic migration plan for remaining functions
- **DYSPROSODY_INMEMORY_IMPLEMENTATION.md**: Detailed implementation notes for lst_dysprosody

### Function Organization and Discovery
- **PKGDOWN_FUNCTION_GROUPING.md**: Comprehensive catalog of all 75+ DSP functions organized by use case
  - 9 main categories: Pitch/F0, Formants, Spectral, Energy, Voice Quality, Prosody, Source-Filter, OpenSMILE, Acoustic Features
  - Performance tiers (C++ > C > Python DL > Python Classical > Parselmouth)
  - Usage recommendations for different scenarios
  - Ready-to-use pkgdown reference YAML configuration
- **FUNCTION_PARAMETERS_REFERENCE.md**: Comprehensive catalog of all 97 unique parameters used across functions
  - Parameter descriptions, default values, and types
  - Functions using each parameter
  - Parameter usage patterns (universal, pitch-specific, spectral-specific, etc.)
  - Development guidelines for parameter standardization

### DSP Extension Audit and Compliance
- **DSP_EXTENSION_AUDIT_REPORT.md**: Complete audit of all 61+ DSP functions (2025-10-29)
  - Extension catalog with 40+ file extensions
  - Function inventory by implementation type (C ASSP, C++ SPTK, Python, etc.)
  - Mismatch analysis and missing attribute identification
  - Testing recommendations and compliance verification
- **DSP_EXTENSION_FIXES_COMPLETED.md**: Detailed implementation log of all fixes
  - 6 functions fixed (lst_eGeMAPS, lst_emobase, trk_creak_union, trk_formants_tvwlp, trk_dv_f0, trk_dv_formants)
  - Step-by-step implementation details with code examples
  - Usage examples and testing recommendations
- **FIXES_IMPLEMENTATION_SUMMARY.md**: Executive summary of DSP extension fixes
  - Statistics: 6 issues fixed, 2 new extensions introduced (dvf, dvfm)
  - 100% backward compatibility maintained
  - Verification status and next steps

### Package Audit (2025-11-01)
- **PACKAGE_AUDIT_2025-11-01.md**: Comprehensive package interface and consistency audit
  - **Critical Fix**: NAMESPACE regenerated to export Phonet functions
  - Interface consistency analysis across 195+ exported functions
  - Deprecation status and recommendations
  - Function categorization by domain and implementation type
  - Grade: A- (excellent consistency, minor Python env parameter naming issue)

## Function Categorization Quick Reference

### By Prefix (Function Naming Conventions)

| Prefix | Count | Purpose | Examples |
|--------|-------|---------|----------|
| `trk_` | 50+ | SSFF track output for time-series | `trk_phonet`, `trk_yin`, `trk_crepe` |
| `lst_` | 15 | List/data.frame output for analysis | `lst_phonet`, `lst_dysprosody`, `lst_voxit` |
| `install_` | 15 | Python dependency installation | `install_phonet`, `install_brouhaha` |
| `*_available` | 13 | Check if tool is installed | `phonet_available()` |
| `*_info` | 12 | Get tool configuration info | `phonet_info()` |
| `*_cpp` | 12 | Low-level C++ functions | `yin_cpp`, `harvest_cpp` |

### By Domain (Speech Analysis Categories)

| Domain | Functions | Primary Implementation | Examples |
|--------|-----------|------------------------|----------|
| **Pitch** | 20 | C++ (7), Python (11), C (1) | `trk_yin`, `trk_crepe`, `trk_sacc` |
| **Formants** | 5 | Python (5) | `trk_formantp`, `trk_deepformants` |
| **Voice Quality** | 11 | Python (6), C++ (1), R (4) | `trk_brouhaha`, `trk_d4c`, `lst_voxit` |
| **Spectral** | 6 | C (3), C++ (1), Python (2) | `trk_mfcc`, `trk_cepstrum` |
| **Features** | 10 | C++ (2), Python (8) | `lst_GeMAPS`, `lst_dysprosody` |
| **Phonology** | 2 | R/Python (2) | `lst_phonet`, `trk_phonet` |
| **Prosody** | 1 | Python (1) | `lst_dysprosody` |

### By Implementation Type (Performance Guide)

| Type | Count | Speed | When to Use | Examples |
|------|-------|-------|-------------|----------|
| **C++** | ~25 | ⚡⚡⚡ Fastest | Production, batch processing | `trk_yin`, `trk_harvest`, `trk_mfcc` |
| **C (ASSP)** | ~10 | ⚡⚡ Fast | Legacy compatibility | `trk_forest`, `trk_cepstrum` |
| **Python** | ~30 | ⚡ Slower | Deep learning, Praat integration | `trk_crepe`, `trk_formantp`, `lst_phonet` |
| **R** | ~130 | ⚡ Variable | Wrappers, utilities, glue code | Various helpers |

### Interface Consistency Status

| Aspect | Status | Standard | Notes |
|--------|--------|----------|-------|
| File input | ✅ Consistent | `listOfFiles` | 70 functions |
| Time windowing | ✅ Consistent | `beginTime`/`endTime` | 70 functions |
| Output control | ✅ Consistent | `toFile` | 50 functions |
| Verbosity | ✅ Consistent | `verbose` | 73 functions |
| Function naming | ✅ Consistent | `trk_`, `lst_`, `install_` | All functions |
| Python env | ⚠️ Minor issue | `envname` vs `conda.env` | Low priority |

**For complete audit details, see PACKAGE_AUDIT_2025-11-01.md**

## Working with This Codebase

### First-Time Setup
```bash
# Clone with submodules
git clone --recursive https://github.com/humlab-speech/superassp.git

# Or if already cloned
cd superassp
git submodule update --init --recursive

# Install in R
Rscript -e "devtools::install()"
```

### Current Development Branch

**Branch**: `cpp_optimization`
**Focus**: C++ performance optimizations and native implementations

When working on this branch:
- Prioritize C++ implementations over Python where possible
- Maintain backwards compatibility with existing APIs
- Run full benchmark suite after major changes
- Test compilation across platforms (macOS, Linux, Windows)

### Development Cycle
```r
# 1. Load package for development
devtools::load_all()

# 2. Make changes to R or C++ code

# 3. If C++ changed, update exports
Rcpp::compileAttributes()

# 4. Regenerate documentation
devtools::document()

# 5. Run tests
devtools::test()

# 6. Check package
devtools::check()

# 7. Build (if needed)
devtools::build()
```

### Before Committing
```r
# Always run these before committing
devtools::document()  # Regenerate docs
devtools::test()      # Run tests
devtools::check()     # Full package check
```

### Release Workflow
1. Update NEWS.md with changes
2. Bump version in DESCRIPTION
3. Regenerate documentation: `devtools::document()`
4. Run full check: `devtools::check()`
5. Build package: `devtools::build()`
6. Create git tag: `git tag v0.X.Y`
7. Push with tags: `git push --tags`

## Critical Files - Do Not Modify Directly

**Auto-generated files (regenerate, don't edit)**:
- `R/RcppExports.R` - Generated by `Rcpp::compileAttributes()`
- `src/RcppExports.cpp` - Generated by `Rcpp::compileAttributes()`
- `NAMESPACE` - Generated by `devtools::document()`
- `man/*.Rd` - Generated by `devtools::document()` from roxygen2 comments

**Submodule code (modify upstream only)**:
- `src/SPTK/` - SPTK submodule
- `src/ESTK/` - Edinburgh Speech Tools submodule
- `src/tcl-snack/` - Snack submodule
- `inst/onnx/swift-f0/` - Swift-F0 submodule
- `inst/python/DeepFormants/` - DeepFormants submodule

**Configuration files (edit carefully)**:
- `src/Makevars` - Compilation configuration (critical for builds)
- `.gitmodules` - Git submodule configuration
- `DESCRIPTION` - Package metadata (note: Depends on wrassp, not Imports)

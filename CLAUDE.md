# CLAUDE.md

**superassp**: R package for speech signal processing. Self-contained (no wrassp required). Unified interface via SSFF/AsspDataObj. Outputs compatible with emuR.

## Quick Reference

```r
# DEV WORKFLOW
devtools::load_all()              # Load for development
Rcpp::compileAttributes()         # After C++ changes
devtools::document()              # Regen docs
devtools::test()                  # Run tests
devtools::check()                 # Full check

# BEFORE COMMITTING
devtools::document()
devtools::test()

# SUBMODULES (first time)
git submodule update --init --recursive
git submodule update --remote       # Update to latest
```

## Adding DSP Functions

### C++ (preferred for speed)
1. Implement `src/myfunction.cpp` with `// [[Rcpp::export]]`
2. `Rcpp::compileAttributes()`
3. Create R wrapper `R/ssff_cpp_myfunction.R`
4. `devtools::document()`
5. Add tests `tests/testthat/test-myfunction.R`

### Python (specialized algorithms)
1. Create `inst/python/mymodule/` or `inst/python/myscript.py`
2. Create `R/install_mymodule.R` (install, check, info functions)
3. Create R wrapper `R/ssff_python_myfunction.R`
4. **Use `av::read_audio_bin()` for audio (NOT librosa)**
5. `devtools::document()` + tests

### ASSP Library (C)
Use unified `processMediaFiles_LoadAndProcess()` helper:
```r
result <- processMediaFiles_LoadAndProcess(
  listOfFiles, beginTime, endTime,
  nativeFiletypes = c("wav", "au"),
  fname = "assp_function_name",
  toFile, verbose, param1 = ...,
  explicitExt, outputDirectory
)
```

### Praat Functions
1. Create script `inst/praat/` (accepts file path, outputs CSV)
2. Create R wrapper via reticulate/Parselmouth
3. Read CSV, convert to AsspDataObj
4. Keep `praat_` prefix

### Summary Functions (`lst_*`)
If writing JSTF (JSON Track Format) files:
```r
lst_function <- function(listOfFiles, ..., toFile = FALSE,
                         explicitExt = "ext", outputDirectory = NULL) {
  results <- your_processing()
  if (toFile) {
    json_obj <- create_json_track_obj(results, ...)
    write_json_track(json_obj, output_path)
    return(invisible(output_path))
  }
  return(results)
}
attr(lst_function, "ext") <- "ext"
attr(lst_function, "outputType") <- "JSTF"
```

## Function Naming & Requirements

**Prefixes**:
- `trk_*` — Time-series tracks (signal following audio, e.g., F0, formants, energy)
- `lst_*` — Summary statistics (scalars/vectors, e.g., jitter, voice quality measures)
- `install_*`, `*_available`, `*_info` — Python module helpers

**Required for all `trk_*` functions**:
- Parameter: `toFile` (default FALSE), `explicitExt`, `outputDirectory`
- Attributes: `ext`, `tracks`, `outputType` ("SSFF"), `nativeFiletypes`
- Logic: Write SSFF if `toFile=TRUE`, return AsspDataObj if `FALSE`

**Standard parameters** (if applicable):
- `listOfFiles` — Input path(s)
- `beginTime`, `endTime` — Time windowing (seconds)
- `windowShift` — Frame shift (ms)
- `minF`, `maxF` — Pitch range
- `verbose` — Progress messages

## File Organization

- `R/ssff_*.R` — Track-based DSP functions
  - `R/ssff_c_assp_*.R` — ASSP C library wrappers
  - `R/ssff_cpp_*.R` — C++ implementations
  - `R/ssff_python_*.R` — Python-based
  - `R/ssff_pladdrr_*.R` — pladdrr (Praat from R)
- `R/list_*.R` — Summary functions
  - `R/list_python_*.R` — Python summaries
  - `R/list_pladdrr_*.R` — pladdrr summaries
- `R/av_helpers.R` — Media loading (`av_to_asspDataObj`, `processMediaFiles_LoadAndProcess`)
- `R/pladdrr_helpers.R` — Praat integration (`av_load_for_pladdrr`, pointer extraction)
- `R/jstf_helpers.R` — JSON Track Format I/O
- `R/s7_avaudio.R`, `R/s7_methods.R` — In-memory audio (S7 class + dispatch)
- `R/install_*.R` — Python module installation
- `src/*.cpp` — C++ implementations
- `src/Makevars` — Compilation config (CRITICAL—maintains include paths, source lists)
- `src/assp/`, `src/SPTK/`, `src/ESTK/` — Bundled libraries
- `src/opensmile/` — OpenSMILE (bundled, editable)
- `inst/python/` — Python modules
- `inst/praat/` — Praat scripts (legacy)
- `tests/testthat/test-*.R` — Tests

## Architecture: Three-Layer DSP

**Layer 1: Core implementations**
- C/C++: ASSP (`src/assp/`), ESTK (`src/ESTK/`), SPTK (`src/SPTK/`)
- Praat scripts (`inst/praat/`) via Parselmouth
- Python modules

**Layer 2: Low-level (`_cpp` functions)**
- Direct C++ bindings: `rapt_cpp()`, `swipe_cpp()`, `reaper_cpp()`
- Require AsspDataObj input
- Exposed in `R/RcppExports.R` (auto-generated)
- Fast but less user-friendly

**Layer 3: High-level wrappers** (recommended)
- User-facing: `rapt()`, `swipe()`, `reaper()`, etc.
- Handle any media format via av package
- Auto-parallelize for 2+ files
- Support `toFile=TRUE/FALSE`

## Media Processing

**Modern pattern (preferred)**:
1. `av_to_asspDataObj()` — Load any format to memory
2. `processMediaFiles_LoadAndProcess()` — Unified DSP + parallelization
3. `.External("performAsspMemory", ...)` — In-memory via C interface

**Key helpers** (`R/av_helpers.R`):
- `av_to_asspDataObj()` — Convert any format (auto-fallback: av → wrassp for niche formats)
- `processMediaFiles_LoadAndProcess()` — Batch + parallel processor

**Parselmouth pattern**:
- Pattern 1: Python creates Sound from numpy (6 functions using this)
- Pattern 2: R creates Sound directly (2 functions migrated)
- Helper: `av_load_for_parselmouth(file, start, end, channels, sample_rate)`

**S7 AVAudio class** (in-memory dispatch):
- S7 generics + method registration
- Automatic dispatch for `trk_*`/`lst_*` functions
- Zero file I/O for preprocessing

## Dependencies

**Required** (Imports):
- av, Rcpp, S7, parallel, cli, rlang
- reticulate (only if using Python-based functions)
- tidyr, assertthat, readr, stringr, tools, digest, logger, uuid, R.matlab, dplyr, purrr

**System**: C++11 compiler (gcc, clang, MSVC)

**Optional Python modules** (install via `install_*()` helpers):
- swift-f0, brouhaha, deepformants, sacc, dysprosody, voice_analysis_python, pysptk, parselmouth

**Bundled**:
- ASSP C library (`src/assp/`)
- SPTK (`src/SPTK/` submodule)
- ESTK (`src/ESTK/` submodule)
- OpenSMILE C++ (`src/opensmile/`, editable)

## Testing Patterns

```r
test_that("function works with single file", {
  skip_if_not_installed("superassp")
  test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
  skip_if(test_wav == "", "Test file not found")
  result <- my_function(test_wav, toFile = FALSE, verbose = FALSE)
  expect_s3_class(result, "AsspDataObj")
  expect_true("track_name" %in% names(result))
})
```

## Troubleshooting

**C++ compilation**: `undefined reference to SPTK::...`
```bash
git submodule update --init --recursive
devtools::clean_dll()
devtools::load_all()
```

**Function not exported**: Add `@export` to roxygen2 comment, then `devtools::document()`

**Python module not found**: `reticulate::py_config()` to check setup. Run `install_yourmodule()`

**av::read_audio_bin fails**: Check av package installed + FFmpeg available (`av::av_video_info()`)

**Tests timeout**: Disable parallel: `trk_rapt(..., parallel = FALSE)`

## Export Policy

User-exportable only:
- `trk_*`, `lst_*` — DSP functions
- `ucnv_*` — Unit conversion
- `read_*`, `write_*` — I/O
- Class definitions: `AVAudio`
- S7 generics: `dur`, `numRecs`, `rate`, `startTime`, `tracks`
- Type predicates: `is.AsspDataObj`

Internal (no `@export`):
- `_cpp` functions
- Helper functions

## Critical Files (DO NOT EDIT DIRECTLY)

Auto-generated (regenerate instead):
- `R/RcppExports.R` — Generated by `Rcpp::compileAttributes()`
- `src/RcppExports.cpp` — Generated by `Rcpp::compileAttributes()`
- `NAMESPACE` — Generated by `devtools::document()`
- `man/*.Rd` — Generated from roxygen2 comments

Submodule code (modify upstream only):
- `src/SPTK/`, `src/ESTK/`, `src/tcl-snack/`, `inst/onnx/swift-f0/`, `inst/python/DeepFormants/`

Configuration (edit carefully):
- `src/Makevars` — Compilation config (critical)
- `.gitmodules` — Submodule config
- `DESCRIPTION` — Package metadata

## Implementation Focus

**Prioritize C++ over Python** (2-3x faster). Use Python only for specialized algorithms (deep learning, Praat integration). **Always use `av::read_audio_bin()` for audio loading**, not librosa. Test with multiple media formats (WAV, MP3, MP4). **Regenerate docs before committing**: `devtools::document()` + `devtools::test()`.

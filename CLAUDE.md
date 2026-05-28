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
1. Implement `src/myfunction.cpp` with `// [[Rcpp::export]]`. Mark internal Rcpp bindings `@keywords internal` + `@noRd` — never `@export`.
2. `Rcpp::compileAttributes()`
3. Create R wrapper `R/ssff_cpp_<origin>_<algo>.R` (origin = `sptk`, `estk`, `snack`, `covarep`, …)
4. Wrapper obtains audio via `assp_load_audio_for_dsp(file, begin, end, framework = "sptk")` (or appropriate framework)
5. `devtools::document()`; add tests `tests/testthat/test-<algo>.R`

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

### pladdrr (Praat from R)
1. Create R wrapper `R/ssff_pladdrr_<algo>.R` (or `R/list_pladdrr_<algo>.R` for summary outputs)
2. Wrapper obtains a pladdrr `Sound` via `assp_load_audio_for_dsp(file, begin, end, framework = "pladdrr")` (which delegates to `av_load_for_pladdrr()` and handles av-fallback transcoding)
3. Convert pladdrr output to AsspDataObj for `trk_*` or to a list for `lst_*`

### Summary Functions (`lst_*`)
If writing JSTF (JSON Track Format) files:
```r
lst_function <- function(listOfFiles, ..., toFile = FALSE,
                         explicitExt = "ext", outputDirectory = NULL) {
  results <- your_processing()
  if (toFile) {
    json_obj <- create_json_track_obj(results, ...)
    write_jstf(json_obj, output_path)
    return(invisible(output_path))
  }
  return(results)
}
attr(lst_function, "ext") <- "ext"
attr(lst_function, "outputType") <- "JSTF"
```

## Function Naming & Requirements

**Allowed user-exported prefixes**:
- `trk_*` — Time-series tracks (signal following audio, e.g., F0, formants, energy)
- `lst_*` — Summary statistics (scalars/vectors, e.g., jitter, voice quality measures)
- `ucnv_*` — Unit conversions (Hz↔Bark/Mel/ERB/semitone, dB↔phon/sone)
- `read_*` / `write_*` — I/O for SSFF, JSTF, audio

Within `trk_*`, names follow `trk_<metric>_<algorithm>` where the metric is non-obvious from the algorithm name. Examples:
- `trk_pitch_rapt`, `trk_pitch_swipe`, `trk_pitch_ac`, `trk_pitch_cc`, … — pitch trackers
- `trk_pitchmark_estk`, `trk_pitchmark_reaper` — pitch-mark detectors
- `trk_formant_burg`, `trk_formant_forest`, `trk_formant_snack`, `trk_formant_cgdzp`, `trk_formant_tvwlp` — formant trackers
- `trk_rms`, `trk_acf`, `trk_zcr`, `trk_lpc`, … — single-metric analyses (no `_ana` suffix)
- `trk_dft_spectrum`, `trk_lps_spectrum`, `trk_css_spectrum` — spectra (snake_case, never camelCase)

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

R-source filenames follow `<output_kind>_<implementation_origin>_<algorithm>.R`.

- **output_kind**: `ssff` (time-series track), `list` (summary), `read`, `write`, `ucnv`, `helpers`
- **implementation_origin**: `c_assp`, `cpp`, `cpp_sptk`, `cpp_estk`, `cpp_snack`, `cpp_covarep`, `cpp_opensmile`, `pladdrr`, `r` (pure R)

Examples:
- `R/ssff_c_assp_rmsana.R` — libassp C, RMS analyser
- `R/ssff_cpp_sptk_rapt.R` — SPTK C++, RAPT pitch
- `R/ssff_cpp_covarep_iaif.R` — COVAREP C++, IAIF
- `R/ssff_cpp_estk_pda.R` — ESTk C++, PDA
- `R/ssff_pladdrr_pitch.R` — pladdrr (Praat) pitch
- `R/list_pladdrr_avqi.R` — pladdrr AVQI summary
- `R/list_cpp_opensmile_eGeMAPS.R` — openSMILE C++ feature set
- `R/list_r_polarity.R`, `R/list_r_voxit.R`, `R/list_r_vowel_space.R` — pure R summaries
- `R/helpers_dysprosody_*.R` — internal helpers for the dysprosody pipeline

Class / format infrastructure:
- `R/assp_dataobj*.R`, `R/assp_generics.R`, `R/assp_dataobj_methods.R` — AsspDataObj class + S3 generics
- `R/jstf_*.R` — JSTF format infrastructure (`read_jstf`, `write_jstf`, validators)
- `R/s7_avaudio.R`, `R/s7_methods.R` — internal S7 AVAudio dispatch
- `R/read_audio.R`, `R/read_ssff.R`, `R/write_ssff.R` — public I/O
- `R/audio_loader.R` — `assp_load_audio_for_dsp()` uniform helper for DSP wrappers
- `R/av_helpers.R`, `R/pladdrr_helpers.R`, `R/sptk_helpers.R`, `R/wav_helpers.R`, `R/prep_recode.R` — internal media plumbing

Bundled C/C++ libraries (do not modify directly):
- `src/assp/`, `src/SPTK/` (submodule), `src/ESTK/` (submodule), `src/tcl-snack/` (submodule)
- `src/opensmile/` — OpenSMILE C++ (bundled, editable)
- `src/Makevars` — Compilation config (CRITICAL—include paths and source lists)
- `tests/testthat/test-*.R` — Tests

## Architecture: Three-Layer DSP

**Layer 1: Core implementations**
- C/C++: ASSP (`src/assp/`), ESTk (`src/ESTK/`), SPTK (`src/SPTK/`), OpenSMILE (`src/opensmile/`), COVAREP (bundled C++)
- pladdrr (Praat from R) for Praat-native algorithms

**Layer 2: Low-level Rcpp bindings (`*_cpp` functions)**
- Direct C++ bindings: `rapt_cpp()`, `swipe_cpp()`, `reaper_cpp()`, `estk_pitchmark_cpp()`, …
- Require AsspDataObj input
- Exposed in `R/RcppExports.R` (auto-generated)
- **Always internal** — never `@export`. Mark with `@keywords internal` + `@noRd` in the C++ roxygen header.

**Layer 3: High-level wrappers** (user-facing)
- `trk_pitch_rapt()`, `trk_pitch_swipe()`, `trk_pitch_reaper()`, `trk_rms()`, … etc.
- Obtain audio via `assp_load_audio_for_dsp()` (uniform fallback contract)
- Auto-parallelize for 2+ files
- Support `toFile=TRUE/FALSE`

## Media Processing

**Audio-loading contract** (every DSP wrapper):
1. Try the framework's native loader (libassp via `processMediaFiles_LoadAndProcess()`, pladdrr via `av_load_for_pladdrr()`, etc.)
2. On failure / unsupported format, fall back to `read_audio()` (which itself does libassp-then-`av`)
3. If the framework can't consume an AsspDataObj, dump to a temp WAV via `avaudio_to_tempfile()` and pass the path

`assp_load_audio_for_dsp(file, begin, end, samples, framework)` (in `R/audio_loader.R`) wraps this contract — call it from new wrappers.

**Public reader**: `read_audio(fname, begin = 0, end = 0, samples = FALSE)` — single user-facing audio loader. Mirrors `wrassp::read.AsspDataObj`. Sample-accurate windowing for variable-rate containers via `av::av_media_info()`.

**Key internal helpers** (all `@keywords internal`):
- `R/av_helpers.R` — `av_to_asspDataObj()`, `processMediaFiles_LoadAndProcess()`
- `R/pladdrr_helpers.R` — `av_load_for_pladdrr()`
- `R/s7_avaudio.R` — `read_avaudio()`, `as_avaudio()`, `avaudio_to_tempfile()` (S7 internal dispatch)
- `R/prep_recode.R` — `prep_recode()` (in-memory transcoding)

**S7 AVAudio class** (internal): used for memory-only dispatch in `s7_methods.R`. The class constructor and helpers are NOT exported; users obtain audio data only via `read_audio()` returning an AsspDataObj.

## Dependencies

**Required** (Imports):
- av, Rcpp, S7, parallel, cli, rlang
- tidyr, assertthat, readr, stringr, tools, digest, logger, uuid, R.matlab, dplyr, purrr
- pladdrr (for Praat-backed algorithms)

**System**: C++11 compiler (gcc, clang, MSVC)

**Bundled**:
- ASSP C library (`src/assp/`)
- SPTK (`src/SPTK/` submodule)
- ESTk (`src/ESTK/` submodule)
- OpenSMILE C++ (`src/opensmile/`, editable)
- COVAREP C++ (bundled)

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

**Function not exported**: Add `@export` to roxygen2 comment (only for `trk_*`, `lst_*`, `ucnv_*`, `read_*`, `write_*` per Export Policy below), then `devtools::document()`

**av::read_audio_bin fails**: Check av package installed + FFmpeg available (`av::av_video_info()`)

**Tests timeout**: Disable parallel: `trk_pitch_rapt(..., parallel = FALSE)`

## Export Policy (strict)

User-exportable **only**:
- `trk_*`, `lst_*` — DSP functions
- `ucnv_*` — Unit conversion
- `read_*`, `write_*` — I/O (`read_ssff`, `read_audio`, `read_jstf`, `read_track`, `write_ssff`, `write_jstf`, `write_track`)
- S3 generics on data classes: `sample_rate`, `n_records`, `signal_duration`, `start_time`, `track_names`, `file_path`, `track_formats`
  (deprecated aliases still exported for 2.8.x compat: `rate`, `numRecs`, `dur`, `startTime`, `tracks`)

**Not exported** (internal — use `inherits()`, internal access via `:::` for tests):
- Class constructors (`AVAudio`)
- Type predicates (`is.AsspDataObj`)
- All `*_cpp` Rcpp bindings
- All audio-loader helpers (`av_to_asspDataObj`, `read_avaudio`, `prep_recode`, …)
- All format-helper / track-helper / validation utilities

The export-policy unit test (`tests/testthat/test-export-policy.R`) enforces this surface — it greps `getNamespaceExports("superassp")` against the allowed regex.

Class methods (S3 `print`, `summary`, `as.data.frame`, `as_tibble`, `cut`, plus the 5 generics above) are documented **with the class** via `@describeIn` / `@rdname`, not on standalone Rd pages.

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

## graphify

This project has a graphify knowledge graph at graphify-out/.

Rules:
- Before answering architecture or codebase questions, read graphify-out/GRAPH_REPORT.md for god nodes and community structure
- If graphify-out/wiki/index.md exists, navigate it instead of reading raw files
- For cross-module "how does X relate to Y" questions, prefer `graphify query "<question>"`, `graphify path "<A>" "<B>"`, or `graphify explain "<concept>"` over grep — these traverse the graph's EXTRACTED + INFERRED edges instead of scanning files
- After modifying code files in this session, run `graphify update .` to keep the graph current (AST-only, no API cost)

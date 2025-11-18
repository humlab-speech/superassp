# DSP Function Naming Convention Analysis

## Proposed Naming Convention

1. **SSFF-outputting functions**: Prefix with `trk_` (for "track")
2. **List-outputting functions**: Prefix with `lst_` 
3. **Framework identifiers**: Remove `praat`, `sptk`, `estk` from names
4. **Exception**: Keep `praat_sauce` as-is (special case)
5. **Parselmouth functions**: Add `p` suffix instead of `praat` prefix
6. **VAT functions**: Deprecate `vat` and `voice_analysis_toolkit` as they will be replaced

## Current Function Inventory

### SSFF Functions (outputType = "SSFF" or "EST_Track")

#### ASSP C Library Functions
| Current Name | Implementation | Proposed Name | Notes |
|-------------|----------------|---------------|-------|
| `acfana` | ASSP C | `trk_acf` | Auto-correlation analysis |
| `arfana` | ASSP C | `trk_arf` | AR formant analysis |
| `cepstrum` | ASSP C | `trk_cepstrum` | Cepstral analysis |
| `cssSpectrum` | ASSP C | `trk_css` | CSS spectrum |
| `dftSpectrum` | ASSP C | `trk_dft` | DFT spectrum |
| `forest` | ASSP C | `trk_forest` | Formant estimation |
| `ksvfo` | ASSP C | `trk_ksvf0` | KSV F0 tracking |
| `larana` | ASSP C | `trk_lar` | LAR analysis |
| `lpcana` | ASSP C | `trk_lpc` | LPC analysis |
| `lpsSpectrum` | ASSP C | `trk_lps` | LPS spectrum |
| `mhspitch` | ASSP C | `trk_mhsf0` | Modified Harmonic Sieve pitch |
| `rfcana` | ASSP C | `trk_rfc` | Reflection coefficient analysis |
| `rmsana` | ASSP C | `trk_rms` | RMS analysis |
| `zcrana` | ASSP C | `trk_zcr` | Zero-crossing rate |

**Legacy ASSP names (aliases):**
- `fo` → Already alias, remove or keep as-is
- `fo_ksv` → Already alias, remove or keep as-is
- `foana` → Already alias, remove or keep as-is
- `pitch` → Already alias, remove or keep as-is
- `pitch_mhs` → Already alias, remove or keep as-is

#### SPTK C++ Functions (Direct C++ Calls)
| Current Name | Implementation | Proposed Name | Notes |
|-------------|----------------|---------------|-------|
| `rapt` | SPTK C++ | `trk_rapt` | RAPT pitch tracker |
| `swipe` | SPTK C++ | `trk_swipe` | SWIPE pitch tracker |
| `reaper` | SPTK C++ | `trk_reaper` | REAPER pitch tracker |
| `dio` | SPTK C++ | `trk_dio` | DIO pitch tracker |
| `harvest` | SPTK C++ | `trk_harvest` | HARVEST pitch tracker |
| `d4c` | SPTK C++ | `trk_d4c` | D4C aperiodicity |
| `sptk_mfcc` | SPTK C++ | `trk_mfcc` | MFCC from SPTK (remove sptk_ prefix) |

#### SPTK Python Functions (Deprecated - via pysptk)
| Current Name | Implementation | Proposed Name | Notes |
|-------------|----------------|---------------|-------|
| `nonopt_rapt` | SPTK Python | **DEPRECATE** | Superseded by C++ version |
| `nonopt_swipe` | SPTK Python | **DEPRECATE** | Superseded by C++ version |
| `nonopt_reaper` | SPTK Python | **DEPRECATE** | Superseded by C++ version |
| `aperiodicities` | SPTK Python | **DEPRECATE** | Superseded by d4c C++ version |

#### ESTK C++ Functions
| Current Name | Implementation | Proposed Name | Notes |
|-------------|----------------|---------------|-------|
| `estk_pitchmark` | ESTK C++ | `trk_pitchmark` | Pitchmark extraction (remove estk_ prefix) |

#### Snack C Functions
| Current Name | Implementation | Proposed Name | Notes |
|-------------|----------------|---------------|-------|
| `snack_pitch` | Snack C | `trk_snackf0` | Keep "snack" to identify algorithm |
| `snack_formant` | Snack C | `trk_snackformant` | Keep "snack" to identify algorithm |

#### Parselmouth/Praat Python Functions
| Current Name | Implementation | Proposed Name | Notes |
|-------------|----------------|---------------|-------|
| `praat_pitch` | Parselmouth | `trk_f0p` | Praat pitch → pitch with 'p' suffix |
| `praat_formant_burg` | Parselmouth | `trk_formantp` | Praat formant → formant with 'p' |
| `praat_formantpath_burg` | Parselmouth | `trk_formantpathp` | Formant path with 'p' |
| `praat_intensity` | Parselmouth | `trk_intensityp` | Intensity with 'p' |
| `praat_spectral_moments` | Parselmouth | `trk_spectralmomentsp` | Spectral moments with 'p' |
| `praat_sauce` | Parselmouth | `praat_sauce` | **KEEP AS-IS** (exception) |

#### Other Python-based Functions
| Current Name | Implementation | Proposed Name | Notes |
|-------------|----------------|---------------|-------|
| `crepe` | Python/CREPE | `trk_crepe` | CREPE pitch tracker |
| `pyin` | Python/pYIN | `trk_pyin` | pYIN pitch tracker |
| `yin` | Python/YIN | `trk_yin` | YIN pitch tracker |
| `yaapt` | Python/YAAPT | `trk_yaapt` | YAAPT pitch tracker |
| `seenc` | Python/SEENC | `trk_seenc` | SEENC pitch tracker |
| `excite` | Python | `trk_excite` | Excitation analysis |
| `harmonics` | Python | `trk_harmonics` | Harmonic analysis |
| `reaper_pm` | Parselmouth | `trk_reaper_pm` | REAPER via Parselmouth |
| `kaldi_pitch` | Torch/Kaldi | `trk_kaldif0` | Kaldi pitch |
| `torch_pitch` | Torch | `trk_torchf0` | Torch pitch |
| `mfcc` | Torch | **DEPRECATE** | Superseded by sptk_mfcc |
| `npy_import` | NumPy import | `trk_npy` | NumPy data import |

### List Functions (outputType = "list")

#### OpenSMILE Functions
| Current Name | Implementation | Proposed Name | Notes |
|-------------|----------------|---------------|-------|
| `ComParE_2016` | OpenSMILE | `lst_ComParE2016` | Feature extraction |
| `eGeMAPS` | OpenSMILE | `lst_eGeMAPS` | Feature extraction |
| `emobase` | OpenSMILE | `lst_emobase` | Feature extraction |
| `GeMAPS` | OpenSMILE | `lst_GeMAPS` | Feature extraction |

#### Parselmouth/Praat Functions
| Current Name | Implementation | Proposed Name | Notes |
|-------------|----------------|---------------|-------|
| `praat_voice_report` | Parselmouth | `lst_voicereportp` | Voice report with 'p' |
| `praat_voice_tremor` | Parselmouth | `lst_voicetremorp` | Voice tremor with 'p' |
| `praat_avqi` | Parselmouth | `lst_avqip` | AVQI with 'p' |
| `praat_dsi` | Parselmouth | `lst_dsip` | DSI with 'p' |

#### MATLAB/Other Functions
| Current Name | Implementation | Proposed Name | Notes |
|-------------|----------------|---------------|-------|
| `voice_analysis_toolkit` | MATLAB | `lst_voiceanalysistoolkit` | Or `lstvat` for brevity |

## Summary of Changes

### SSFF Functions (→ `trk` prefix)
- **Total**: ~50 functions
- **Requiring rename**: ~47 functions
- **Keep as-is**: `praat_sauce` (1 function)
- **To deprecate**: 4 functions (nonopt_* and aperiodicities)

### List Functions (→ `lst` prefix)
- **Total**: 9 functions
- **Requiring rename**: 9 functions

### Key Patterns

1. **Remove framework prefixes**: `praat_`, `sptk_`, `estk_` removed (except praat_sauce)
2. **Add output type prefix**: `trk` or `lst`
3. **Parselmouth functions**: Add `p` suffix (e.g., `trkf0p` for praat_pitch)
4. **Keep algorithm identifiers**: `snack`, `crepe`, `kaldi`, etc. preserved for clarity
5. **F0 notation**: Consider standardizing pitch functions with `f0` suffix (e.g., `trkrapt` could be `trkraptf0`)

## Implementation Impact

### Files to Modify
- All R files in `R/ssff_*.R` and `R/list_*.R`
- Documentation in `man/*.Rd`
- Tests in `tests/testthat/test-*.R`
- Benchmarking scripts in `benchmarking/`
- README.md and vignettes
- NAMESPACE file (exports)

### Backward Compatibility Options

1. **Hard break**: Remove old names entirely (cleanest, but breaks existing code)
2. **Deprecation period**: Keep old names as deprecated aliases
3. **Dual naming**: Support both old and new names indefinitely

### Recommended Approach

Create **deprecated aliases** with warnings for 1-2 package versions before removal:

```r
# Example for rapt
rapt <- function(...) {
  .Deprecated("trkrapt", package = "superassp",
              msg = "rapt() is deprecated. Use trkrapt() instead.")
  trkrapt(...)
}
```

## Questions for Decision

1. **F0 vs pitch**: Should pitch trackers use `f0` in name? (e.g., `trkraptf0` vs `trkrapt`)
2. **Abbreviations**: How aggressive? (e.g., `trkvoicereportp` vs `lstvrepp`)
3. **Snack**: Keep "snack" prefix or remove? (Recommendation: keep for algorithm identification)
4. **Deprecation timeline**: How many versions to support old names?
5. **voice_analysis_toolkit**: Rename to `lstvoiceanalysistoolkit` or abbreviate to `lstvat`?

## Estimated Effort

- **Code changes**: ~60 function definitions + tests + docs
- **Test updates**: All test files referencing these functions
- **Documentation**: All .Rd files + README + vignettes
- **Time estimate**: 4-6 hours for systematic renaming + testing

## Quick Reference Table

### Pitch/F0 Tracking Functions

| Current Name | New Name | Framework | Status |
|--------------|----------|-----------|--------|
| `rapt` | `trk_rapt` | SPTK C++ | Rename |
| `swipe` | `trk_swipe` | SPTK C++ | Rename |
| `reaper` | `trk_reaper` | SPTK C++ | Rename |
| `dio` | `trk_dio` | SPTK C++ | Rename |
| `harvest` | `trk_harvest` | SPTK C++ | Rename |
| `mhspitch` | `trk_mhsf0` | ASSP C | Rename |
| `ksvfo` | `trk_ksvf0` | ASSP C | Rename |
| `praat_pitch` | `trk_f0p` | Parselmouth | Rename |
| `snack_pitch` | `trk_snackf0` | Snack C | Rename |
| `crepe` | `trk_crepe` | Python | Rename |
| `pyin` | `trk_pyin` | Python | Rename |
| `yin` | `trk_yin` | Python | Rename |
| `yaapt` | `trk_yaapt` | Python | Rename |
| `seenc` | `trk_seenc` | Python | Rename |
| `kaldi_pitch` | `trk_kaldif0` | Torch | Rename |
| `torch_pitch` | `trk_torchf0` | Torch | Rename |
| `reaper_pm` | `trk_reaper_pm` | Parselmouth | Rename |

### Formant Functions

| Current Name | New Name | Framework | Status |
|--------------|----------|-----------|--------|
| `forest` | `trk_forest` | ASSP C | Rename |
| `arfana` | `trk_arf` | ASSP C | Rename |
| `praat_formant_burg` | `trk_formantp` | Parselmouth | Rename |
| `praat_formantpath_burg` | `trk_formantpathp` | Parselmouth | Rename |
| `snack_formant` | `trk_snackformant` | Snack C | Rename |

### Spectral Analysis Functions

| Current Name | New Name | Framework | Status |
|--------------|----------|-----------|--------|
| `cepstrum` | `trk_cepstrum` | ASSP C | Rename |
| `dftSpectrum` | `trk_dft` | ASSP C | Rename |
| `lpsSpectrum` | `trk_lps` | ASSP C | Rename |
| `cssSpectrum` | `trk_css` | ASSP C | Rename |
| `sptk_mfcc` | `trk_mfcc` | SPTK C++ | Rename |
| `praat_spectral_moments` | `trk_spectralmomentsp` | Parselmouth | Rename |

### Other Signal Analysis

| Current Name | New Name | Framework | Status |
|--------------|----------|-----------|--------|
| `acfana` | `trk_acf` | ASSP C | Rename |
| `lpcana` | `trk_lpc` | ASSP C | Rename |
| `larana` | `trk_lar` | ASSP C | Rename |
| `rfcana` | `trk_rfc` | ASSP C | Rename |
| `rmsana` | `trk_rms` | ASSP C | Rename |
| `zcrana` | `trk_zcr` | ASSP C | Rename |
| `praat_intensity` | `trk_intensityp` | Parselmouth | Rename |
| `d4c` | `trk_d4c` | SPTK C++ | Rename |
| `estk_pitchmark` | `trk_pitchmark` | ESTK C++ | Rename |
| `excite` | `trk_excite` | Python | Rename |
| `harmonics` | `trk_harmonics` | Python | Rename |
| `npy_import` | `trk_npy` | NumPy | Rename |
| `praat_sauce` | `praat_sauce` | Parselmouth | **KEEP** |

### Feature Extraction (List Output)

| Current Name | New Name | Framework | Status |
|--------------|----------|-----------|--------|
| `ComParE_2016` | `lst_ComParE2016` | OpenSMILE | Rename |
| `eGeMAPS` | `lst_eGeMAPS` | OpenSMILE | Rename |
| `emobase` | `lst_emobase` | OpenSMILE | Rename |
| `GeMAPS` | `lst_GeMAPS` | OpenSMILE | Rename |
| `praat_voice_report` | `lst_voicereportp` | Parselmouth | Rename |
| `praat_voice_tremor` | `lst_voicetremorp` | Parselmouth | Rename |
| `praat_avqi` | `lst_avqip` | Parselmouth | Rename |
| `praat_dsi` | `lst_dsip` | Parselmouth | Rename |
| `voice_analysis_toolkit` | `lst_voiceanalysistoolkit` | MATLAB | Rename |

### Functions to Deprecate

| Current Name | Reason | Superseded By |
|--------------|--------|---------------|
| `nonopt_rapt` | Slow Python version | `trk_rapt` (C++) |
| `nonopt_swipe` | Slow Python version | `trk_swipe` (C++) |
| `nonopt_reaper` | Slow Python version | `trk_reaper` (C++) |
| `aperiodicities` | Slow Python version | `trk_d4c` (C++) |
| `mfcc` | Slow Torch version | `trk_mfcc` (SPTK C++) |

## Implementation Checklist

### Phase 1: Preparation
- [ ] Create deprecation helper function
- [ ] Set up automated testing for renamed functions
- [ ] Document migration path in vignette
- [ ] Update version to 2.0.0 in DESCRIPTION

### Phase 2: Core Renaming (SSFF Functions)
- [ ] Rename ASSP functions (14 functions)
- [ ] Rename SPTK functions (7 functions)
- [ ] Rename Parselmouth functions (6 functions)
- [ ] Rename Python functions (10+ functions)
- [ ] Rename Snack functions (2 functions)
- [ ] Rename ESTK functions (1 function)

### Phase 3: List Functions
- [ ] Rename OpenSMILE functions (4 functions)
- [ ] Rename Parselmouth list functions (4 functions)
- [ ] Rename other list functions (1 function)

### Phase 4: Deprecation
- [ ] Add deprecation warnings to old function names
- [ ] Mark deprecated functions in documentation
- [ ] Update NEWS.md with deprecation notices

### Phase 5: Documentation
- [ ] Update all .Rd files in man/
- [ ] Update README.md
- [ ] Update vignettes
- [ ] Update CLAUDE.md
- [ ] Create migration guide

### Phase 6: Testing
- [ ] Update all test files
- [ ] Run full test suite
- [ ] Update benchmarking scripts
- [ ] Verify backward compatibility with deprecated aliases

### Phase 7: Release
- [ ] Update NEWS.md with complete change list
- [ ] Tag release as 2.0.0
- [ ] Update package website
- [ ] Announce breaking changes to users


# Function Renaming Implementation Plan

## Overview

This document outlines the comprehensive renaming of DSP functions in the superassp package according to the new naming convention:

- **SSFF output functions**: `trk_` prefix
- **List output functions**: `lst_` prefix
- **Remove framework names**: Remove "praat", "sptk", "estk" from function names (except `praat_sauce`)
- **Parselmouth functions**: Add `p` suffix (e.g., `trk_pitchp`)
- **Deprecate redundant functions**: Functions with more efficient implementations

## Summary Statistics

- **Total functions to rename**: 44
- **SSFF functions (trk_)**: 36
- **List functions (lst_)**: 8
- **Functions to deprecate**: 8
  - 6 SSFF functions (Python versions replaced by C++)
  - 2 list functions (VAT toolkit)

## Detailed Mapping

### SSFF Functions (trk_ prefix)

#### ASSP Library Functions
| Old Name | New Name | Status |
|----------|----------|--------|
| `acfana` | `trk_acfana` | Rename |
| `cepstrum` | `trk_cepstrum` | Rename |
| `cssSpectrum` | `trk_cssSpectrum` | Rename |
| `dftSpectrum` | `trk_dftSpectrum` | Rename |
| `forest` | `trk_forest` | Rename |
| `ksvfo` | `trk_ksvfo` | Rename |
| `lpsSpectrum` | `trk_lpsSpectrum` | Rename |
| `mhspitch` | `trk_mhspitch` | Rename |
| `rmsana` | `trk_rmsana` | Rename |
| `zcrana` | `trk_zcrana` | Rename |

#### SPTK C++ Functions (remove sptk prefix)
| Old Name | New Name | Status |
|----------|----------|--------|
| `rapt` | `trk_rapt` | Rename |
| `swipe` | `trk_swipe` | Rename |
| `reaper` | `trk_reaper` | Rename |
| `dio` | `trk_dio` | Rename |
| `harvest` | `trk_harvest` | Rename |
| `d4c` | `trk_d4c` | Rename |
| `sptk_mfcc` | `trk_mfcc` | Rename |

#### ESTK Functions (remove estk prefix)
| Old Name | New Name | Status |
|----------|----------|--------|
| `estk_pitchmark` | `trk_pitchmark` | Rename |

#### Parselmouth Functions (add p suffix, remove praat prefix)
| Old Name | New Name | Status |
|----------|----------|--------|
| `praat_pitch` | `trk_pitchp` | Rename |
| `praat_formant_burg` | `trk_formantp` | Rename |
| `praat_formantpath_burg` | `trk_formantpathp` | Rename |
| `praat_intensity` | `trk_intensityp` | Rename |
| `praat_spectral_moments` | `trk_spectral_momentsp` | Rename |
| `praat_sauce` | `trk_praat_sauce` | Rename (keep praat_sauce) |

#### Python-Based Pitch Trackers
| Old Name | New Name | Status |
|----------|----------|--------|
| `crepe` | `trk_crepe` | Rename |
| `pyin` | `trk_pyin` | Rename |
| `yin` | `trk_yin` | Rename |
| `yaapt` | `trk_yaapt` | Rename |
| `kaldi_pitch` | `trk_kaldi_pitch` | Rename |
| `torch_pitch` | `trk_torch_pitch` | Rename |
| `seenc` | `trk_seenc` | Rename |
| `excite` | `trk_excite` | Rename |

#### Snack Functions
| Old Name | New Name | Status |
|----------|----------|--------|
| `snack_pitch` | `trk_snackp` | Rename |
| `snack_formant` | `trk_snackf` | Rename |

#### Other SSFF Functions
| Old Name | New Name | Status |
|----------|----------|--------|
| `aperiodicities` | `trk_aperiodicities` | Rename |
| `npy_import` | `trk_npy_import` | Rename |

#### SSFF Functions to Deprecate
| Old Name | Reason | Replacement |
|----------|--------|-------------|
| `nonopt_rapt` | Python version slower | `trk_rapt` |
| `nonopt_swipe` | Python version slower | `trk_swipe` |
| `nonopt_reaper` | Python version slower | `trk_reaper` |
| `reaper_pm` | Duplicate implementation | `trk_reaper` |
| `dio_python` | Python version slower | `trk_dio` |
| `mfcc` | Torch version slower | `trk_mfcc` |

### List Functions (lst_ prefix)

#### OpenSMILE Feature Sets
| Old Name | New Name | Status |
|----------|----------|--------|
| `ComParE_2016` | `lst_ComParE_2016` | Rename |
| `GeMAPS` | `lst_GeMAPS` | Rename |
| `eGeMAPS` | `lst_eGeMAPS` | Rename |
| `emobase` | `lst_emobase` | Rename |

#### Parselmouth List Functions (add p suffix)
| Old Name | New Name | Status |
|----------|----------|--------|
| `praat_avqi` | `lst_avqip` | Rename |
| `praat_dsi` | `lst_dsip` | Rename |
| `praat_voice_report` | `lst_voice_reportp` | Rename |
| `praat_voice_tremor` | `lst_voice_tremorp` | Rename |

#### List Functions to Deprecate
| Old Name | Reason |
|----------|--------|
| `voice_analysis_toolkit` | Will be replaced with better implementation |
| `vat` | Alias for voice_analysis_toolkit |

### Low-Level C++ Functions (No Changes)

These internal functions maintain their `_cpp` suffix and are not renamed:
- `rapt_cpp`, `swipe_cpp`, `reaper_cpp`, `dio_cpp`
- `harvest_cpp`, `d4c_cpp`, `sptk_mfcc_cpp`
- `estk_pitchmark_cpp`, `estk_pda_cpp`

## Implementation Steps

### Phase 1: Function Implementations (R/)
For each renamed function:
1. Copy function definition with new name
2. Add deprecation wrapper for old name
3. Update internal references
4. Update function attributes (ext, tracks, outputType)

### Phase 2: NAMESPACE
1. Export new function names
2. Keep old names exported temporarily with deprecation warnings
3. Update any S3 method exports if needed

### Phase 3: Documentation (man/)
1. Create new .Rd files for renamed functions
2. Add alias links from old to new names
3. Add deprecation notes to old function docs
4. Update cross-references

### Phase 4: Tests (tests/testthat/)
Files to update:
- `test-sptk-pitch.R` → Update to new names
- `test-snack.R` → Update to new names
- `test-kaldi-torch.R` → Update to new names
- `test-*.R` → Update all function calls

### Phase 5: Benchmarks (benchmarking/)
Files to update:
- `benchmark_suite.R` → Update function calls
- `benchmark_pitch.R` → Update pitch function names
- `benchmark_formant.R` → Update formant function names
- Individual benchmark scripts

### Phase 6: README.md
Update all documentation:
- Function name tables
- Code examples
- Benchmark results tables
- Quick reference sections

## Deprecation Strategy

For deprecated functions, add:

```r
#' @export
old_function_name <- function(...) {
  .Deprecated("new_function_name", 
              msg = "old_function_name is deprecated. Use new_function_name instead (C++ implementation is faster).")
  new_function_name(...)
}
```

## Impact Analysis

### Breaking Changes
- All 44 function names change
- Existing user code will need updates
- Package documentation fully regenerated

### Backward Compatibility
- Old function names maintained with deprecation warnings
- Gradual transition period before removal
- Clear migration messages

### Benefits
- Consistent naming convention
- Clear distinction between SSFF and list outputs
- Simpler names (framework prefixes removed)
- Clear identification of Parselmouth functions

## Testing Strategy

After renaming:
1. Run full test suite: `devtools::test()`
2. Check package: `devtools::check()`
3. Build documentation: `devtools::document()`
4. Run benchmarks to verify performance unchanged
5. Test deprecated function warnings

## Files Requiring Changes

### R Source Files (36 files)
- All `ssff_*.R` files with renamed functions
- All `list_*.R` files with renamed functions
- Update cross-references

### Test Files (~15 files)
- All `test-*.R` files in `tests/testthat/`

### Documentation (~50 files)
- All `.Rd` files in `man/`
- New files for renamed functions

### Benchmark Files (~10 files)
- All scripts in `benchmarking/`

### Main Documentation
- `README.md`
- Vignettes if any

### Metadata
- `NAMESPACE`
- `DESCRIPTION` (version bump)

## Timeline Estimate

- **Phase 1** (R implementations): 2-3 hours
- **Phase 2** (NAMESPACE): 30 minutes
- **Phase 3** (Documentation): 2-3 hours
- **Phase 4** (Tests): 2-3 hours
- **Phase 5** (Benchmarks): 1-2 hours
- **Phase 6** (README): 1 hour
- **Testing & Validation**: 2 hours

**Total**: ~12-15 hours

## Success Criteria

- [ ] All 44 functions renamed
- [ ] All tests pass with new names
- [ ] All benchmarks run with new names
- [ ] Documentation generated without errors
- [ ] `R CMD check` passes with no errors
- [ ] README updated with new names
- [ ] Deprecation warnings work correctly
- [ ] Old functions still work (with warnings)

## Post-Implementation

1. Create comprehensive migration guide for users
2. Update all examples in package documentation
3. Add NEWS.md entry documenting changes
4. Consider deprecation timeline (e.g., 2-3 releases before removal)

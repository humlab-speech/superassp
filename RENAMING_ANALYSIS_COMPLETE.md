# Naming Convention Analysis - Complete Report

## Executive Summary

The superassp package contains 52 DSP functions that require renaming according to new naming conventions. This analysis provides a complete mapping and implementation strategy.

## Naming Convention Rules

### Core Principles

1. **Output Type Prefixes**
   - `trk_` prefix for functions outputting SSFF files / AsspDataObj (track data)
   - `lst_` prefix for functions outputting R lists (feature measurements)

2. **Framework Name Removal**
   - Remove `praat_`, `sptk_`, `estk_` prefixes from function names
   - **Exception**: Keep `praat_sauce` as-is (explicit requirement)

3. **Implementation Suffixes**
   - Add `p` suffix for Parselmouth-based functions
   - Example: `praat_pitch` → `trk_pitchp`

4. **Deprecation**
   - Deprecate functions with more efficient implementations
   - Deprecate `voice_analysis_toolkit` (VAT) - to be replaced

## Detailed Analysis

### Category 1: ASSP C Library Functions (10 functions)

Native C implementations from the ASSP library. High performance, well-tested.

| Current Name | New Name | Type | File |
|--------------|----------|------|------|
| acfana | trk_acfana | Autocorrelation | ssff_c_assp_acfana.R |
| cepstrum | trk_cepstrum | Cepstral analysis | ssff_c_assp_cepstrum.R |
| cssSpectrum | trk_cssSpectrum | CSS spectrum | ssff_c_assp_cssSpectrum.R |
| dftSpectrum | trk_dftSpectrum | DFT spectrum | ssff_c_assp_dftSpectrum.R |
| forest | trk_forest | Formant tracking | ssff_c_assp_forest.R |
| ksvfo | trk_ksvfo | Pitch (Kiel algorithm) | ssff_c_assp_ksvfo.R |
| lpsSpectrum | trk_lpsSpectrum | LPC spectrum | ssff_c_assp_lpsSpectrum.R |
| mhspitch | trk_mhspitch | Pitch (MHS algorithm) | ssff_c_assp_mhspitch.R |
| rmsana | trk_rmsana | RMS energy | ssff_c_assp_rmsana.R |
| zcrana | trk_zcrana | Zero-crossing rate | ssff_c_assp_zcrana.R |

**Action**: Rename by adding `trk_` prefix. Simple find-replace operation.

### Category 2: SPTK C++ Functions (7 functions)

Direct C++ bindings to SPTK library via Rcpp. Highest performance tier.

| Current Name | New Name | Algorithm | File |
|--------------|----------|-----------|------|
| rapt | trk_rapt | RAPT pitch | ssff_cpp_sptk_rapt.R |
| swipe | trk_swipe | SWIPE pitch | ssff_cpp_sptk_swipe.R |
| reaper | trk_reaper | REAPER pitch | ssff_cpp_sptk_reaper.R |
| dio | trk_dio | DIO pitch | ssff_cpp_sptk_dio.R |
| harvest | trk_harvest | Harvest pitch | ssff_cpp_sptk_harvest.R |
| d4c | trk_d4c | D4C aperiodicity | ssff_cpp_sptk_d4c.R |
| sptk_mfcc | trk_mfcc | MFCC features | ssff_cpp_sptk_mfcc.R |

**Action**: Rename by adding `trk_` prefix and removing `sptk_` prefix.

### Category 3: ESTK C++ Functions (1 function)

Edinburgh Speech Tools Kit implementations.

| Current Name | New Name | Purpose | File |
|--------------|----------|---------|------|
| estk_pitchmark | trk_pitchmark | Pitch marking | ssff_cpp_estk_pitchmark.R |

**Action**: Rename by adding `trk_` prefix and removing `estk_` prefix.

### Category 4: Parselmouth Functions (6 SSFF + 4 list)

Python Parselmouth library wrapping Praat functionality.

#### SSFF Output (6 functions)
| Current Name | New Name | Purpose | File |
|--------------|----------|---------|------|
| praat_pitch | trk_pitchp | Pitch tracking | ssff_python_pm_ppitch.R |
| praat_formant_burg | trk_formantp | Formant tracking (Burg) | ssff_python_pm_pformantb.R |
| praat_formantpath_burg | trk_formantpathp | Formant path tracking | ssff_python_pm_pformantpathb.R |
| praat_intensity | trk_intensityp | Intensity contour | ssff_python_pm_pintensity.R |
| praat_spectral_moments | trk_spectral_momentsp | Spectral moments | ssff_python_pm_pspectral_moments.R |
| praat_sauce | trk_praat_sauce | SAUCE analysis | ssff_python_pm_psauce.R |

#### List Output (4 functions)
| Current Name | New Name | Purpose | File |
|--------------|----------|---------|------|
| praat_avqi | lst_avqip | Acoustic Voice Quality Index | list_python_pm_pavqi.R |
| praat_dsi | lst_dsip | Dysphonia Severity Index | list_python_pm_pdsi.R |
| praat_voice_report | lst_voice_reportp | Voice quality report | list_python_pm_pvoice_report.R |
| praat_voice_tremor | lst_voice_tremorp | Voice tremor analysis | list_python_pm_pvoice_tremor.R |

**Action**: Remove `praat_` prefix, add `trk_`/`lst_` prefix, add `p` suffix (except `praat_sauce`).

### Category 5: Python Pitch Trackers (8 functions)

Various Python-based pitch tracking algorithms.

| Current Name | New Name | Library | File |
|--------------|----------|---------|------|
| crepe | trk_crepe | CREPE (deep learning) | ssff_python_crepe.R |
| pyin | trk_pyin | pYIN (probabilistic YIN) | ssff_python_pyin.R |
| yin | trk_yin | YIN algorithm | ssff_python_yin.R |
| yaapt | trk_yaapt | YAAPT algorithm | ssff_python_yaapt.R |
| kaldi_pitch | trk_kaldi_pitch | Kaldi (via torch) | ssff_python_kaldi_pitch.R |
| torch_pitch | trk_torch_pitch | PyTorch pitch | ssff_python_torch_pitch.R |
| seenc | trk_seenc | SEENC pitch | ssff_python_seenc.R |
| excite | trk_excite | Excite toolkit | ssff_python_excite.R |

**Action**: Add `trk_` prefix.

### Category 6: Snack Functions (2 functions)

Tcl/Snack library implementations for comparison purposes.

| Current Name | New Name | Purpose | File |
|--------------|----------|---------|------|
| snack_pitch | trk_snackp | Pitch tracking | ssff_python_snack_pitch.R |
| snack_formant | trk_snackf | Formant tracking | ssff_python_snack_formant.R |

**Action**: Remove `snack_` prefix, add `trk_` prefix, add descriptive suffix (`p` for pitch, `f` for formant).

### Category 7: Other SSFF Functions (2 functions)

| Current Name | New Name | Purpose | File |
|--------------|----------|---------|------|
| aperiodicities | trk_aperiodicities | Aperiodicity measure | ssff_python_aperiodicities.R |
| npy_import | trk_npy_import | NumPy array import | ssff_python_npy_import.R |

**Action**: Add `trk_` prefix.

### Category 8: OpenSMILE Feature Sets (4 functions)

OpenSMILE-based feature extraction returning lists.

| Current Name | New Name | Feature Set | File |
|--------------|----------|-------------|------|
| ComParE_2016 | lst_ComParE_2016 | Computational Paralinguistics Challenge 2016 | list_python_opensmile_ComParE_2016.R |
| GeMAPS | lst_GeMAPS | Geneva Minimalistic Acoustic Parameter Set | list_python_opensmile_GeMAPS.R |
| eGeMAPS | lst_eGeMAPS | Extended GeMAPS | list_python_opensmile_eGeMAPS.R |
| emobase | lst_emobase | Emotion recognition baseline | list_python_opensmile_emobase.R |

**Action**: Add `lst_` prefix.

### Category 9: Functions to Deprecate (8 functions)

#### Python SPTK Versions (replaced by C++)
| Function | Reason | Replacement | File |
|----------|--------|-------------|------|
| nonopt_rapt | Slower Python implementation | trk_rapt | ssff_python_sptk_rapt.R |
| nonopt_swipe | Slower Python implementation | trk_swipe | ssff_python_sptk_swipe.R |
| nonopt_reaper | Slower Python implementation | trk_reaper | ssff_python_sptk_reaper.R |
| reaper_pm | Duplicate Parselmouth version | trk_reaper | ssff_python_reaper_pm.R |
| dio_python | Slower Python implementation | trk_dio | ssff_python_world_dio.R |

#### Torch MFCC (replaced by SPTK)
| Function | Reason | Replacement | File |
|----------|--------|-------------|------|
| mfcc | Slower Torch implementation | trk_mfcc | ssff_python_torch_mfcc.R |

#### Voice Analysis Toolkit
| Function | Reason | Replacement | File |
|----------|--------|-------------|------|
| voice_analysis_toolkit | Will be replaced | TBD | matlab_slicefunctions.R |
| vat | Alias for VAT | TBD | matlab.R |

**Action**: Add `.Deprecated()` calls directing to replacements.

## Impact Assessment

### Files to Modify

1. **R Source Files**: 44 files
   - Function definitions
   - Internal function calls
   - Function attributes

2. **NAMESPACE**: 1 file
   - 44 export changes
   - Maintain old names temporarily

3. **Documentation (man/)**: ~50 files
   - Create new .Rd files
   - Add aliases for old names
   - Update examples

4. **Tests (tests/testthat/)**: ~15 files
   - Update all function calls
   - Test deprecated warnings

5. **Benchmarks (benchmarking/)**: ~10 files
   - Update benchmark scripts
   - Regenerate plots

6. **README.md**: 1 file
   - Update all tables
   - Update examples
   - Update benchmark results

### Affected Code Patterns

Functions are called in various contexts:

```r
# Direct calls
result <- rapt("audio.wav")

# In pipelines
files %>% map(~rapt(.x, toFile = FALSE))

# With parameters
pitch <- praat_pitch(audio, minPitch = 75, maxPitch = 500)

# Batch processing
results <- rapt(c("file1.wav", "file2.wav"))
```

All need updating to new names.

## Implementation Complexity

### Low Complexity (Simple Rename)
- ASSP functions: Just add `trk_` prefix
- Python pitch trackers: Just add `trk_` prefix
- OpenSMILE: Just add `lst_` prefix

### Medium Complexity (Prefix + Suffix)
- Parselmouth SSFF: Remove `praat_`, add `trk_`, add `p` suffix
- Parselmouth list: Remove `praat_`, add `lst_`, add `p` suffix
- Snack: Remove `snack_`, add `trk_`, add descriptive suffix

### Higher Complexity (Prefix Changes)
- SPTK: Remove `sptk_`, add `trk_`
- ESTK: Remove `estk_`, add `trk_`

### Special Cases
- `praat_sauce`: Keep as `trk_praat_sauce` (name preserved)
- Low-level `_cpp` functions: No changes (internal only)

## Testing Requirements

### Unit Tests
- [ ] All 44 renamed functions pass existing tests
- [ ] Deprecated functions issue warnings
- [ ] Deprecated functions still work (call new versions)

### Integration Tests
- [ ] Batch processing with new names
- [ ] Parallel processing with new names
- [ ] File I/O with new names

### Documentation Tests
- [ ] All examples in .Rd files work
- [ ] README examples work
- [ ] Vignette examples work (if any)

### Performance Tests
- [ ] Benchmarks produce same results
- [ ] No performance regression
- [ ] Plots regenerate correctly

## Risk Assessment

### Low Risk
- Function implementations unchanged (only names)
- Deprecation wrappers prevent breaking changes
- Comprehensive test coverage

### Medium Risk
- Documentation updates (many files)
- Potential for missing cross-references
- User code migration needed

### Mitigation Strategies
1. Keep old names working with warnings
2. Multi-release deprecation period
3. Clear migration guide
4. Automated testing of all functions
5. Examples in documentation use new names

## Success Metrics

1. ✅ All 44 functions renamed
2. ✅ All 8 deprecated functions issue warnings
3. ✅ All tests pass with new names
4. ✅ All documentation generated without errors
5. ✅ `R CMD check` passes with 0 errors, 0 warnings
6. ✅ Benchmarks run successfully
7. ✅ README fully updated
8. ✅ Migration guide created

## Next Steps

1. **Review and Approve**: Confirm mapping is correct
2. **Implement Phase 1**: Rename function implementations
3. **Update Tests**: Modify test files to use new names
4. **Update Documentation**: Regenerate all .Rd files
5. **Update Benchmarks**: Modify benchmark scripts
6. **Update README**: Replace all function names
7. **Validate**: Run full test suite and check
8. **Document Migration**: Create user-facing migration guide

## Appendix: Complete Mapping Table

See `NAMING_CONVENTION_MAPPING.csv` for the complete CSV with all 52 functions.

## Document History

- **Created**: 2025-10-18
- **Purpose**: Comprehensive analysis of function renaming requirements
- **Status**: Analysis complete, ready for implementation

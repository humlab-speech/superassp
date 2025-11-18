# Function Renaming Quick Reference

## Quick Lookup Tables

### Pitch Tracking Functions

| Old Name | New Name | Implementation | Notes |
|----------|----------|----------------|-------|
| `rapt` | `trk_rapt` | SPTK C++ | ✓ Recommended |
| `swipe` | `trk_swipe` | SPTK C++ | ✓ Recommended |
| `reaper` | `trk_reaper` | SPTK C++ | ✓ Recommended |
| `dio` | `trk_dio` | SPTK C++ | ✓ Recommended |
| `harvest` | `trk_harvest` | SPTK C++ | ✓ Recommended |
| `mhspitch` | `trk_mhspitch` | ASSP C | ✓ Recommended |
| `ksvfo` | `trk_ksvfo` | ASSP C | ✓ Recommended |
| `praat_pitch` | `trk_pitchp` | Parselmouth | ✓ Recommended |
| `snack_pitch` | `trk_snackp` | Snack/Tcl | For comparison |
| `crepe` | `trk_crepe` | Python/Deep Learning | High accuracy |
| `pyin` | `trk_pyin` | Python | Probabilistic |
| `yin` | `trk_yin` | Python | Fast |
| `yaapt` | `trk_yaapt` | Python | Robust |
| `kaldi_pitch` | `trk_kaldi_pitch` | Torch/Kaldi | Speech recognition |
| `torch_pitch` | `trk_torch_pitch` | PyTorch | Experimental |
| `seenc` | `trk_seenc` | Python | Research |
| `excite` | `trk_excite` | Python | Research |
| `nonopt_rapt` | **DEPRECATED** | Use `trk_rapt` | |
| `nonopt_swipe` | **DEPRECATED** | Use `trk_swipe` | |
| `nonopt_reaper` | **DEPRECATED** | Use `trk_reaper` | |
| `reaper_pm` | **DEPRECATED** | Use `trk_reaper` | |
| `dio_python` | **DEPRECATED** | Use `trk_dio` | |

### Formant Tracking Functions

| Old Name | New Name | Implementation | Notes |
|----------|----------|----------------|-------|
| `forest` | `trk_forest` | ASSP C | ✓ Recommended |
| `praat_formant_burg` | `trk_formantp` | Parselmouth | ✓ Recommended |
| `praat_formantpath_burg` | `trk_formantpathp` | Parselmouth | Advanced tracking |
| `snack_formant` | `trk_snackf` | Snack/Tcl | For comparison |

### Spectral Analysis Functions

| Old Name | New Name | Implementation |
|----------|----------|----------------|
| `cepstrum` | `trk_cepstrum` | ASSP C |
| `cssSpectrum` | `trk_cssSpectrum` | ASSP C |
| `dftSpectrum` | `trk_dftSpectrum` | ASSP C |
| `lpsSpectrum` | `trk_lpsSpectrum` | ASSP C |
| `praat_spectral_moments` | `trk_spectral_momentsp` | Parselmouth |

### Energy/Intensity Functions

| Old Name | New Name | Implementation |
|----------|----------|----------------|
| `rmsana` | `trk_rmsana` | ASSP C |
| `praat_intensity` | `trk_intensityp` | Parselmouth |
| `zcrana` | `trk_zcrana` | ASSP C |

### Other Signal Processing Functions

| Old Name | New Name | Implementation | Purpose |
|----------|----------|----------------|---------|
| `acfana` | `trk_acfana` | ASSP C | Autocorrelation |
| `estk_pitchmark` | `trk_pitchmark` | ESTK C++ | Pitch marking |
| `d4c` | `trk_d4c` | SPTK C++ | Aperiodicity |
| `aperiodicities` | `trk_aperiodicities` | Python | Aperiodicity |
| `praat_sauce` | `trk_praat_sauce` | Parselmouth | SAUCE voice analysis |
| `npy_import` | `trk_npy_import` | Python | Import NumPy arrays |

### Feature Extraction (MFCC)

| Old Name | New Name | Implementation | Notes |
|----------|----------|----------------|-------|
| `sptk_mfcc` | `trk_mfcc` | SPTK C++ | ✓ Recommended |
| `mfcc` | **DEPRECATED** | Use `trk_mfcc` | Torch version slower |

### OpenSMILE Feature Sets (List Output)

| Old Name | New Name | Purpose |
|----------|----------|---------|
| `ComParE_2016` | `lst_ComParE_2016` | Computational Paralinguistics |
| `GeMAPS` | `lst_GeMAPS` | Geneva Minimalistic Acoustic Parameter Set |
| `eGeMAPS` | `lst_eGeMAPS` | Extended GeMAPS |
| `emobase` | `lst_emobase` | Emotion recognition baseline |

### Voice Quality Analysis (List Output)

| Old Name | New Name | Implementation |
|----------|----------|----------------|
| `praat_avqi` | `lst_avqip` | Parselmouth/AVQI |
| `praat_dsi` | `lst_dsip` | Parselmouth/DSI |
| `praat_voice_report` | `lst_voice_reportp` | Parselmouth |
| `praat_voice_tremor` | `lst_voice_tremorp` | Parselmouth |
| `voice_analysis_toolkit` | **DEPRECATED** | Will be replaced |
| `vat` | **DEPRECATED** | Alias for VAT |

## Migration Examples

### Before (Old Names)
```r
# Pitch tracking
pitch_data <- rapt("audio.wav", toFile = FALSE)
pitch_praat <- praat_pitch("audio.wav", toFile = FALSE)

# Formants
formants <- forest("audio.wav", toFile = FALSE)
formants_p <- praat_formant_burg("audio.wav", toFile = FALSE)

# MFCC
mfcc_data <- sptk_mfcc("audio.wav", toFile = FALSE)

# Voice quality
avqi_results <- praat_avqi("audio.wav")
gemaps <- GeMAPS("audio.wav")
```

### After (New Names)
```r
# Pitch tracking
pitch_data <- trk_rapt("audio.wav", toFile = FALSE)
pitch_praat <- trk_pitchp("audio.wav", toFile = FALSE)

# Formants
formants <- trk_forest("audio.wav", toFile = FALSE)
formants_p <- trk_formantp("audio.wav", toFile = FALSE)

# MFCC
mfcc_data <- trk_mfcc("audio.wav", toFile = FALSE)

# Voice quality
avqi_results <- lst_avqip("audio.wav")
gemaps <- lst_GeMAPS("audio.wav")
```

## Naming Convention Summary

### Prefix Rules
- **`trk_`**: Functions returning SSFF files or AsspDataObj (track data)
- **`lst_`**: Functions returning R lists (feature sets, measurements)

### Suffix Rules
- **`p`**: Parselmouth-based implementations (e.g., `trk_pitchp`, `lst_avqip`)
- **No suffix**: Direct C/C++ or Python implementations

### Framework Names Removed
- ~~`praat_`~~ → Add `p` suffix instead (except `praat_sauce`)
- ~~`sptk_`~~ → Remove prefix
- ~~`estk_`~~ → Remove prefix
- ~~`snack_`~~ → Remove prefix, add descriptive suffix

### Special Cases
- `praat_sauce` → `trk_praat_sauce` (keep full name per rules)
- Low-level C++ functions keep `_cpp` suffix (not renamed)

## Implementation Status

- [x] Analysis completed
- [x] Mapping created (NAMING_CONVENTION_MAPPING.csv)
- [x] Implementation plan created
- [ ] Function implementations updated
- [ ] Tests updated
- [ ] Documentation updated
- [ ] Benchmarks updated
- [ ] README updated

## Notes for Developers

1. **Old names still work** (with deprecation warnings)
2. **Transition period**: Keep both names available for 2-3 releases
3. **Update examples**: All documentation should use new names
4. **Low-level functions**: `_cpp` functions unchanged (internal use)
5. **Backward compatibility**: Deprecation wrappers ensure code doesn't break

## Search & Replace Patterns

For bulk updating user code:

```r
# Pitch functions
rapt\\( → trk_rapt(
swipe\\( → trk_swipe(
reaper\\( → trk_reaper(
dio\\( → trk_dio(
harvest\\( → trk_harvest(
praat_pitch\\( → trk_pitchp(

# Formant functions
forest\\( → trk_forest(
praat_formant_burg\\( → trk_formantp(

# List functions
GeMAPS\\( → lst_GeMAPS(
praat_avqi\\( → lst_avqip(
```

## Questions or Issues

If you encounter any issues with the renaming:
1. Check if the function is in the deprecation list
2. Look up the new name in this guide
3. Update your code to use the new name
4. Report any missing or incorrect mappings to package maintainers

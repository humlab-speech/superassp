# Quick Reference: Function Renaming

## Naming Pattern Summary

| Old Pattern | New Pattern | Example |
|-------------|-------------|---------|
| ASSP functions (`*ana`, etc.) | `trk_*` | `rmsana` → `trk_rms` |
| SPTK functions (`rapt`, etc.) | `trk_*` | `rapt` → `trk_rapt` |
| `sptk_*` | `trk_*` | `sptk_mfcc` → `trk_mfcc` |
| `estk_*` | `trk_*` | `estk_pitchmark` → `trk_pitchmark` |
| `praat_*` (SSFF) | `trk_*p` | `praat_pitch` → `trk_pitchp` |
| `praat_*` (list) | `lst_*p` | `praat_avqi` → `lst_avqip` |
| OpenSmile | `lst_*` | `eGeMAPS` → `lst_egemaps` |
| `snack_*` | `trk_snack_*` | `snack_pitch` → `trk_snack_pitch` |

## Most Commonly Used Functions

### Pitch Tracking
```r
# Before
rapt(file)
swipe(file)
reaper(file)
dio(file)
harvest(file)
praat_pitch(file)
ksvfo(file)
mhspitch(file)

# After
trk_rapt(file)
trk_swipe(file)
trk_reaper(file)
trk_dio(file)
trk_harvest(file)
trk_pitchp(file)
trk_ksvf0(file)
trk_mhsf0(file)
```

### Formant Analysis
```r
# Before
forest(file)
praat_formant_burg(file)
praat_formantpath_burg(file)
snack_formant(file)

# After
trk_forest(file)
trk_formant_burgp(file)
trk_formantpath_burgp(file)
trk_snack_formant(file)
```

### Spectral Analysis
```r
# Before
sptk_mfcc(file)
praat_spectral_moments(file)
cepstrum(file)
dftSpectrum(file)

# After
trk_mfcc(file)
trk_spectral_momentsp(file)
trk_cepstrum(file)
trk_dft_spectrum(file)
```

### Energy/Intensity
```r
# Before
rmsana(file)
praat_intensity(file)

# After
trk_rms(file)
trk_intensityp(file)
```

### Voice Quality Features
```r
# Before
praat_voice_report(file)
praat_avqi(sv_df)
praat_dsi(soft_df)
d4c(file)

# After
lst_voice_reportp(file)
lst_avqip(sv_df)
lst_dsip(soft_df)
trk_d4c(file)
```

### Feature Extraction
```r
# Before
eGeMAPS(file)
GeMAPS(file)
emobase(file)
ComParE_2016(file)

# After
lst_egemaps(file)
lst_gemaps(file)
lst_emobase(file)
lst_compare_2016(file)
```

## Functions to Stop Using (Deprecated)

These have faster replacements:

```r
# DEPRECATED - Use C++ versions instead
nonopt_rapt(file)    # Use: trk_rapt(file)
nonopt_swipe(file)   # Use: trk_swipe(file)
nonopt_reaper(file)  # Use: trk_reaper(file)
aperiodicities(file) # Use: trk_d4c(file)
mfcc(file)           # Use: trk_mfcc(file)

# DEPRECATED - Will be replaced
vat(file)
voice_analysis_toolkit(file)
```

## Complete Alphabetical List

### A-C
| Old | New |
|-----|-----|
| `acfana` | `trk_acf` |
| `aperiodicities` | **DEPRECATED** → use `trk_d4c` |
| `arfana` | `trk_arf` |
| `cepstrum` | `trk_cepstrum` |
| `ComParE_2016` | `lst_compare_2016` |
| `crepe` | `trk_crepe` |
| `cssSpectrum` | `trk_css_spectrum` |

### D-E
| Old | New |
|-----|-----|
| `d4c` | `trk_d4c` |
| `dftSpectrum` | `trk_dft_spectrum` |
| `dio` | `trk_dio` |
| `eGeMAPS` | `lst_egemaps` |
| `emobase` | `lst_emobase` |
| `estk_pitchmark` | `trk_pitchmark` |
| `excite` | `trk_excite` |

### F-H
| Old | New |
|-----|-----|
| `forest` | `trk_forest` |
| `GeMAPS` | `lst_gemaps` |
| `harmonics` | `trk_harmonics` |
| `harvest` | `trk_harvest` |

### K-L
| Old | New |
|-----|-----|
| `kaldi_pitch` | `trk_kaldi_pitch` |
| `ksvfo` | `trk_ksvf0` |
| `larana` | `trk_lar` |
| `lpcana` | `trk_lpc` |
| `lpsSpectrum` | `trk_lps_spectrum` |

### M-N
| Old | New |
|-----|-----|
| `mfcc` | **DEPRECATED** → use `trk_mfcc` |
| `mhspitch` | `trk_mhsf0` |
| `nonopt_rapt` | **DEPRECATED** → use `trk_rapt` |
| `nonopt_reaper` | **DEPRECATED** → use `trk_reaper` |
| `nonopt_swipe` | **DEPRECATED** → use `trk_swipe` |
| `npy_import` | `trk_npy_import` |

### P
| Old | New |
|-----|-----|
| `praat_avqi` | `lst_avqip` |
| `praat_dsi` | `lst_dsip` |
| `praat_formant_burg` | `trk_formant_burgp` |
| `praat_formantpath_burg` | `trk_formantpath_burgp` |
| `praat_intensity` | `trk_intensityp` |
| `praat_pitch` | `trk_pitchp` |
| `praat_sauce` | `praat_sauce` *(unchanged)* |
| `praat_spectral_moments` | `trk_spectral_momentsp` |
| `praat_voice_report` | `lst_voice_reportp` |
| `praat_voice_tremor` | `lst_voice_tremorp` |
| `pyin` | `trk_pyin` |

### R-S
| Old | New |
|-----|-----|
| `rapt` | `trk_rapt` |
| `reaper` | `trk_reaper` |
| `reaper_pm` | `trk_reaperp` |
| `rfcana` | `trk_rfc` |
| `rmsana` | `trk_rms` |
| `seenc` | `trk_seenc` |
| `snack_formant` | `trk_snack_formant` |
| `snack_pitch` | `trk_snack_pitch` |
| `sptk_mfcc` | `trk_mfcc` |
| `swipe` | `trk_swipe` |

### T-Z
| Old | New |
|-----|-----|
| `torch_pitch` | `trk_torch_pitch` |
| `vat` | **DEPRECATED** |
| `voice_analysis_toolkit` | **DEPRECATED** |
| `yaapt` | `trk_yaapt` |
| `yin` | `trk_yin` |
| `zcrana` | `trk_zcr` |

## Search & Replace Patterns

For automated migration of your code:

```r
# Pitch functions
s/\brapt\(/trk_rapt(/g
s/\bswipe\(/trk_swipe(/g
s/\breaper\(/trk_reaper(/g
s/\bdio\(/trk_dio(/g
s/\bharvest\(/trk_harvest(/g

# Praat functions (SSFF)
s/\bpraat_pitch\(/trk_pitchp(/g
s/\bpraat_formant_burg\(/trk_formant_burgp(/g
s/\bpraat_intensity\(/trk_intensityp(/g

# Praat functions (list)
s/\bpraat_voice_report\(/lst_voice_reportp(/g
s/\bpraat_avqi\(/lst_avqip(/g

# ASSP functions
s/\brmsana\(/trk_rms(/g
s/\bacfana\(/trk_acf(/g
s/\bzcrana\(/trk_zcr(/g
s/\bksvfo\(/trk_ksvf0(/g
s/\bmhspitch\(/trk_mhsf0(/g

# OpenSmile
s/\beGeMAPS\(/lst_egemaps(/g
s/\bGeMAPS\(/lst_gemaps(/g
s/\bemobase\(/lst_emobase(/g
s/\bComParE_2016\(/lst_compare_2016(/g

# SPTK
s/\bsptk_mfcc\(/trk_mfcc(/g
s/\bd4c\(/trk_d4c(/g

# Deprecations (with warnings)
s/\bnonopt_rapt\(/trk_rapt(/g   # Update to faster version
s/\bnonopt_swipe\(/trk_swipe(/g
s/\bnonopt_reaper\(/trk_reaper(/g
s/\baperiodicities\(/trk_d4c(/g
s/\bmfcc\(/trk_mfcc(/g  # Use SPTK version
```

## Migration Checklist

- [ ] Update function calls in R scripts
- [ ] Update function calls in Rmarkdown documents
- [ ] Update tests
- [ ] Update benchmarks
- [ ] Update documentation
- [ ] Update citations/references
- [ ] Test that everything still works
- [ ] Update version numbers if needed

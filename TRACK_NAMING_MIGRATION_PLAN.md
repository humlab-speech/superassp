# Track Naming Migration Plan: Titze 2015 + Nylén 2024 Standards

## Executive Summary

This document outlines a comprehensive plan to migrate all DSP function output track names in the superassp package to comply with:

1. **Titze et al. (2015)**: "Toward a consensus on symbolic notation of harmonics, resonances, and formants in vocalization" (JASA)
2. **Nylén et al. (2024)**: Extensions for corrected formants in gender/voice quality analysis (JASA)
3. **Existing standards**: VoiceSauce notation, ISO units integration

**Goal**: Standardize all `attr(result, "tracks")` across ~149 track definitions without breaking existing code.

## Key Standards from Titze 2015

### 1. Fundamental Frequency (fo)
- **Standard notation**: `fo` (lowercase) for frequency of oscillation
- **Unit**: Hz
- **Symbol**: f_o or f0
- **Current usage**: Mixed (`fo`, `F0`, `fo`, `fo[Hz]`)

### 2. Formants (Resonances)
- **Standard notation**: `Fi` where i = 1, 2, 3, 4, ... (uppercase F, subscript number)
- **Unit**: Hz
- **Examples**: F₁, F₂, F₃, F₄
- **Current usage**: Mixed (`F1`, `F2`, `fm`, `F[Hz]`)

### 3. Formant Bandwidths
- **Standard notation**: `Bi` where i = 1, 2, 3, 4, ...
- **Unit**: Hz
- **Examples**: B₁, B₂, B₃, B₄
- **Current usage**: Mixed (`B1`, `B2`, `bw`, `B[Hz]`)

### 4. Harmonics
- **Standard notation**: `Hi` or `Ai` where i = harmonic number or formant proximity
- **Unit**: dB (amplitude)
- **Examples**:
  - H₁ = first harmonic (usually = fo)
  - H₂ = second harmonic
  - A₁ = harmonic nearest F₁
  - A₂ = harmonic nearest F₂
  - A₃ = harmonic nearest F₃

### 5. Harmonic Difference Measures
- **Standard notation**: `Hi-Hj` or `Hi-Aj`
- **Unit**: dB
- **Examples**:
  - H1-H2 (open quotient correlate)
  - H1-A1 (spectral tilt at F1)
  - H1-A2 (spectral tilt at F2)
  - H1-A3 (spectral tilt at F3)

## Key Extensions from Nylén 2024

### Corrected vs Uncorrected Formants

**Problem**: Raw formant measurements are affected by vocal tract filtering. Correction removes formant influence from harmonic measures.

**Notation**:
- **Uncorrected**: Append `*` (asterisk) or use `u` suffix
  - Examples: `H1*`, `H1u`, `F1*`, `F1u`
- **Corrected**: No asterisk/suffix OR append `c`
  - Examples: `H1`, `H1c`, `F1`, `F1c`

**Recommended for superassp**:
- Use `c` suffix for corrected: `H1c`, `F1c`, `A1c`
- Use `u` suffix for uncorrected: `H1u`, `F1u`, `A1u`
- This is clearer than asterisk in R variable names

### Gender-Specific Formants (Nylén Extension)

**Notation for gender-corrected formants**:
- `F1g`, `F2g`, `F3g` where `g` indicates gender-normalized
- Or: `F1_fem`, `F1_masc` for gender-specific ranges

## Current State Analysis

### Track Naming Patterns Found (149 definitions)

#### Category 1: Fundamental Frequency (fo)
**Current patterns**:
- `fo` (lowercase) - Most Python implementations
- `F0` (uppercase) - Swift-F0, some others
- `fo` (no zero) - ASSP library functions
- `fo[Hz]` (with unit) - ASSP library with explicit unit

**Affected functions** (~30):
- `trk_rapt`, `trk_swipe`, `trk_dio`, `trk_reaper`, `trk_harvest`
- `trk_yin`, `trk_pyin`, `trk_crepe`, `trk_yaapt`, `trk_kaldi_pitch`
- `trk_torch_pitch`, `trk_swiftf0`, `trk_snackp`, `trk_ksvfo`
- `fo`, `foana`, `fo_ksv`, `trk_mhspitch`
- VoiceSauce: `strF0`, `sF0`, `pF0`, `shrF0`

**Proposed standard**: `fo[Hz]`
- Lowercase `fo` (Titze standard)
- Explicit unit `[Hz]` (superassp convention)
- Enables automatic unit conversion via units package

#### Category 2: Formants
**Current patterns**:
- `F1`, `F2`, `F3`, `F4` (uppercase, VoiceSauce)
- `fm` (generic formant frequency)
- `F[Hz]` (ASSP forest - multiple formants in matrix)

**Affected functions** (~15):
- `trk_praat_sauce`: `F1`, `F2`, `F3`
- `trk_forest`: `F[Hz]` (returns matrix with multiple formants)
- `trk_formantp`, `trk_formantpathp`: `fm`
- `trk_snackf`: dynamic based on `n` parameter

**Proposed standards**:
- Single formant: `F1[Hz]`, `F2[Hz]`, `F3[Hz]`, `F4[Hz]`
- Multi-column tracks: `F1[Hz]`, `F2[Hz]`, `F3[Hz]`, `F4[Hz]` (separate columns)
- Generic (unknown number): `Fi[Hz]` with i as column suffix

#### Category 3: Formant Bandwidths
**Current patterns**:
- `B1`, `B2`, `B3` (VoiceSauce)
- `bw` (generic bandwidth)
- `B[Hz]` (ASSP forest - multiple bandwidths)

**Affected functions** (~10):
- `trk_praat_sauce`: `B1`, `B2`, `B3`
- `trk_forest`: `B[Hz]`
- `trk_formantp`, `trk_formantpathp`: `bw`

**Proposed standards**:
- `B1[Hz]`, `B2[Hz]`, `B3[Hz]`, `B4[Hz]`
- Multi-column: separate columns per formant bandwidth

#### Category 4: Harmonics (VoiceSauce)
**Current patterns**:
- `H1`, `H2`, `H4` (first, second, fourth harmonics)
- `H2K` (harmonic at ~2 kHz)
- `H5K` (harmonic at ~5 kHz)
- `A1`, `A2`, `A3` (harmonics nearest F1, F2, F3)

**Affected functions**:
- `trk_praat_sauce`: `H1`, `H2`, `H4`, `A1`, `A2`, `A3`, `H2K`, `H5K`

**Proposed standards**:
- Keep notation but add units: `H1[dB]`, `H2[dB]`, `A1[dB]`, etc.
- Distinguish corrected: `H1c[dB]`, `H1u[dB]` when applicable

#### Category 5: Harmonic Differences (VoiceSauce)
**Current patterns**:
- `H1H2`, `H2H4`, `H1A1`, `H1A2`, `H1A3`
- `H4H2K`, `H2KH5K`
- Corrected versions: `H1H2c`, `H1A1c`, etc.
- Uncorrected: `H1H2u`, `H1A1u`, etc.

**Affected functions**:
- `trk_praat_sauce`: Full set of difference measures

**Proposed standards**:
- Use hyphen: `H1-H2[dB]`, `H1-A1[dB]`, `H1-A2[dB]`
- Corrected: `H1-H2c[dB]`, `H1-A1c[dB]`
- Uncorrected: `H1-H2u[dB]`, `H1-A1u[dB]`
- Clearer notation, matches Titze 2015

#### Category 6: Voice Quality Measures
**Current patterns**:
- `CPP` (Cepstral Peak Prominence)
- `HNR`, `HNR05`, `HNR15`, `HNR25` (Harmonics-to-Noise Ratio)
- `SHR` (Subharmonic-to-Harmonic Ratio)
- `Jitter (local)`, `Shimmer (local)` (Praat-style with parentheses)

**Affected functions**:
- `trk_praat_sauce`: `CPP`, `HNR`, `SHR`, `Energy`
- `lst_voice_reportp`: Jitter/shimmer with descriptive names
- `trk_covarep_srh`: `SRH`

**Proposed standards**:
- `CPP[dB]` (Cepstral Peak Prominence)
- `HNRi[dB]` where i = frequency band (05, 15, 25)
  - `HNR05[dB]`, `HNR15[dB]`, `HNR25[dB]`
- `SHR[dB]` (Subharmonic-to-Harmonic Ratio)
- Jitter: `Jitter_local[%]`, `Jitter_RAP[%]`, `Jitter_PPQ5[%]`
- Shimmer: `Shimmer_local[dB]`, `Shimmer_APQ3[dB]`, `Shimmer_APQ5[dB]`

#### Category 7: Spectral Measures
**Current patterns**:
- `mfcc_1`, `mfcc_2`, ..., `mfcc_13` (Mel-Frequency Cepstral Coefficients)
- `ACF` (Autocorrelation Function)
- `RMS[dB]` (Root Mean Square)
- `Energy` (generic energy)

**Affected functions**:
- `mfcc`: `mfcc_1` through `mfcc_n`
- `trk_acfana`: `ACF`
- `trk_rmsana`: `RMS[dB]`

**Proposed standards**:
- `MFCCi` where i = 1, 2, 3, ... (uppercase, Titze-compatible)
  - Or keep `mfcc_i` for backwards compatibility
- `ACF` (dimensionless correlation)
- `RMS[dB]` (already correct)
- `Energy[dB]` (add unit)

#### Category 8: Other Measures
**Current patterns**:
- `voicing`, `voiced`, `vprob` (voicing probability)
- `periodicity`, `confidence` (pitch confidence)
- `corr` (correlation)
- `pm` (pitch marks)
- `VUV` (Voiced/Unvoiced/Mixed)

**Affected functions**:
- `trk_pyin`: `voiced`, `vprob`
- `trk_swiftf0`: `confidence`, `voicing`
- `trk_crepe`: `periodicity`
- `reaper_pm`: `pm`
- `trk_covarep_srh`: `VUV`

**Proposed standards**:
- `voicing_prob` (probability, 0-1)
- `voiced` (binary, 0/1)
- `periodicity` (0-1)
- `confidence` (0-1)
- `corr` (correlation, -1 to 1)
- `PM[s]` (pitch marks in seconds)
- `VUV` (categorical: 0/1/2)

## Migration Strategy

### Phase 1: Documentation & Mapping (No Code Changes)

**Goal**: Create comprehensive mapping tables without breaking existing code.

**Tasks**:
1. **Create master mapping table** (`TRACK_NAMES_MAPPING.csv`):
   ```csv
   function,old_name,new_name,unit,category,notes
   trk_rapt,f0,fo[Hz],Hz,f0,Titze 2015 standard
   trk_ksvfo,fo[Hz],fo[Hz],Hz,f0,Normalize to lowercase
   trk_swiftf0,F0,fo[Hz],Hz,f0,Normalize to lowercase
   trk_forest,F[Hz],Fi[Hz],Hz,formant,Multi-column formants
   trk_praat_sauce,H1H2c,H1-H2c[dB],dB,harmonic_diff,Add hyphen and unit
   ...
   ```

2. **Document all affected functions** with current→proposed mappings

3. **Create validation script** to verify all track names in codebase

4. **Impact analysis**: Identify downstream dependencies
   - emuR compatibility
   - User scripts that reference track names
   - Vignettes and examples

### Phase 2: Deprecation Period (Backwards Compatible)

**Goal**: Support both old and new names simultaneously for 1-2 major versions.

**Implementation**:
1. **Add alias system** in `AsspDataObj`:
   ```r
   # Internal function to map old names to new names
   .track_name_aliases <- function() {
     list(
       "fo" = "fo[Hz]",
       "F0" = "fo[Hz]",
       "fo" = "fo[Hz]",
       "fo[Hz]" = "fo[Hz]",
       "F1" = "F1[Hz]",
       "B1" = "B1[Hz]",
       "H1H2c" = "H1-H2c[dB]",
       # ... full mapping
     )
   }
   ```

2. **Update track access methods** to support both:
   ```r
   # When user accesses old name, return new name with deprecation warning
   `[[.AsspDataObj` <- function(x, name) {
     aliases <- .track_name_aliases()
     if (name %in% names(aliases)) {
       warning("Track name '", name, "' is deprecated. ",
               "Use '", aliases[[name]], "' instead.",
               call. = FALSE)
       name <- aliases[[name]]
     }
     NextMethod()
   }
   ```

3. **Update all function outputs** to use new names:
   ```r
   # Example: trk_rapt
   attr(trk_rapt, "tracks") <- c("fo[Hz]")  # Was: c("fo")
   ```

4. **Add `.track_names_version` attribute**:
   ```r
   attr(result, ".track_names_version") <- "2.0_titze2015"
   ```

### Phase 3: Full Migration (Breaking Change)

**Goal**: Remove old names, full compliance with standards.

**Timeline**: Version 1.0 or 2.0 (major version bump)

**Implementation**:
1. Remove all aliases and deprecation warnings
2. Update all documentation
3. Update all vignettes and examples
4. Release migration guide for users

## Detailed Mapping Tables

### Table 1: Fundamental Frequency (f0)

| Function | Old Name | New Name | Notes |
|----------|----------|----------|-------|
| `trk_rapt` | `fo` | `fo[Hz]` | Add unit |
| `trk_swipe` | `fo` | `fo[Hz]` | Add unit |
| `trk_dio` | `fo` | `fo[Hz]` | Add unit |
| `trk_reaper` | `fo` | `fo[Hz]` | Add unit |
| `trk_harvest` | `fo` | `fo[Hz]` | Add unit |
| `trk_yin` | `fo` | `fo[Hz]` | Add unit |
| `trk_pyin` | `fo` | `fo[Hz]` | Add unit |
| `trk_crepe` | `fo` | `fo[Hz]` | Add unit |
| `trk_yaapt` | `fo` | `fo[Hz]` | Add unit |
| `trk_torch_pitch` | `fo` | `fo[Hz]` | Add unit |
| `trk_kaldi_pitch` | `fo` | `fo[Hz]` | Add unit |
| `trk_snackp` | `fo` | `fo[Hz]` | Add unit |
| `trk_swiftf0` | `F0` | `fo[Hz]` | Normalize case + unit |
| `trk_ksvfo` | `fo[Hz]` | `fo[Hz]` | Add zero |
| `fo` | `fo[Hz]` | `fo[Hz]` | Add zero |
| `foana` | `fo[Hz]` | `fo[Hz]` | Add zero |
| `fo_ksv` | `fo[Hz]` | `fo[Hz]` | Add zero |
| `trk_mhspitch` | `fo[Hz]` | `fo[Hz]` | Add zero |
| `trk_covarep_srh` | `F0[Hz]` | `fo[Hz]` | Normalize case |

### Table 2: Formants

| Function | Old Name(s) | New Name(s) | Notes |
|----------|-------------|-------------|-------|
| `trk_praat_sauce` | `F1`, `F2`, `F3` | `F1[Hz]`, `F2[Hz]`, `F3[Hz]` | Add units |
| `trk_forest` | `F[Hz]` | `F1[Hz]`, `F2[Hz]`, `F3[Hz]`, `F4[Hz]` | Separate columns |
| `trk_formantp` | `fm` | `F1[Hz]`, `F2[Hz]`, `F3[Hz]`, `F4[Hz]` | Expand to specific formants |
| `trk_formantpathp` | `fm` | `F1[Hz]`, `F2[Hz]`, `F3[Hz]`, `F4[Hz]` | Expand to specific formants |
| `trk_snackf` | `F1`, `F2`, ... | `F1[Hz]`, `F2[Hz]`, ... | Add units |

### Table 3: Formant Bandwidths

| Function | Old Name(s) | New Name(s) | Notes |
|----------|-------------|-------------|-------|
| `trk_praat_sauce` | `B1`, `B2`, `B3` | `B1[Hz]`, `B2[Hz]`, `B3[Hz]` | Add units |
| `trk_forest` | `B[Hz]` | `B1[Hz]`, `B2[Hz]`, `B3[Hz]`, `B4[Hz]` | Separate columns |
| `trk_formantp` | `bw` | `B1[Hz]`, `B2[Hz]`, `B3[Hz]`, `B4[Hz]` | Expand to specific bandwidths |
| `trk_formantpathp` | `bw` | `B1[Hz]`, `B2[Hz]`, `B3[Hz]`, `B4[Hz]` | Expand to specific bandwidths |

### Table 4: Harmonics (VoiceSauce)

| Function | Old Name | New Name | Notes |
|----------|----------|----------|-------|
| `trk_praat_sauce` | `H1` | `H1[dB]` | Add unit |
| `trk_praat_sauce` | `H2` | `H2[dB]` | Add unit |
| `trk_praat_sauce` | `H4` | `H4[dB]` | Add unit |
| `trk_praat_sauce` | `A1` | `A1[dB]` | Add unit |
| `trk_praat_sauce` | `A2` | `A2[dB]` | Add unit |
| `trk_praat_sauce` | `A3` | `A3[dB]` | Add unit |
| `trk_praat_sauce` | `H2K` | `H2k[dB]` | Add unit, lowercase k |
| `trk_praat_sauce` | `H5K` | `H5k[dB]` | Add unit, lowercase k |

### Table 5: Harmonic Differences (VoiceSauce)

| Function | Old Name | New Name | Notes |
|----------|----------|----------|-------|
| `trk_praat_sauce` | `H1H2` | `H1-H2[dB]` | Add hyphen + unit |
| `trk_praat_sauce` | `H2H4` | `H2-H4[dB]` | Add hyphen + unit |
| `trk_praat_sauce` | `H1A1` | `H1-A1[dB]` | Add hyphen + unit |
| `trk_praat_sauce` | `H1A2` | `H1-A2[dB]` | Add hyphen + unit |
| `trk_praat_sauce` | `H1A3` | `H1-A3[dB]` | Add hyphen + unit |
| `trk_praat_sauce` | `H4H2K` | `H4-H2k[dB]` | Add hyphen + unit |
| `trk_praat_sauce` | `H2KH5K` | `H2k-H5k[dB]` | Add hyphen + unit |
| `trk_praat_sauce` | `H1H2c` | `H1-H2c[dB]` | Corrected version |
| `trk_praat_sauce` | `H1A1c` | `H1-A1c[dB]` | Corrected version |
| `trk_praat_sauce` | `H1A2c` | `H1-A2c[dB]` | Corrected version |
| `trk_praat_sauce` | `H1A3c` | `H1-A3c[dB]` | Corrected version |
| `trk_praat_sauce` | `H1H2u` | `H1-H2u[dB]` | Uncorrected version |
| `trk_praat_sauce` | `H1A1u` | `H1-A1u[dB]` | Uncorrected version |
| `trk_praat_sauce` | `H1A2u` | `H1-A2u[dB]` | Uncorrected version |
| `trk_praat_sauce` | `H1A3u` | `H1-A3u[dB]` | Uncorrected version |

### Table 6: Voice Quality Measures

| Function | Old Name | New Name | Notes |
|----------|----------|----------|-------|
| `trk_praat_sauce` | `CPP` | `CPP[dB]` | Add unit |
| `trk_praat_sauce` | `HNR05` | `HNR05[dB]` | Add unit |
| `trk_praat_sauce` | `HNR15` | `HNR15[dB]` | Add unit |
| `trk_praat_sauce` | `HNR25` | `HNR25[dB]` | Add unit |
| `trk_praat_sauce` | `HNR35` | `HNR35[dB]` | Add unit |
| `trk_praat_sauce` | `SHR` | `SHR[dB]` | Add unit |
| `trk_praat_sauce` | `Energy` | `Energy[dB]` | Add unit |
| `trk_covarep_srh` | `SRH` | `SRH[dB]` | Add unit |
| `lst_voice_reportp` | `Jitter (local)` | `Jitter_local[%]` | Remove parentheses, add unit |
| `lst_voice_reportp` | `Jitter (rap)` | `Jitter_RAP[%]` | Remove parentheses, add unit |
| `lst_voice_reportp` | `Shimmer (local)` | `Shimmer_local[dB]` | Remove parentheses, add unit |

### Table 7: Spectral Measures

| Function | Old Name | New Name | Notes |
|----------|----------|----------|-------|
| `mfcc` | `mfcc_1`, ... | `MFCC1`, `MFCC2`, ... | Uppercase (Titze) OR keep lowercase |
| `trk_acfana` | `ACF` | `ACF` | No change (dimensionless) |
| `trk_rmsana` | `RMS[dB]` | `RMS[dB]` | Already correct |
| `trk_zcrana` | `zcr` | `ZCR[Hz]` | Add unit |

### Table 8: Other Measures

| Function | Old Name | New Name | Notes |
|----------|----------|----------|-------|
| `trk_pyin` | `voiced` | `voiced` | Binary, no unit |
| `trk_pyin` | `vprob` | `voicing_prob` | Clearer name |
| `trk_swiftf0` | `confidence` | `confidence` | 0-1, no unit |
| `trk_swiftf0` | `voicing` | `voicing` | Binary |
| `trk_crepe` | `periodicity` | `periodicity` | 0-1, no unit |
| `reaper_pm` | `pm` | `PM[s]` | Pitch marks in seconds |
| `trk_covarep_srh` | `VUV` | `VUV` | Categorical |
| `trk_snackp` | `rms` | `RMS[dB]` | Normalize name + unit |

## Implementation Checklist

### Preparation Phase
- [ ] Read and analyze Titze 2015 paper in detail
- [ ] Read and analyze Nylén 2024 paper for extensions
- [ ] Create comprehensive inventory of all track names
- [ ] Create master mapping CSV file
- [ ] Identify all functions that use each track name
- [ ] Document VoiceSauce parameter equivalents
- [ ] Create validation script to check all track names

### Development Phase
- [ ] Implement `.track_name_aliases()` function
- [ ] Update `[[.AsspDataObj` with deprecation support
- [ ] Create helper function to convert old→new names
- [ ] Update all `attr(*, "tracks")` definitions
- [ ] Add `.track_names_version` attribute to outputs
- [ ] Update `as.data.frame.AsspDataObj` to handle aliases
- [ ] Update `as_tibble.AsspDataObj` to handle aliases

### Testing Phase
- [ ] Update all existing tests with new names
- [ ] Add backwards compatibility tests
- [ ] Test deprecation warnings
- [ ] Test unit conversion integration
- [ ] Verify emuR compatibility
- [ ] Test all example code

### Documentation Phase
- [ ] Update all function documentation
- [ ] Update vignettes with new names
- [ ] Create migration guide for users
- [ ] Update CHANGELOG with breaking changes
- [ ] Update README examples
- [ ] Create visual comparison table for old→new

### Release Phase
- [ ] Version bump (major for breaking changes)
- [ ] Announcement of deprecation period
- [ ] Blog post explaining rationale
- [ ] Update website documentation
- [ ] Notify emuR maintainers

## Benefits of Migration

1. **Scientific Compliance**: Align with JASA-published standards (Titze 2015, Nylén 2024)
2. **Cross-Software Compatibility**: Match notation used in Praat, VoiceSauce, phonetics literature
3. **Clarity**: Explicit units prevent confusion (Hz vs kHz, dB vs linear)
4. **Automatic Unit Conversion**: Integration with `units` package (already implemented)
5. **Gender/Voice Analysis**: Support for corrected formants (Nylén 2024 extensions)
6. **Consistency**: Unified notation across all 149 track definitions
7. **Citability**: Can reference JASA papers in documentation

## Risks & Mitigation

### Risk 1: Breaking User Code
**Mitigation**:
- Long deprecation period (1-2 versions)
- Alias system maintains backwards compatibility
- Clear migration guide
- Deprecation warnings guide users

### Risk 2: emuR Compatibility
**Mitigation**:
- Coordinate with emuR maintainers
- Ensure SSFF files remain compatible
- Test extensively

### Risk 3: Confusion During Transition
**Mitigation**:
- Clear version marking (`.track_names_version`)
- Documentation explains both systems
- Examples show new notation

### Risk 4: Performance Overhead
**Mitigation**:
- Alias lookup is O(1) hash table
- Only triggers on deprecated names
- Minimal impact

## Timeline Proposal

### Version 0.7.0 (Deprecation Start)
- Implement alias system
- Update all outputs to new names
- Deprecation warnings active
- Documentation shows both notations

### Version 0.8.0-0.9.0 (Transition Period)
- Refine based on user feedback
- Update vignettes progressively
- Encourage migration

### Version 1.0.0 (Clean Break)
- Remove aliases
- Only new names supported
- Full documentation update
- Migration complete

## References

1. Titze, I. R., et al. (2015). "Toward a consensus on symbolic notation of harmonics, resonances, and formants in vocalization." *The Journal of the Acoustical Society of America*, 137(5), 3005–3007. https://doi.org/10.1121/1.4919349

2. Nylén, F., et al. (2024). "Acoustic cues to femininity and masculinity in spontaneous speech." *The Journal of the Acoustical Society of America*, 155(2), 1373-1387. https://doi.org/10.1121/10.0024751

3. VoiceSauce Documentation: https://www.phonetics.ucla.edu/voicesauce/documentation/parameters.html

4. ISO/IEC 80000-3:2019 - Quantities and units (acoustics)

## Appendix: Complete Function Inventory

### F0 Tracking Functions (19)
- trk_rapt, trk_swipe, trk_dio, trk_reaper, trk_harvest
- trk_yin, trk_pyin, trk_crepe, trk_yaapt, trk_torch_pitch
- trk_kaldi_pitch, trk_snackp, trk_swiftf0, trk_ksvfo
- fo, foana, fo_ksv, trk_mhspitch, trk_covarep_srh

### Formant Tracking Functions (6)
- trk_forest, trk_formantp, trk_formantpathp, trk_snackf
- trk_praat_sauce (combined F0+formants+voice quality)

### Voice Quality Functions (5)
- trk_praat_sauce, lst_voice_reportp, lst_voice_tremorp
- lst_avqip, trk_covarep_srh

### Spectral Functions (8)
- mfcc, trk_acfana, trk_rmsana, trk_zcrana
- trk_cepstrum, trk_dftSpectrum, trk_lpsSpectrum, trk_cssSpectrum

### Pitch Mark Functions (2)
- reaper_pm, trk_estk_pm

### Other (10+)
- Various specialized DSP functions with unique track names

**Total**: ~149 track name definitions across ~50+ functions

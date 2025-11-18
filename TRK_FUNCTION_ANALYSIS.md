# Analysis of superassp trk_ Functions for Potential Migration to protoscribe

## Overview

Analyzed all 46 trk_/lst_ functions in superassp to identify candidates for migration to protoscribe based on output type: **EVENT annotations** (discrete time points) vs **DSP measurements** (regular intervals).

## Analysis Criteria

Functions should be migrated to protoscribe if they:
1. Output **discrete EVENT annotations** (sparse time points marking specific moments)
2. Mark **phonetic events** (boundaries, onsets, closures, etc.)
3. Do NOT output regular-interval measurements

Functions should STAY in superassp if they:
1. Output **continuous measurements** at regular time intervals
2. Provide DSP analysis (formants, spectrograms, MFCC, etc.)
3. Generate tracks with fixed windowShift

## Function Inventory

### Total Functions: 46
- **trk_ functions**: 34 (track/measurement functions)
- **lst_ functions**: 12 (list/summary functions)

## Analysis Results

### ✅ MIGRATED: Event Annotation Functions

#### 1. trk_pitchmark() → draft_pitchmark() ✓ COMPLETE
**Output**: Discrete glottal closure instants (EVENT)
**Rationale**: Marks specific phonetic moments (glottal closures)
**Status**: **MIGRATED** to protoscribe@17c0649
**Commit**: superassp@ba790f3 (deprecated)

---

### ⚠️ CANDIDATES: Functions with EVENT Components

#### 2. trk_reaper()
**Primary Output**: F0 track (continuous, SSFF)
**Secondary Output**: Epochs (discrete GCI events) as attribute

**Analysis**:
- Returns F0 track at regular intervals (DSP measurement)
- Also provides `attr(result, "epochs")` with GCI times
- **Primary purpose**: F0 estimation (DSP)
- **Secondary feature**: Epoch detection (EVENT)

**Recommendation**: **KEEP in superassp**
- Primary function is F0 tracking (regular intervals)
- Epochs are auxiliary data, not the main output
- If users need only epochs for annotation, use `draft_pitchmark()` or `draft_hgci()`

**Alternative**: Could add `draft_reaper_epochs()` in protoscribe if there's demand for REAPER-specific epoch extraction without F0 tracking.

---

### ✅ STAY in superassp: DSP Measurement Functions

All remaining functions output **continuous measurements at regular intervals**:

#### Pitch/F0 Tracking (15 functions)
Regular-interval F0 estimates (SSFF tracks):
- trk_crepe
- trk_dio
- trk_harvest
- trk_pitchp
- trk_pyin
- trk_rapt
- trk_reaper (primary output)
- trk_snackp
- trk_swiftf0
- trk_swipe
- trk_yaapt
- trk_yin
- trk_acfana (autocorrelation)
- trk_forest (formant estimation with F0)
- trk_ksvf0

**Rationale**: All output F0 values at fixed time intervals (windowShift parameter)

#### Formant Tracking (3 functions)
Regular-interval formant frequencies:
- trk_formantp
- trk_formantpathp (formant path tracking)
- trk_forest (combined F0/formant)

**Rationale**: Formants measured every N milliseconds, continuous DSP measurement

#### Spectral Analysis (7 functions)
Regular-interval spectral measurements:
- trk_cepstrum
- trk_cssSpectrum
- trk_lpsSpectrum
- trk_dftSpectrum
- trk_mfcc (mel-frequency cepstral coefficients)
- trk_spectral_momentsp
- trk_praat_sauce

**Rationale**: All generate spectral features at fixed time windows

#### Energy/Amplitude (3 functions)
Regular-interval amplitude measurements:
- trk_rmsana (RMS energy)
- trk_intensityp
- trk_zcrana (zero crossing rate)

**Rationale**: Energy measured at regular intervals

#### Voice Source Analysis (5 functions)
Regular-interval voice quality parameters:
- trk_covarep_iaif (glottal flow)
- trk_covarep_srh (spectral/cepstral harmonics)
- trk_aperiodicities (aperiodicity measures)
- trk_d4c (aperiodicity for WORLD vocoder)
- trk_excite (excitation signal)

**Rationale**: All continuous measurements of voice source characteristics

#### Spectro-Temporal (1 function)
- trk_seenc (spectro-temporal energy)

**Rationale**: Continuous spectro-temporal measurements

#### Utility (1 function)
- trk_npy_import (import numpy arrays as tracks)

**Rationale**: Generic import, preserves input format (typically regular-interval)

---

### ✅ STAY in superassp: List Functions

All lst_ functions output **summary statistics** or **feature sets**, not time-series data:

#### OpenSMILE Features (4 functions)
Summary statistics from entire recording:
- lst_ComParE_2016
- lst_eGeMAPS
- lst_emobase
- lst_GeMAPS

**Rationale**: Single-value summaries per recording, not annotations

#### Voice Quality Summaries (6 functions)
Clinical voice metrics:
- lst_avqip (Acoustic Voice Quality Index)
- lst_covarep_vq (COVAREP voice quality)
- lst_dsip (Dysphonia Severity Index)
- lst_vat (Voice Activity Time)
- lst_voice_reportp (voice report)
- lst_voice_tremorp (tremor analysis)

**Rationale**: Summary metrics, not time-based annotations

#### VoiceSauce (1 function)
Comprehensive voice metrics:
- lst_voice_sauce

**Rationale**: Feature extraction, not event annotation

#### Prosody (1 function)
Prosodic assessments:
- lst_dysprosody

**Rationale**: Clinical assessment scores, not temporal annotations

---

## Summary by Category

| Category | Keep in superassp | Migrate to protoscribe |
|----------|-------------------|------------------------|
| **Pitch/F0** | 15 functions | 0 |
| **Formants** | 3 functions | 0 |
| **Spectral** | 7 functions | 0 |
| **Energy** | 3 functions | 0 |
| **Voice Source** | 5 functions | 0 |
| **Other DSP** | 2 functions | 0 |
| **Summaries** | 12 functions | 0 |
| **Events** | 0 functions | 1 (pitchmark) ✓ |
| **TOTAL** | 47 functions | 1 function |

---

## Key Findings

### 1. trk_pitchmark() Was the Only EVENT Function
- **Only function** in superassp that output discrete EVENT annotations
- All other trk_ functions output continuous measurements at regular intervals
- Migration was **correct and necessary**

### 2. REAPER is Hybrid but Stays
- Primarily an F0 tracker (regular intervals)
- Epochs are auxiliary output (attribute)
- Primary use case is F0 estimation, not event annotation
- **Recommendation**: Stay in superassp

### 3. Clear Package Boundaries
The analysis confirms clear functional boundaries:

**superassp** = DSP measurements (regular intervals)
- Pitch tracking every N ms
- Formant tracking every N ms
- Spectral analysis every N ms
- Energy measurements every N ms
- Summary statistics

**protoscribe** = Event annotations (discrete points)
- Glottal closure instants (pitchmarks)
- VOT boundaries
- Pitch targets (MOMEL)
- Period markers
- Other phonetic events

### 4. No Additional Migrations Needed
The packages are now correctly organized. No other functions need migration.

---

## Potential Future Additions to Protoscribe

While no current superassp functions need migration, future EVENT annotation functions could be added to protoscribe:

### Voice Activity Detection (VAD)
**Function**: draft_vad()
**Output**: Start/end times of voiced segments
**Type**: SEGMENT boundaries (EVENT)

### Syllable Segmentation
**Function**: draft_syllables()
**Output**: Syllable boundary times
**Type**: BOUNDARY events

### Prosodic Boundary Detection
**Function**: draft_boundaries()
**Output**: Intonational phrase boundaries
**Type**: BOUNDARY events

### Burst Detection
**Function**: draft_bursts()
**Output**: Stop burst release times
**Type**: POINT events

These would be **NEW functions**, not migrations from superassp.

---

## Conclusion

### Migration Assessment: COMPLETE ✓

**Only 1 function needed migration**: `trk_pitchmark()` → `draft_pitchmark()`

**Status**: ✅ Successfully migrated (protoscribe@17c0649, superassp@ba790f3)

### Package Organization: CORRECT ✓

- **superassp**: 46 functions providing DSP measurements at regular intervals
- **protoscribe**: EVENT annotation generators (draft_ functions)

### Recommendation: NO FURTHER MIGRATIONS NEEDED

The packages are correctly organized:
- superassp handles continuous DSP analysis
- protoscribe handles discrete event annotation
- Clear functional boundaries
- No overlap or misplaced functions

---

## REAPER Special Case

### Question: Should REAPER epochs be extracted separately?

**Current**: trk_reaper() returns F0 track + epochs as attribute
**Alternative**: Add draft_reaper_epochs() to protoscribe?

**Analysis**:
- REAPER's primary purpose is F0 estimation
- Epochs are auxiliary data for voice source analysis
- Most users want F0, not just epochs
- Users needing only GCI events can use:
  - `draft_pitchmark()` (ESTK algorithm, optimized for EGG)
  - `draft_hgci()` (Hilbert transform, works with acoustic)

**Recommendation**: **NO** - Keep REAPER in superassp
- Primary function is F0 tracking (DSP)
- Epoch extraction is secondary feature
- Better dedicated options exist for GCI annotation

If demand emerges for REAPER-specific epoch extraction without F0, could add later as a new function in protoscribe.

---

## References

**Edinburgh Speech Tools pitchmark**: 
- Macon, M. W., & Taylor, P. (1997)
- https://www.cstr.ed.ac.uk/projects/speech_tools/

**REAPER (Robust Epoch And Pitch EstimatoR)**:
- David Talkin, Google
- https://github.com/google/REAPER

**Migration commits**:
- protoscribe@17c0649 (implementation)
- superassp@ba790f3 (deprecation)

---

## Date: 2025-10-25

**Analyst**: Claude (Anthropic)
**Review Status**: Complete
**Migration Status**: 1 of 1 candidates migrated (100%)

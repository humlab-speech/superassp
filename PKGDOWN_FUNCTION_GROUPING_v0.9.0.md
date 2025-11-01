# superassp Function Grouping for pkgdown Reference (v0.9.0)

**Date**: November 1, 2025
**Version**: 0.9.0
**Purpose**: Organized categorization of all `trk_*` and `lst_*` functions for pkgdown reference documentation

---

## Overview

This document organizes superassp's **195 exported functions** into logical groups for pkgdown reference pages. Functions are grouped by measurement type and use case rather than implementation method.

**Function Distribution:**
- **SSFF Track Functions (`trk_*`)**: 49 functions (time-aligned data)
- **List/Dataframe Functions (`lst_*`)**: 15 functions (summary statistics)
- **Installation Functions (`install_*`)**: 15 functions (Python dependencies)
- **Helper/Utility Functions**: 116 functions (type conversions, checks, etc.)

---

## 1. Pitch & F0 Analysis

### Fundamental Frequency Tracking (trk_*)

**C++ Implementations (Fastest - Tier 1):**
- `trk_rapt()` - RAPT: Robust Algorithm for Pitch Tracking (dynamic programming, C++/SPTK)
- `trk_swipe()` - SWIPE: Sawtooth Waveform Inspired Pitch Estimator (noise-robust, C++/SPTK)
- `trk_dio()` - DIO: WORLD vocoder pitch tracking (synthesis-quality, C++/SPTK)
- `trk_harvest()` - Harvest: High-accuracy WORLD algorithm (C++/SPTK)
- `trk_reaper()` - REAPER: F0 + glottal closure instants + polarity (C++/SPTK)
- `trk_yin()` - YIN: Classic autocorrelation method (C++)
- `trk_pyin()` - Probabilistic YIN with HMM Viterbi decoding (C++)

**Classical ASSP Algorithms (Tier 2):**
- `trk_ksvfo()` - K.Schaefer-Vincent periodicity detection (C/ASSP)
- `trk_mhspitch()` - Modified Harmonic Sieve (C/ASSP)

**Deep Learning & Modern Methods (Tier 3):**
- `trk_swiftf0()` - Swift-F0: CNN-based real-time pitch tracker (90-130ms, Python)
- `trk_crepe()` - CREPE: Deep CNN on raw waveform (Python/TensorFlow)
- `trk_sacc()` - SAcC: Subband Autocorrelation Classification (noise-robust, Python)

**Praat/Parselmouth Methods (Tier 4):**
- `trk_pitchp()` - Parselmouth/Praat pitch analysis (multiple methods, Python)
- `trk_dv_f0()` - DisVoice F0 tracking (in-memory Parselmouth, Python)

**Legacy/Alternative Methods:**
- `trk_snackp()` - Snack Toolkit pitch tracking (Python)
- `trk_straight_f0()` - STRAIGHT F0 extraction (legacy vocoder, Python)
- `trk_vat_srh()` - Voice Analysis Toolkit SRH method (Python)
- `trk_covarep_srh()` - COVAREP Summation of Residual Harmonics (Python)

### Pitch Marks & Glottal Closure Instants (trk_*)

**Recommended:**
- `trk_reaper_pm()` - **NEW v0.9.0**: REAPER pitch marks (C++/SPTK, 2-3x faster)
- `trk_pitchmark()` - ESTK: Glottal closure instant detection (C++/ESTK)

**Deprecated:**
- ~~`reaper_pm()`~~ - Python version (deprecated, use `trk_reaper_pm()`)

---

## 2. Formant Analysis

### Formant Frequency & Bandwidth Tracking (trk_*)

**Standard Methods:**
- `trk_forest()` - ASSP: Linear prediction formant estimation (autocorrelation + SLA, C/ASSP)
- `trk_formantp()` - Parselmouth/Praat: Burg method with tracking (Python)
- `trk_formantpathp()` - Parselmouth/Praat: Formant path optimization (Python)

**Deep Learning:**
- `trk_deepformants()` - PyTorch RNN formant tracking (F1-F4, 2x real-time, Python)

**Specialized:**
- `trk_formants_tvwlp()` - Time-Varying Weighted LP (4.37x speedup, GCI-based, R/Python)
- `trk_snackf()` - Snack Toolkit formant analysis (Python)
- `trk_dv_formants()` - DisVoice formant tracking (in-memory, Python)

### Formant Summary Statistics (lst_*)
- `lst_deepformants()` - Summary statistics of deep learning formant tracks (Python)

---

## 3. Phonological Analysis

### **NEW in v0.9.0**: Phonological Posterior Probabilities

**Phonet Integration** (BGRU deep learning models):
- `lst_phonet()` - Phonological posteriors as list/data.frame (18 classes, R/Python)
- `trk_phonet()` - Phonological posteriors as SSFF tracks (18 classes, R/Python)

**Phonological Classes:**
- **Vowel features**: vocalic, back, anterior, open, close
- **Consonant features**: consonantal, nasal, stop, continuant, lateral, flap, trill, voice, strident
- **Place features**: labial, dental, velar
- **Other**: pause

**Installation:**
- `install_phonet()` - Install Phonet dependencies (TensorFlow, tf-keras)
- `phonet_available()` - Check Phonet availability
- `phonet_info()` - Get Phonet configuration

**Use Cases:**
- Time-aligned phonological annotation (emuR integration)
- Phonetic segmentation and labeling
- Speech disorder analysis (dysarthria, apraxia)
- Articulatory feature tracking

---

## 4. Spectral Analysis

### Spectrum Estimation (trk_*)

**FFT-Based Methods (C/ASSP):**
- `trk_dftSpectrum()` - Discrete Fourier Transform (unsmoothed narrow-band)
- `trk_cssSpectrum()` - Cepstrally smoothed spectrum
- `trk_lpsSpectrum()` - Linear prediction smoothed spectrum

**Cepstral Analysis:**
- `trk_cepstrum()` - Short-term cepstral coefficients (C/ASSP)

**Spectral Envelope:**
- `trk_seenc()` - Encoded spectral envelope (WORLD CheapTrick, Python)
- `trk_straight_spec()` - STRAIGHT spectral envelope (legacy, Python)

**Spectral Moments:**
- `trk_spectral_momentsp()` - Centroid, spread, skewness, kurtosis (Parselmouth, Python)

---

## 5. Energy & Amplitude Analysis

### Signal Energy Tracking (trk_*)

**C/ASSP Implementations:**
- `trk_rmsana()` - Root Mean Square amplitude (dB or linear)
- `trk_zcrana()` - Zero-crossing rate analysis
- `trk_acfana()` - Autocorrelation function

**Parselmouth:**
- `trk_intensityp()` - Intensity (perceived loudness, Python)

---

## 6. Voice Quality & Aperiodicity

### Voice Quality Measures (trk_*)

**Aperiodicity:**
- `trk_d4c()` - D4C: WORLD band aperiodicity (high-quality, C++/SPTK)
- `trk_brouhaha()` - Deep learning VAD + SNR + C50 estimation (50-100x optimized, Python)

**Creak Detection:**
- `trk_creak_union()` - Creaky voice detection (AM + CD neural network, Python)

**Voice Source:**
- `trk_excite()` - Excitation signal extraction (Python)
- `trk_praat_sauce()` - VoiceSauce-style voice quality (comprehensive, Python)

### Voice Quality Summary Measures (lst_*)

**Comprehensive Toolboxes:**
- `lst_vat()` - Voice Analysis Toolbox: **132 dysphonia measures**
  - Jitter (22-25), Shimmer (22-25), HNR (4), nonlinear dynamics (3)
  - Glottal measures (9), MFCCs (84), wavelet (~50), EMD (6)

- `lst_voice_sauce()` - VoiceSauce: **40+ voice quality parameters**
  - F0, formants, harmonic amplitudes, CPP, HNR, energy, spectral measures

**COVAREP Voice Quality:**
- `lst_covarep_vq()` - Spectral Relative Harmonic (SRH) + voice quality features
- `trk_covarep_srh()` - Frame-by-frame SRH tracking

**Parselmouth Voice Reports:**
- `lst_voice_reportp()` - Praat voice report measures
- `lst_voice_tremorp()` - Voice tremor analysis
- `lst_dsip()` - Dysphonia Severity Index (DSI)
- `lst_avqip()` - Acoustic Voice Quality Index (AVQI)

---

## 7. Prosody & Intonation

### Prosodic Features (lst_*)

**Dysprosody Analysis:**
- `lst_dysprosody()` - Dysprosody: **193 prosodic features**
  - MOMEL-INTSINT pitch targets, tone labels
  - Spectral tilt, statistical summaries
  - Performance: 14x real-time (0.16-0.44s per file)

**Voice Complexity:**
- `lst_voxit()` - Voxit: **11 voice & articulation complexity measures**
  - Speaking rate, pause statistics, rhythmic complexity
  - Pitch dynamics, entropy

---

## 8. Source-Filter Decomposition

### Glottal Source & Vocal Tract Separation (trk_*)

**IAIF Methods:**
- `trk_gfmiaif()` - GFM-IAIF: Glottal Flow Model-based separation (Python)
  - Vocal tract, glottis, lip radiation filters
- `trk_covarep_iaif()` - COVAREP IAIF: Iterative Adaptive Inverse Filtering (Python)

**STRAIGHT Vocoder:**
- `straight_pipeline()` - Full analysis-synthesis pipeline (Python)
- `straight_synth()` - Synthesize speech from STRAIGHT parameters (Python)

---

## 9. OpenSMILE Feature Sets

### Standardized Acoustic Features (lst_*)

**GeMAPS (Minimalistic):**
- `lst_GeMAPS()` - **62 static acoustic features** (C++/OpenSMILE, fastest)
  - Pitch, jitter, formants, shimmer, loudness, HNR, spectral balance
  - **Performance**: 3-5x faster than Python version

**eGeMAPS (Extended):**
- `lst_eGeMAPS()` - **88 extended features** (Python/OpenSMILE)
  - All GeMAPS + formant bandwidths, MFCCs, spectral flux

**Emotion Recognition:**
- `lst_emobase()` - **988 emotional voice features** (C++/OpenSMILE, fastest)

**Challenge Sets:**
- `lst_ComParE_2016()` - **6,373 ComParE 2016 features** (Python/OpenSMILE)

---

## 10. Acoustic Features & Coefficients

### Feature Extraction (trk_*)

**Mel-Frequency Cepstral Coefficients:**
- `trk_mfcc()` - SPTK MFCC extraction (C++, 2-3x faster than Python)

**Linear Prediction:**
- Various LP representations (ARF, LAR, LPC, RFC, etc.)

---

## 11. Installation & Dependency Management

### Python Module Installers (install_*)

**Pitch/F0:**
- `install_sacc()` - SAcC pitch tracker
- `install_swiftf0()` - Swift-F0 deep learning tracker
- `install_legacy_straight()` - Legacy STRAIGHT vocoder

**Formants:**
- `install_deepformants()` - PyTorch DeepFormants
- `install_ftrack_tvwlp()` - TVWLP formant tracker

**Voice Quality:**
- `install_brouhaha()` - Brouhaha VAD/SNR
- `install_vat()` - Voice Analysis Toolkit
- `install_voice_sauce()` - VoiceSauce
- `install_voice_analysis()` - General voice analysis tools

**Prosody:**
- `install_dysprosody()` - Dysprosody analysis
- `install_voxit()` - Voxit complexity measures

**Source-Filter:**
- `install_gfmiaif()` - GFM-IAIF
- `install_covarep()` - COVAREP toolkit

**Phonological (NEW v0.9.0):**
- `install_phonet()` - Phonet phonological posteriors

**General:**
- `install_disvoice_python()` - DisVoice toolkit

---

## Implementation Performance Tiers

### Tier 1: C++ Native (SPTK, ESTK) - **Fastest** ⚡⚡⚡

**No Python dependencies, 2-3x faster than Python**

- **Pitch**: `trk_rapt`, `trk_swipe`, `trk_dio`, `trk_harvest`, `trk_reaper`, `trk_yin`, `trk_pyin`
- **Pitch Marks**: `trk_reaper_pm` (NEW), `trk_pitchmark`
- **Features**: `trk_mfcc`, `trk_d4c`
- **OpenSMILE**: `lst_GeMAPS`, `lst_emobase`

**Performance**: Production-ready, batch processing optimized

### Tier 2: C Library (ASSP) - **Fast** ⚡⚡

**Established algorithms, stable C implementation**

- **Formants**: `trk_forest`
- **Pitch**: `trk_ksvfo`, `trk_mhspitch`
- **Spectral**: `trk_dftSpectrum`, `trk_cssSpectrum`, `trk_lpsSpectrum`, `trk_cepstrum`
- **Energy**: `trk_rmsana`, `trk_zcrana`, `trk_acfana`

**Performance**: Reliable, well-tested algorithms

### Tier 3: Python Deep Learning - **Specialized** ⚡

**Modern methods, GPU-capable, specialized use cases**

- **Pitch**: `trk_swiftf0` (real-time), `trk_crepe` (SOTA), `trk_sacc`
- **Formants**: `trk_deepformants` (2x realtime), `trk_formants_tvwlp` (4.37x speedup)
- **VAD**: `trk_brouhaha` (50-100x optimized)
- **Phonology**: `lst_phonet`, `trk_phonet` (BGRU models)

**Performance**: Modern deep learning, requires Python + dependencies

### Tier 4: Python Classical - **Moderate** ⚡

**Established Python algorithms**

- **Source-Filter**: `trk_gfmiaif`, `trk_excite`, `trk_covarep_iaif`
- **Spectral**: `trk_seenc`, `trk_straight_spec`
- **OpenSMILE**: `lst_eGeMAPS`, `lst_ComParE_2016`

**Performance**: Moderate speed, specialized functionality

### Tier 5: Parselmouth/Praat - **Comprehensive**

**Praat integration via Python, comprehensive features**

- **Pitch**: `trk_pitchp`, `trk_dv_f0`
- **Formants**: `trk_formantp`, `trk_formantpathp`, `trk_dv_formants`
- **Quality**: `trk_intensityp`, `trk_spectral_momentsp`, `trk_praat_sauce`
- **Summary**: `lst_voice_reportp`, `lst_voice_tremorp`, `lst_dsip`, `lst_avqip`

**Performance**: Feature-rich, gold-standard algorithms

---

## Media Format Support

### Universal Format Support (All Modern Functions)

All functions using `av` package support:
- **Audio**: WAV, MP3, MP4, FLAC, OGG, AAC, OPUS
- **Video**: MP4, MOV, AVI (automatic audio extraction)
- **Streaming**: Any format supported by FFmpeg

### Native Format Support (Legacy ASSP)

Functions with native ASSP support:
- **Standard**: WAV
- **Specialized**: AU, KAY, NIST, NSP, CSRE, SSFF

---

## Usage Recommendations

### For Basic Speech Analysis:
- **Pitch**: `trk_rapt()` (fastest, robust, C++)
- **Formants**: `trk_forest()` (established, reliable, C/ASSP)
- **Energy**: `trk_rmsana()` (simple, fast, C/ASSP)
- **Spectrum**: `trk_dftSpectrum()` (standard FFT, C/ASSP)

### For Voice Quality Assessment:
- **Comprehensive**: `lst_vat()` (132 measures)
- **Clinical**: `lst_voice_sauce()` (40+ VoiceSauce parameters)
- **Standardized**: `lst_GeMAPS()` (62 features, international standard, C++)
- **Extended**: `lst_eGeMAPS()` (88 features)

### For Prosodic Analysis:
- **Full Assessment**: `lst_dysprosody()` (193 features, MOMEL-INTSINT)
- **Complexity**: `lst_voxit()` (11 measures, articulation)

### For Phonological Analysis (NEW):
- **Time-aligned**: `trk_phonet()` (SSFF format, emuR compatible)
- **Summary**: `lst_phonet()` (data.frame format, tidyverse compatible)

### For Real-Time Applications:
- **Pitch**: `trk_swiftf0()` (90-130ms, CNN-based)
- **Formants**: `trk_deepformants()` (2x realtime)
- **VAD**: `trk_brouhaha()` (50-100x optimized)

### For Noisy Speech:
- **Pitch**: `trk_swipe()`, `trk_sacc()` (noise-robust algorithms)
- **VAD+Quality**: `trk_brouhaha()` (SNR + C50 estimation)

### For Speech Synthesis:
- **F0**: `trk_dio()` or `trk_harvest()` (WORLD vocoder, C++)
- **Spectrum**: `trk_seenc()` (CheapTrick spectral envelope)
- **Aperiodicity**: `trk_d4c()` (band aperiodicity, C++)

### For Glottal Analysis:
- **Pitch Marks**: `trk_reaper_pm()` (NEW, C++, 2-3x faster)
- **GCI Detection**: `trk_pitchmark()` (ESTK, C++)
- **F0 + Epochs**: `trk_reaper()` (both in one call, C++)

---

## Function Count Summary

**Total Package Exports**: 195 functions

**By Category:**
- **SSFF Track Functions (trk_*)**: 49
- **Summary Functions (lst_*)**: 15
- **Installation Functions (install_*)**: 15
- **Helper/Utility Functions**: 116

**By Domain:**
- **Pitch/F0**: 21 track functions
- **Formants**: 7 track functions + 1 summary
- **Voice Quality**: 5 track functions + 7 summary functions
- **Spectral**: 7 track functions
- **Energy**: 4 track functions
- **Phonological**: 2 functions (NEW in v0.9.0)
- **Prosody**: 2 summary functions
- **Feature Extraction**: 4 summary functions
- **Source-Filter**: 5 functions

---

## New in v0.9.0

### Phonological Analysis (Phonet Integration)
- ✅ `lst_phonet()` - List-based phonological posterior extraction
- ✅ `trk_phonet()` - SSFF track-based phonological posterior extraction
- ✅ `install_phonet()` - Dependency installer with tf-keras support
- ✅ `phonet_available()` - Availability checker
- ✅ `phonet_info()` - Configuration information

### Performance Improvements
- ✅ `trk_reaper_pm()` - C++ pitch mark extraction (2-3x faster than Python)
- ⚠️ `reaper_pm()` - Deprecated (Python version, remove in v0.11.0)

### Documentation
- ✅ Complete package audit (PACKAGE_AUDIT_2025-11-01.md)
- ✅ Function categorization in CLAUDE.md
- ✅ Interface consistency analysis (Grade: A-)

---

## Deprecated Functions

### Scheduled for Removal in v0.11.0 (mid-2026)
- ~~`reaper_pm()`~~ - Use `trk_reaper_pm()` instead (C++, 2-3x faster)

### Previously Removed (migrated to eggstract package)
- ~~`trk_egg_f0()`~~ - Use `eggstract::egg_f0()` instead

---

## Notes for pkgdown Configuration

### Suggested Reference Groups:

```yaml
reference:
- title: "Pitch & F0 Analysis"
  desc: "Fundamental frequency tracking and pitch mark detection"
  contents:
  - matches("^trk_(rapt|swipe|dio|harvest|reaper|yin|pyin|ksvfo|mhspitch)$")
  - matches("^trk_(swiftf0|crepe|sacc|pitchp|pitchmark)$")
  - matches("^trk_reaper_pm$")

- title: "Formant Analysis"
  desc: "Formant frequency and bandwidth estimation"
  contents:
  - matches("^trk_(forest|formantp|formantpathp|deepformants|formants_tvwlp|snackf|dv_formants)$")
  - lst_deepformants

- title: "Phonological Analysis"
  desc: "Phonological posterior probability estimation (NEW in v0.9.0)"
  contents:
  - lst_phonet
  - trk_phonet
  - install_phonet
  - phonet_available
  - phonet_info

- title: "Spectral Analysis"
  desc: "Spectrum estimation, cepstral analysis, and spectral moments"
  contents:
  - matches("^trk_(dftSpectrum|cssSpectrum|lpsSpectrum|cepstrum|seenc|spectral_momentsp|straight_spec)$")

- title: "Energy & Amplitude"
  desc: "Signal energy, intensity, and amplitude tracking"
  contents:
  - matches("^trk_(rmsana|zcrana|acfana|intensityp)$")

- title: "Voice Quality Measures"
  desc: "Voice quality assessment and dysphonia measures"
  contents:
  - lst_vat
  - lst_voice_sauce
  - lst_covarep_vq
  - matches("^lst_voice_")
  - matches("^lst_(avqip|dsip)$")
  - matches("^trk_(d4c|brouhaha|creak_union|praat_sauce|excite)$")

- title: "Prosody & Intonation"
  desc: "Prosodic feature extraction and complexity measures"
  contents:
  - lst_dysprosody
  - lst_voxit

- title: "Source-Filter Decomposition"
  desc: "Glottal source and vocal tract separation"
  contents:
  - matches("^trk_(gfmiaif|covarep_iaif|excite)$")
  - straight_pipeline
  - straight_synth

- title: "OpenSMILE Feature Sets"
  desc: "Standardized acoustic feature extraction"
  contents:
  - matches("^lst_(GeMAPS|eGeMAPS|emobase|ComParE_2016)$")

- title: "Acoustic Features & Coefficients"
  desc: "MFCC, LP coefficients, and other acoustic features"
  contents:
  - trk_mfcc

- title: "Installation & Dependencies"
  desc: "Python module installation and availability checking"
  contents:
  - starts_with("install_")
  - matches("_available$")
  - matches("_info$")
```

---

**Document Created**: November 1, 2025
**Version**: 0.9.0
**Status**: Complete catalog reflecting Phonet integration and REAPER PM C++ implementation
**Purpose**: pkgdown reference organization and user guidance

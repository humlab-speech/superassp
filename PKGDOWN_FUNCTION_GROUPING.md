# superassp Function Grouping for pkgdown References

**Date**: 2025-10-28
**Purpose**: Organized categorization of all `trk_*` and `lst_*` functions for pkgdown reference documentation

---

## Overview

This document organizes superassp's 75+ DSP functions into logical groups for pkgdown reference pages. Functions are grouped by measurement type and use case rather than implementation method.

---

## 1. Pitch & F0 Analysis

### Fundamental Frequency Tracking (trk_*)

**Recommended C++ Implementations (Fastest - 2-3x faster than Python):**
- `trk_rapt()` - RAPT: Robust Algorithm for Pitch Tracking (dynamic programming)
- `trk_swipe()` - SWIPE: Sawtooth Waveform Inspired Pitch Estimator (noise-robust)
- `trk_dio()` - DIO: WORLD vocoder pitch tracking (synthesis-quality)
- `trk_harvest()` - Harvest: High-accuracy WORLD algorithm
- `trk_reaper()` - REAPER: F0 + glottal closure instants

**Classical Algorithms:**
- `trk_ksvfo()` - K.Schaefer-Vincent periodicity detection (ASSP C library)
- `trk_mhspitch()` - Modified Harmonic Sieve (ASSP C library)

**Deep Learning & Modern Methods:**
- `trk_swiftf0()` - Swift-F0: CNN-based real-time pitch tracker (90-130ms)
- `trk_crepe()` - CREPE: Deep CNN on raw waveform
- `trk_sacc()` - SAcC: Subband Autocorrelation Classification (noise-robust)

**Alternative Methods:**
- `trk_pyin()` - Probabilistic YIN with HMM Viterbi decoding
- `trk_yin()` - Classic YIN autocorrelation method
- `trk_yaapt()` - Yet Another Algorithm for Pitch Tracking (NCC + DP)
- `trk_snackp()` - Snack Toolkit pitch tracking
- `trk_pitchp()` - Parselmouth/Praat pitch analysis (multiple methods)
- `trk_dv_f0()` - DisVoice F0 tracking (in-memory Parselmouth)

**Specialized:**
- `trk_pitchmark()` - ESTK: Glottal closure instants (laryngograph)

---

## 2. Formant Analysis

### Formant Frequency & Bandwidth Tracking (trk_*)

**Standard Methods:**
- `trk_forest()` - ASSP: Linear prediction formant estimation (autocorrelation + SLA)
- `trk_formantp()` - Parselmouth/Praat: Burg method with tracking
- `trk_formantpathp()` - Parselmouth/Praat: Formant path tracking
- `trk_dv_formants()` - DisVoice formant tracking (in-memory)

**Advanced Methods:**
- `trk_deepformants()` - PyTorch RNN formant tracking (2x real-time)
- `trk_formants_tvwlp()` - Time-Varying Weighted LP (4.37x speedup, GCI-based)
- `trk_snackf()` - Snack Toolkit formant analysis

**Summary Statistics:**
- `lst_deepformants()` - Summary statistics of deep learning formant tracks

---

## 3. Spectral Analysis

### Spectrum Estimation (trk_*)

**FFT-Based Methods:**
- `trk_dftSpectrum()` - Discrete Fourier Transform (unsmoothed narrow-band)
- `trk_cssSpectrum()` - Cepstrally smoothed spectrum
- `trk_lpsSpectrum()` - Linear prediction smoothed spectrum

**Cepstral Analysis:**
- `trk_cepstrum()` - Short-term cepstral coefficients

**Spectral Envelope:**
- `trk_seenc()` - Encoded spectral envelope (WORLD CheapTrick)

**Spectral Moments:**
- `trk_spectral_momentsp()` - Centroid, spread, skewness, kurtosis

---

## 4. Energy & Amplitude Analysis

### Signal Energy Tracking (trk_*)

- `trk_rmsana()` - Root Mean Square amplitude (dB or linear)
- `trk_zcrana()` - Zero-crossing rate analysis
- `trk_acfana()` - Autocorrelation function
- `trk_intensityp()` - Parselmouth intensity (perceived loudness)

---

## 5. Voice Quality & Aperiodicity

### Voice Quality Measures (trk_*)

**Aperiodicity:**
- `trk_d4c()` - D4C: WORLD band aperiodicity (high-quality)
- `trk_brouhaha()` - Deep learning VAD + SNR + C50 estimation (50-100x optimized)

**Creak Detection:**
- `trk_creak_union()` - Creaky voice detection (AM + CD neural network)

**Voice Source:**
- `trk_excite()` - Excitation signal extraction
- `trk_praat_sauce()` - VoiceSauce-style voice quality

### Voice Quality Summary Measures (lst_*)

**Comprehensive Toolboxes:**
- `lst_vat()` - Voice Analysis Toolbox: 132 dysphonia measures
  - Jitter (22-25), Shimmer (22-25), HNR (4), nonlinear dynamics (3)
  - Glottal measures (9), MFCCs (84), wavelet (~50), EMD (6)

- `lst_voice_sauce()` - VoiceSauce: 40+ voice quality parameters
  - F0, formants, harmonic amplitudes, CPP, HNR, energy, spectral measures

**COVAREP Voice Quality:**
- `lst_covarep_vq()` - Spectral Relative Harmonic (SRH) + voice quality
- `trk_covarep_srh()` - Frame-by-frame SRH tracking

**Parselmouth Voice Reports:**
- `lst_voice_reportp()` - Praat voice report measures
- `lst_voice_tremorp()` - Voice tremor analysis
- `lst_dsip()` - Dysphonia Severity Index
- `lst_avqip()` - Acoustic Voice Quality Index

---

## 6. Prosody & Intonation

### Prosodic Features (lst_*)

- `lst_dysprosody()` - Dysprosody: 193 prosodic features
  - MOMEL-INTSINT pitch targets, tone labels
  - Spectral tilt, statistical summaries
  - Performance: 14x real-time (0.16-0.44s per file)

- `lst_voxit()` - Voxit: 11 voice & articulation complexity measures
  - Speaking rate, pause statistics, rhythmic complexity
  - Pitch dynamics, entropy

---

## 7. Source-Filter Decomposition

### Glottal Source & Vocal Tract Separation (trk_*)

- `trk_gfmiaif()` - GFM-IAIF: Glottal Flow Model-based separation
  - Vocal tract, glottis, lip radiation filters
- `trk_covarep_iaif()` - COVAREP IAIF: Iterative Adaptive Inverse Filtering

---

## 8. OpenSMILE Feature Sets

### Standardized Acoustic Features (lst_*)

**GeMAPS (Minimalistic):**
- `lst_GeMAPS()` / `lst_GeMAPS_cpp()` / `lst_GeMAPS_python()` - 62 static acoustic features
  - Pitch, jitter, formants, shimmer, loudness, HNR, spectral balance
  - C++ version: 3-5x faster than Python

**eGeMAPS (Extended):**
- `lst_eGeMAPS()` / `lst_eGeMAPS_cpp()` / `lst_eGeMAPS_python()` - 88 extended features
  - All GeMAPS + formant bandwidths, MFCCs, spectral flux

**Emotion & Challenge Sets:**
- `lst_emobase()` / `lst_emobase_cpp()` / `lst_emobase_python()` - Emotional voice features
- `lst_ComParE_2016()` / `lst_ComParE_2016_cpp()` / `lst_ComParE_2016_python()` - Challenge features

**Custom:**
- `lst_GenericOpenSmile()` - Configurable OpenSMILE feature extraction

---

## 9. Acoustic Features & Coefficients

### Feature Extraction (trk_*)

**Mel-Frequency Cepstral Coefficients:**
- `trk_mfcc()` - SPTK MFCC extraction (C++, 2-3x faster than Python)

**Linear Prediction:**
- `trk_lp_analysis()` - LP coefficient representations (ARF, LAR, LPC, RFC)

**Utilities:**
- `trk_npy_import()` - Import Python-generated numpy arrays

---

## Implementation Performance Tiers

### By Speed (Fastest to Slowest)

**Tier 1: C++ Native (SPTK, ESTK)** - Fastest, no Python dependencies
- Pitch: `trk_rapt()`, `trk_swipe()`, `trk_dio()`, `trk_harvest()`, `trk_reaper()`
- Features: `trk_mfcc()`, `trk_d4c()`
- Epochs: `trk_pitchmark()`

**Tier 2: C Library (ASSP)** - Medium speed, established algorithms
- Formants: `trk_forest()`
- Pitch: `trk_ksvfo()`, `trk_mhspitch()`
- Spectral: `trk_dftSpectrum()`, `trk_cssSpectrum()`, `trk_lpsSpectrum()`, `trk_cepstrum()`
- Energy: `trk_rmsana()`, `trk_zcrana()`, `trk_acfana()`
- LP: `trk_lp_analysis()`

**Tier 3: Python Deep Learning** - Modern methods, GPU-capable
- Pitch: `trk_swiftf0()` (real-time), `trk_crepe()`, `trk_sacc()`
- Formants: `trk_deepformants()` (2x realtime), `trk_formants_tvwlp()` (4.37x speedup)
- VAD: `trk_brouhaha()` (50-100x optimized)

**Tier 4: Python Classical** - Established algorithms, moderate speed
- Pitch: `trk_yin()`, `trk_pyin()`, `trk_yaapt()`
- Source: `trk_gfmiaif()`, `trk_excite()`
- Spectral: `trk_seenc()`

**Tier 5: Parselmouth/Praat** - Comprehensive features, recently optimized
- Pitch: `trk_pitchp()`, `trk_dv_f0()`
- Formants: `trk_formantp()`, `trk_formantpathp()`, `trk_dv_formants()`
- Quality: `trk_intensityp()`, `trk_spectral_momentsp()`, `trk_praat_sauce()`
- Summary: `lst_voice_reportp()`, `lst_voice_tremorp()`, `lst_dsip()`, `lst_avqip()`

---

## Media Format Support

### Universal Format Support (All Modern Functions)
All functions using `av` package support:
- Audio: WAV, MP3, MP4, FLAC, OGG, AAC, OPUS
- Video: MP4, MOV, AVI (automatic audio extraction)
- Streaming: Any format supported by FFmpeg

### Native Format Support (Legacy)
Functions with native ASSP support:
- Standard: WAV
- Specialized: AU, KAY, NIST, NSP, CSRE, SSFF

---

## Usage Recommendations

### For Basic Speech Analysis:
- **Pitch**: `trk_rapt()` (fastest, robust)
- **Formants**: `trk_forest()` (established, reliable)
- **Energy**: `trk_rmsana()` (simple, fast)
- **Spectrum**: `trk_dftSpectrum()` (standard FFT)

### For Voice Quality Assessment:
- **Comprehensive**: `lst_vat()` (132 measures)
- **Clinical**: `lst_voice_sauce()` (40+ VoiceSauce parameters)
- **Standardized**: `lst_GeMAPS()` (62 features, international standard)
- **Extended**: `lst_eGeMAPS()` (88 features)

### For Prosodic Analysis:
- **Full Assessment**: `lst_dysprosody()` (193 features, MOMEL-INTSINT)
- **Complexity**: `lst_voxit()` (11 measures, articulation)

### For Real-Time Applications:
- **Pitch**: `trk_swiftf0()` (90-130ms, CNN-based)
- **Formants**: `trk_deepformants()` (2x realtime)
- **VAD**: `trk_brouhaha()` (50-100x optimized)

### For Noisy Speech:
- **Pitch**: `trk_swipe()`, `trk_sacc()` (noise-robust algorithms)
- **VAD+Quality**: `trk_brouhaha()` (SNR + C50 estimation)

### For Speech Synthesis:
- **F0**: `trk_dio()` or `trk_harvest()` (WORLD vocoder)
- **Spectrum**: `trk_seenc()` (CheapTrick spectral envelope)
- **Aperiodicity**: `trk_d4c()` (band aperiodicity)

---

## Function Count Summary

- **Track Functions (trk_*)**: 54 functions
- **Summary Functions (lst_*)**: 62+ functions (including implementation variants)
- **Total Documented Functions**: 75+

---

## Notes for pkgdown Configuration

### Suggested Reference Groups:

```yaml
reference:
- title: "Pitch & F0 Analysis"
  desc: "Fundamental frequency tracking algorithms"
  contents:
  - starts_with("trk_") & matches("rapt|swipe|dio|harvest|reaper|pitch|f0|swiftf0|crepe|sacc|yin|yaapt")

- title: "Formant Analysis"
  desc: "Formant frequency and bandwidth estimation"
  contents:
  - starts_with("trk_") & matches("forest|formant|deepformant|tvwlp|snackf")
  - lst_deepformants

- title: "Spectral Analysis"
  desc: "Spectrum estimation and cepstral analysis"
  contents:
  - starts_with("trk_") & matches("Spectrum|cepstrum|seenc|spectral_moments")

- title: "Energy & Amplitude"
  desc: "Signal energy and amplitude tracking"
  contents:
  - starts_with("trk_") & matches("rms|zcr|acf|intensity")

- title: "Voice Quality Measures"
  desc: "Voice quality assessment and dysphonia measures"
  contents:
  - lst_vat
  - lst_voice_sauce
  - lst_covarep_vq
  - starts_with("lst_voice_")
  - starts_with("lst_") & matches("avqi|dsi")
  - starts_with("trk_") & matches("d4c|brouhaha|creak|praat_sauce")

- title: "Prosody & Intonation"
  desc: "Prosodic feature extraction and complexity measures"
  contents:
  - lst_dysprosody
  - lst_voxit

- title: "Source-Filter Decomposition"
  desc: "Glottal source and vocal tract separation"
  contents:
  - starts_with("trk_") & matches("gfmiaif|iaif|excite")

- title: "OpenSMILE Feature Sets"
  desc: "Standardized acoustic feature extraction"
  contents:
  - starts_with("lst_") & matches("GeMAPS|eGeMAPS|emobase|ComParE|OpenSmile")

- title: "Acoustic Features & Coefficients"
  desc: "MFCC, LP coefficients, and other acoustic features"
  contents:
  - trk_mfcc
  - trk_lp_analysis
  - trk_npy_import

- title: "Pitch Marking & Epochs"
  desc: "Glottal closure instant detection"
  contents:
  - trk_pitchmark
```

---

**Document Created**: 2025-10-28
**Status**: Complete catalog of 75+ superassp DSP functions
**Purpose**: pkgdown reference organization and user guidance

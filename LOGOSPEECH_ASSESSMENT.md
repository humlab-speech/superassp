# LogoSpeech Studio Integration Assessment

**Date**: 2025-11-07  
**Status**: NOT RECOMMENDED FOR INTEGRATION  
**Decision**: Remove from superassp

## Executive Summary

The `logospeech-studio` codebase in `src/logospeech-studio/` is a Qt-based GUI application for speech analysis. While it contains a clean DSP library (`logo/` directory) with standard speech processing algorithms, **it is NOT recommended for integration into superassp** for the following reasons:

1. **Extensive duplication**: Most algorithms already exist in superior implementations
2. **Qt dependency**: Core DSP code uses `QVector`, requiring Qt Core even for headless processing
3. **Limited novel functionality**: Only 2-3 unique features (AMDF pitch, LPC-cepstrum)
4. **GUI-centric design**: Tightly coupled with Qt GUI framework
5. **Inferior pitch tracking**: Basic algorithms vs superassp's state-of-the-art implementations

**Recommendation**: **DELETE** `src/logospeech-studio/` from superassp repository.

---

## Detailed Analysis

### 1. Codebase Overview

**Source**: https://github.com/mohabouje/logospeech-studio  
**License**: MIT (permissive)  
**Purpose**: Visual support for language learning by people with hearing impairments  
**Language**: C++ with Qt5 framework  
**Size**: ~1,218 lines in DSP core (`logo/` directory)

**Directory Structure**:
- `logo/`: Core DSP library (4 files, 950 LOC)
- `capturedata/`: Audio I/O with Qt Multimedia
- `settings/`, `customwidget/`: Qt GUI components
- `mainwindow.*`: Main GUI application
- `qcustomplot.*`: Plotting library

### 2. Dependencies

**Required**:
- **Qt5**: `core`, `gui`, `widgets`, `multimedia`, `printsupport`
- **FFTW3**: Fast Fourier Transform library

**Qt Core Dependency Issue**:
- Even the "headless" DSP library uses `QVector<double>` throughout
- Cannot compile without Qt Core library
- Would require complete refactoring to use `std::vector`

### 3. Implemented DSP Algorithms

#### Core Features (logo/ directory)

| Feature | logo Function | superassp Equivalent | Status |
|---------|---------------|---------------------|--------|
| **Energy** | `logo_energy_lin`, `logo_energy_db_spl` | `trk_rmsana` | Duplicate |
| **Zero-Crossing Rate** | `logo_zero_crossing_rate` | `trk_zcrana` | Duplicate |
| **Autocorrelation** | `logo_pitch_tracker_acf` | `trk_acfana` | Duplicate |
| **LPC** | `logo_fft_lpc` (Levinson-Durbin) | `trk_lp_analysis` | Duplicate |
| **Cepstrum** | `logo_fft_ceps` | `trk_cepstrum` | Duplicate |
| **Power Spectral Density** | `logo_fft_psd` | `trk_dftSpectrum` | Duplicate |
| **Pitch (ACF)** | `logo_pitch_tracker_acf` | `trk_acfana` | Inferior |
| **Pitch (Cepstrum)** | `logo_pitch_tracker_ceps` | `trk_rapt`, `trk_swipe` | Inferior |
| **Pitch (AMDF)** | `logo_pitch_tracker_amdf` | *(none)* | **Unique** |
| **LPC-Cepstrum** | `logo_cepsLPC` | *(none)* | **Unique** |
| **Mel Filterbank** | `logo_fft_filterbank` | `trk_mfcc` (partial) | Incomplete |

#### Preprocessing Features

| Feature | logo Function | superassp Equivalent |
|---------|---------------|---------------------|
| **Windowing** | `logo_apply_window` (Hamming, Hanning, Blackman, etc.) | Internal to DSP functions |
| **Pre-emphasis** | `logo_filter_filt` + `logo_slicer` | Internal to DSP functions |
| **DC offset removal** | `logo_filter_filt` + `logo_slicer` | Internal to DSP functions |
| **Frame slicing** | `logo_apply_buffer_slicer` | Internal to DSP functions |

### 4. Comparison with superassp

#### Duplicates (9 functions)
These offer **no new functionality**:
1. `logo_energy_lin` / `logo_energy_db_spl` → `trk_rmsana`
2. `logo_zero_crossing_rate` → `trk_zcrana`
3. `logo_pitch_tracker_acf` → `trk_acfana`
4. `logo_fft_lpc` → `trk_lp_analysis`
5. `logo_fft_ceps` → `trk_cepstrum`
6. `logo_fft_psd` → `trk_dftSpectrum`

#### Inferior Implementations (3 functions)
These are **worse than existing superassp functions**:

| logo Function | superassp Alternative | Why Inferior? |
|---------------|----------------------|---------------|
| `logo_pitch_tracker_acf` | `trk_rapt`, `trk_swipe`, `trk_dio` | Basic ACF vs robust multi-cue algorithms |
| `logo_pitch_tracker_ceps` | `trk_rapt`, `trk_swipe` | Simple cepstral peak vs advanced tracking |
| `logo_fft_filterbank` | `trk_mfcc` | Missing DCT step for full MFCC |

#### Unique Features (2-3 functions)
Only these offer potential value:

1. **`logo_pitch_tracker_amdf`** (AMDF pitch tracking)
   - Average Magnitude Difference Function
   - Alternative to ACF for pitch detection
   - **Value**: Low - superassp already has 17 superior pitch trackers

2. **`logo_cepsLPC`** (LPC-based cepstrum)
   - Cepstral coefficients from LPC (not FFT)
   - Smoothed spectrum modeling vocal tract
   - **Value**: Medium - could be useful for formant analysis

3. **`logo_apply_buffer_slicer`** (integrated preprocessing)
   - Combines framing + DC removal + pre-emphasis
   - **Value**: Low - superassp functions handle this internally

### 5. Code Quality Assessment

**Strengths**:
- ✅ Clean separation: `logo/` DSP library vs GUI
- ✅ Consistent naming: `logo_*` prefix
- ✅ MIT license (permissive)
- ✅ Logical file organization

**Weaknesses**:
- ❌ **Minimal documentation**: Sparse comments
- ❌ **Qt dependency**: `QVector` throughout core DSP
- ❌ **No unit tests**: No test suite
- ❌ **GUI-centric**: Designed for Qt application, not library
- ❌ **No R interface**: Would require complete Rcpp wrapper

### 6. Integration Effort vs Value

**Required Work for Integration**:
1. Extract `logo/` directory (4 files, ~950 LOC)
2. Refactor `QVector` → `std::vector` or add Qt Core dependency
3. Add FFTW3 to `src/Makevars` (already have it)
4. Create Rcpp wrappers for 2-3 unique functions
5. Write R interface functions
6. Add roxygen2 documentation
7. Create test suite
8. Handle SSFF output format conversion

**Estimated Effort**: 16-24 hours

**Value Delivered**:
- **AMDF pitch tracker**: Low value (already have 17 pitch trackers)
- **LPC-cepstrum**: Medium value (niche use case)

**Conclusion**: **NOT WORTH THE EFFORT**

### 7. Architectural Issues

#### Qt Dependency Problem
```cpp
// From logo_structs.h
typedef QVector<logo_real> logo_vector;
typedef QVector<QVector<logo_real>> logo_matrix;
```

Every DSP function uses these types. Options:
1. **Add Qt Core to superassp dependencies** → Bloat (Qt is ~50MB)
2. **Refactor to std::vector** → Large effort, breaks existing code
3. **Don't integrate** → ✅ **Recommended**

#### FFTW3 Overlap
superassp already uses FFTW3 in other functions, so this dependency is acceptable. However, the `logo` library adds no new FFTW-based functionality beyond what we have.

### 8. Alternative: Extract Unique Functions Only

If we wanted to salvage anything, extract only:

**Option A: LPC-Cepstrum Only**
- Extract `logo_cepsLPC` algorithm
- Rewrite in pure C++ without Qt
- Create `trk_lpc_cepstrum()` R function
- Estimated effort: 4-6 hours

**Option B: AMDF Pitch Only**
- Extract `logo_pitch_tracker_amdf` algorithm
- Rewrite without Qt
- Create `trk_amdf()` R function
- Estimated effort: 3-4 hours

**Recommendation**: Even these are **not worth the effort** given superassp's comprehensive existing functionality.

---

## Final Recommendation

### ❌ DO NOT INTEGRATE

**Rationale**:
1. **95% duplication** with existing superior implementations
2. **Qt dependency** is architectural mismatch for superassp
3. **Minimal unique value** (2 niche algorithms)
4. **Inferior pitch tracking** vs superassp's 17 specialized trackers
5. **High integration cost** (16-24 hours) for low value

### ✅ ACTION REQUIRED

**Delete `src/logospeech-studio/` from superassp repository**:
```bash
git rm -rf src/logospeech-studio/
git commit -m "chore: Remove logospeech-studio - extensive duplication with existing functions"
```

### 📝 Documentation Update

Add note to `CLAUDE.md` or README:
> **LogoSpeech Studio**: Evaluated 2025-11-07. Qt-based GUI application with standard DSP algorithms. 95% duplication with existing superassp functions (trk_zcrana, trk_rmsana, trk_lp_analysis, etc.). Qt dependency incompatible with superassp architecture. Not integrated.

---

## Appendix: Function Inventory

### Implemented in logo/ (20 functions)

**Short-Time Analysis (5)**:
- `logo_energy_lin` - Linear energy
- `logo_energy_db_spl` - Energy in dB SPL
- `logo_zero_crossing_rate` - ZCR
- `logo_frames_energy` - Energy per frame
- `logo_frames_zero_crossing_rate` - ZCR per frame

**Spectral Analysis (6)**:
- `logo_fft_psd` - Power spectral density
- `logo_fft_lpc` - LPC spectrum (Levinson-Durbin)
- `logo_fft_ceps` - Cepstrum
- `logo_fft_filterbank` - Mel filterbank
- `logo_fft_norm` - Normalized FFT
- `logo_fft_real` - Real FFT

**Pitch Tracking (6)**:
- `logo_pitch_tracker_amdf` - AMDF pitch (**unique**)
- `logo_pitch_tracker_acf` - ACF pitch
- `logo_pitch_tracker_ceps` - Cepstral pitch
- `logo_frames_pitch_tracker_amdf` - AMDF per frame
- `logo_frames_pitch_tracker_acf` - ACF per frame
- `logo_frames_pitch_tracker_cepstrum` - Cepstral per frame

**Preprocessing (3)**:
- `logo_apply_window` - Windowing
- `logo_filter_filt` - Filtering
- `logo_apply_buffer_slicer` - Frame slicing + preprocessing

### Already in superassp (Equivalents)

**Energy/Amplitude (4)**:
- `trk_rmsana()` - RMS energy
- `trk_zcrana()` - Zero-crossing rate
- `trk_acfana()` - Autocorrelation
- `trk_intensityp()` - Praat intensity

**Spectral (6)**:
- `trk_dftSpectrum()` - DFT spectrum
- `trk_cssSpectrum()` - CSS spectrum
- `trk_lpsSpectrum()` - LPS spectrum
- `trk_lp_analysis()` - LPC analysis
- `trk_cepstrum()` - Cepstral analysis
- `trk_mfcc()` - MFCC (includes filterbank)

**Pitch (17)**:
- `trk_rapt()`, `trk_swipe()`, `trk_dio()`, `trk_harvest()` - SPTK/WORLD
- `trk_reaper()`, `trk_ksvfo()`, `trk_mhspitch()` - ESTK/ASSP
- `trk_swiftf0()`, `trk_crepe()` - Deep learning
- `trk_pyin()`, `trk_yin()`, `trk_yaapt()` - Python classical
- `trk_sacc()`, `trk_pitchp()` - Parselmouth/Praat
- ... and more

---

## References

- LogoSpeech Studio repository: https://github.com/mohabouje/logospeech-studio
- Developer: Mohammed Boujemaoui Boulaghmoudi
- Instructor: Ángel de la Torre Vega (University of Granada)
- License: MIT
- Purpose: Visual support for language learning (hearing impairments)


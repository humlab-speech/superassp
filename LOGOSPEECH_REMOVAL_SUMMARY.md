# LogoSpeech Studio Removal Summary

**Date**: 2025-11-07  
**Action**: Removed from superassp repository  
**Commit**: 8c8a8eb

## Decision

**NOT INTEGRATED** - Removed `src/logospeech-studio/` from superassp.

## Rationale

The logospeech-studio codebase was assessed for integration into superassp and found to be unsuitable:

### Key Issues

1. **95% Functional Duplication**
   - Most algorithms already exist in superior implementations
   - Energy, ZCR, LPC, Cepstrum, Autocorrelation all duplicates

2. **Architectural Incompatibility**
   - Qt dependency throughout core DSP (`QVector<double>`)
   - Would require Qt Core library even for headless processing
   - Incompatible with superassp's lightweight C++/R architecture

3. **Inferior Implementations**
   - Basic pitch tracking (ACF, Cepstrum) vs superassp's 17 specialized algorithms
   - No advantage over existing `trk_rapt`, `trk_swipe`, `trk_dio`, etc.

4. **Minimal Unique Value**
   - Only 2 unique features: AMDF pitch tracking, LPC-cepstrum
   - Not worth 16-24 hour integration effort

## What LogoSpeech Offered

### DSP Algorithms (20 functions)
- Energy estimation (linear, dB SPL)
- Zero-crossing rate
- FFT-based spectral analysis (PSD, LPC, Cepstrum)
- Pitch tracking (ACF, AMDF, Cepstrum)
- Mel filterbank
- Windowing and preprocessing

### superassp Equivalents (Already Implemented)

| LogoSpeech Function | superassp Equivalent | Status |
|---------------------|----------------------|--------|
| `logo_energy_lin/db_spl` | `trk_rmsana()` | Superior |
| `logo_zero_crossing_rate` | `trk_zcrana()` | Equivalent |
| `logo_pitch_tracker_acf` | `trk_acfana()` | Equivalent |
| `logo_fft_lpc` | `trk_lp_analysis()` | Equivalent |
| `logo_fft_ceps` | `trk_cepstrum()` | Equivalent |
| `logo_fft_psd` | `trk_dftSpectrum()` | Equivalent |
| `logo_pitch_tracker_ceps` | `trk_rapt/swipe/dio()` | **Inferior** |
| `logo_fft_filterbank` | `trk_mfcc()` | Incomplete (no DCT) |

## Unique Features (Not Worth Integration)

1. **AMDF Pitch Tracker** (`logo_pitch_tracker_amdf`)
   - Average Magnitude Difference Function
   - Alternative to ACF
   - **Value**: Low (superassp has 17 pitch trackers already)

2. **LPC-Cepstrum** (`logo_cepsLPC`)
   - Cepstral coefficients from LPC (not FFT)
   - Vocal tract transfer function modeling
   - **Value**: Medium (niche use case, not worth effort)

## Technical Assessment

**Language**: C++ with Qt5  
**License**: MIT (permissive)  
**Size**: ~1,218 lines in DSP core  
**Dependencies**: Qt5 Core/GUI, FFTW3  
**Code Quality**: Clean structure, minimal documentation, no tests  

**Integration Effort**: 16-24 hours  
**Value Delivered**: Minimal (2 niche algorithms)  
**Conclusion**: **Not cost-effective**

## Comparison with superassp

superassp already provides comprehensive coverage:

### Pitch Tracking (17 algorithms)
- **C++ SPTK/WORLD**: `trk_rapt`, `trk_swipe`, `trk_dio`, `trk_harvest`, `trk_reaper`
- **C ASSP/ESTK**: `trk_ksvfo`, `trk_mhspitch`, `trk_estk_pitchmark`
- **Deep Learning**: `trk_swiftf0`, `trk_crepe`
- **Python Classical**: `trk_pyin`, `trk_yin`, `trk_yaapt`
- **Parselmouth**: `trk_sacc`, `trk_pitchp`

### Spectral Analysis (6 algorithms)
- `trk_dftSpectrum`, `trk_cssSpectrum`, `trk_lpsSpectrum`
- `trk_lp_analysis`, `trk_cepstrum`, `trk_mfcc`

### Energy/Amplitude (4 algorithms)
- `trk_rmsana`, `trk_zcrana`, `trk_acfana`, `trk_intensityp`

LogoSpeech offers nothing beyond what superassp already has in superior form.

## Documentation

Full assessment in **LOGOSPEECH_ASSESSMENT.md** includes:
- Detailed function comparison
- Code quality analysis
- Integration effort estimation
- Architectural compatibility analysis

## References

- **Repository**: https://github.com/mohabouje/logospeech-studio
- **Developer**: Mohammed Boujemaoui Boulaghmoudi
- **Instructor**: Ángel de la Torre Vega (University of Granada)
- **Purpose**: Visual support for language learning (hearing impairments)
- **License**: MIT

## Recommendation for Similar Codebases

When evaluating external DSP libraries for integration:

1. ✅ Check for unique algorithms not in superassp
2. ✅ Assess code quality and documentation
3. ✅ Verify licensing compatibility (MIT, BSD, Apache 2.0)
4. ✅ Consider dependencies (avoid Qt, prefer standard C++)
5. ✅ Compare performance vs existing implementations
6. ✅ Estimate integration effort vs value delivered

**LogoSpeech failed criteria**: Extensive duplication, Qt dependency, minimal unique value.


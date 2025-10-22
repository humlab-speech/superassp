# Track Name Inventory

Generated: 2025-10-22 17:06:50.893938

## Summary Statistics

- **Total track definitions**: 111
- **Unique old names**: 74
- **Unique new names**: 71
- **Functions affected**: 53
- **Tracks needing changes**: 36 (32.4%)
- **Tracks unchanged**: 75 (67.6%)

## Changes Required

- **Add units**: 35 tracks
- **Case normalization**: 0 tracks
- **Add hyphens**: 0 tracks (harmonic differences)

## Tracks by Category

| Category | Count | Percentage |
|----------|-------|------------|
| other | 45 | 40.5% |
| f0 | 24 | 21.6% |
| metadata | 16 | 14.4% |
| voice_quality | 11 | 9.9% |
| spectral | 10 | 9.0% |
| bandwidth | 2 | 1.8% |
| formant | 2 | 1.8% |
| pitch_mark | 1 | 0.9% |

## Functions with Most Track Definitions (Top 20)

| Function | Track Count | File |
|----------|-------------|------|
| `lst_voice_reportp` | 26 | list_python_pm_pvoice_report.R |
| `trk_spectral_momentsp` | 4 | ssff_python_pm_pspectral_moments.R |
| `arfana` | 3 | ssff_c_assp_lp_analysis.R |
| `larana` | 3 | ssff_c_assp_lp_analysis.R |
| `lpcana` | 3 | ssff_c_assp_lp_analysis.R |
| `rfcana` | 3 | ssff_c_assp_lp_analysis.R |
| `trk_covarep_srh` | 3 | covarep_srh.R |
| `trk_formantp` | 3 | ssff_python_pm_pformantb.R |
| `trk_formantpathp` | 3 | ssff_python_pm_pformantpathb.R |
| `trk_pyin` | 3 | ssff_python_pyin.R |
| `trk_snackp` | 3 | ssff_python_snack_pitch.R |
| `trk_swiftf0` | 3 | ssff_python_swiftf0.R |
| `nonopt_rapt` | 2 | ssff_python_sptk_rapt.R |
| `nonopt_reaper` | 2 | ssff_python_sptk_reaper.R |
| `nonopt_swipe` | 2 | ssff_python_sptk_swipe.R |
| `trk_covarep_iaif` | 2 | covarep_iaif.R |
| `trk_crepe` | 2 | ssff_python_crepe.R |
| `trk_dio` | 2 | ssff_cpp_sptk_dio.R |
| `trk_forest` | 2 | ssff_c_assp_forest.R |
| `trk_harvest` | 2 | ssff_python_world_harvest.R |

## Category: f0 (24 tracks)

| Old Name | New Name | Unit | Notes |
|----------|----------|------|-------|
| `f0` | `fo[Hz]` | Hz | add unit |
| `F0` | `fo[Hz]` | Hz | lowercase, add unit |
| `F0[Hz]` | `fo[Hz]` | Hz | lowercase |
| `fo[Hz]` | `fo[Hz]` | Hz | Already correct (Titze 2015) |

**Functions** (22): `fo`, `fo_ksv`, `foana`, `nonopt_rapt`, `nonopt_reaper`, `nonopt_swipe`, `trk_covarep_srh`, `trk_crepe`, `trk_dio`, `trk_harvest`, `trk_kaldi_pitch`, `trk_ksvfo`, `trk_npy_import`, `trk_pyin`, `trk_rapt`, `trk_reaper`, `trk_snackp`, `trk_swiftf0`, `trk_swipe`, `trk_torch_pitch`, `trk_yaapt`, `trk_yin`

## Category: formant (2 tracks)

| Old Name | New Name | Unit | Notes |
|----------|----------|------|-------|
| `fm` | `Fi[Hz]` | Hz | add unit |

**Functions** (2): `trk_formantp`, `trk_formantpathp`

## Category: bandwidth (2 tracks)

| Old Name | New Name | Unit | Notes |
|----------|----------|------|-------|
| `bw` | `Bi[Hz]` | Hz | add unit |

**Functions** (2): `trk_formantp`, `trk_formantpathp`

## Category: voice_quality (11 tracks)

| Old Name | New Name | Unit | Notes |
|----------|----------|------|-------|
| `Jitter (ddp)` | `Jitter_DDP[%]` | % | add unit, remove parentheses |
| `Jitter (local, absolute)` | `Jitter_local_abs[us]` | us | add unit, remove parentheses |
| `Jitter (local)` | `Jitter_local[%]` | % | add unit, remove parentheses |
| `Jitter (ppq5)` | `Jitter_PPQ5[%]` | % | add unit, remove parentheses |
| `Jitter (rap)` | `Jitter_RAP[%]` | % | add unit, remove parentheses |
| `Shimmer (apq11)` | `Shimmer_APQ11[%]` | % | add unit, remove parentheses |
| `Shimmer (apq3)` | `Shimmer_APQ3[%]` | % | add unit, remove parentheses |
| `Shimmer (apq5)` | `Shimmer_APQ5[%]` | % | add unit, remove parentheses |
| `Shimmer (dda)` | `Shimmer_DDA[%]` | % | add unit, remove parentheses |
| `Shimmer (local, dB)` | `Shimmer_local[dB]` | dB | add unit, remove parentheses |
| `Shimmer (local)` | `Shimmer_local[%]` | % | add unit, remove parentheses |

**Functions** (1): `lst_voice_reportp`

## Appendix: Complete Mapping Table

| Function | File | Old Name | New Name | Category | Unit | Notes |
|----------|------|----------|----------|----------|------|-------|
| `trk_formantp` | ssff_python_pm_pformantb.R | `bw` | `Bi[Hz]` | bandwidth | Hz | add unit |
| `trk_formantpathp` | ssff_python_pm_pformantpathb.R | `bw` | `Bi[Hz]` | bandwidth | Hz | add unit |
| `fo` | ssff_c_assp_ksvfo.R | `fo[Hz]` | `fo[Hz]` | f0 | Hz | Already correct (Titze 2015) |
| `fo_ksv` | ssff_c_assp_ksvfo.R | `fo[Hz]` | `fo[Hz]` | f0 | Hz | Already correct (Titze 2015) |
| `foana` | ssff_c_assp_ksvfo.R | `fo[Hz]` | `fo[Hz]` | f0 | Hz | Already correct (Titze 2015) |
| `nonopt_rapt` | ssff_python_sptk_rapt.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `nonopt_reaper` | ssff_python_sptk_reaper.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `nonopt_swipe` | ssff_python_sptk_swipe.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `trk_covarep_srh` | covarep_srh.R | `F0[Hz]` | `fo[Hz]` | f0 | Hz | lowercase |
| `trk_crepe` | ssff_python_crepe.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `trk_dio` | ssff_cpp_sptk_dio.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `trk_dio` | ssff_python_world_dio.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `trk_harvest` | ssff_python_world_harvest.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `trk_harvest` | ssff_cpp_sptk_harvest.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `trk_kaldi_pitch` | ssff_python_kaldi_pitch.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `trk_ksvfo` | ssff_c_assp_ksvfo.R | `fo[Hz]` | `fo[Hz]` | f0 | Hz | Already correct (Titze 2015) |
| `trk_npy_import` | ssff_python_npy_import.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `trk_pyin` | ssff_python_pyin.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `trk_rapt` | ssff_cpp_sptk_rapt.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `trk_reaper` | ssff_cpp_sptk_reaper.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `trk_snackp` | ssff_python_snack_pitch.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `trk_swiftf0` | ssff_python_swiftf0.R | `F0` | `fo[Hz]` | f0 | Hz | lowercase, add unit |
| `trk_swipe` | ssff_cpp_sptk_swipe.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `trk_torch_pitch` | ssff_python_torch_pitch.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `trk_yaapt` | ssff_python_yaapt.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `trk_yin` | ssff_python_yin.R | `f0` | `fo[Hz]` | f0 | Hz | add unit |
| `trk_formantp` | ssff_python_pm_pformantb.R | `fm` | `Fi[Hz]` | formant | Hz | add unit |
| `trk_formantpathp` | ssff_python_pm_pformantpathb.R | `fm` | `Fi[Hz]` | formant | Hz | add unit |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Degree of voice breaks` | `Degree of voice breaks` | metadata | — | No change needed |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Fraction of locally unvoiced frames` | `Fraction of locally unvoiced frames` | metadata | — | No change needed |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Mean autocorrelation` | `Mean autocorrelation` | metadata | — | No change needed |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Mean period` | `Mean period` | metadata | — | No change needed |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Number of periods` | `Number of periods` | metadata | — | No change needed |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Number of voice breaks` | `Number of voice breaks` | metadata | — | No change needed |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Standard deviation of period` | `Standard deviation of period` | metadata | — | No change needed |
| `nonopt_reaper` | ssff_python_sptk_reaper.R | `corr` | `corr` | metadata | — | No change needed |
| `trk_aperiodicities` | ssff_python_aperiodicities.R | `aperiod` | `aperiod` | metadata | — | No change needed |
| `trk_crepe` | ssff_python_crepe.R | `periodicity` | `periodicity` | metadata | — | No change needed |
| `trk_d4c` | ssff_cpp_sptk_d4c.R | `aperiodicity` | `aperiodicity` | metadata | — | No change needed |
| `trk_pyin` | ssff_python_pyin.R | `voiced` | `voiced` | metadata | — | No change needed |
| `trk_snackp` | ssff_python_snack_pitch.R | `voicing` | `voicing` | metadata | — | No change needed |
| `trk_swiftf0` | ssff_python_swiftf0.R | `confidence` | `confidence` | metadata | — | No change needed |
| `trk_swiftf0` | ssff_python_swiftf0.R | `voicing` | `voicing` | metadata | — | No change needed |
| `trk_yaapt` | ssff_python_yaapt.R | `voiced` | `voiced` | metadata | — | No change needed |
| `arfana` | ssff_c_assp_lp_analysis.R | `ARF` | `ARF` | other | — | No change needed |
| `arfana` | ssff_c_assp_lp_analysis.R | `gain[dB]` | `gain[dB]` | other | dB | No change needed |
| `larana` | ssff_c_assp_lp_analysis.R | `gain[dB]` | `gain[dB]` | other | dB | No change needed |
| `larana` | ssff_c_assp_lp_analysis.R | `LAR` | `LAR` | other | — | No change needed |
| `lpcana` | ssff_c_assp_lp_analysis.R | `gain[dB]` | `gain[dB]` | other | dB | No change needed |
| `lpcana` | ssff_c_assp_lp_analysis.R | `LPC` | `LPC` | other | — | No change needed |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Maximum pitch` | `Maximum pitch` | other | — | No change needed |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Mean harmonics-to-noise ratio` | `Mean harmonics-to-noise ratio` | other | — | No change needed |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Mean noise-to-harmonics ratio` | `Mean noise-to-harmonics ratio` | other | — | No change needed |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Mean pitch` | `Mean pitch` | other | — | No change needed |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Median pitch` | `Median pitch` | other | — | No change needed |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Minimum pitch` | `Minimum pitch` | other | — | No change needed |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Number of pulses` | `Number of pulses` | other | — | No change needed |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Standard deviation` | `Standard deviation` | other | — | No change needed |
| `nonopt_rapt` | ssff_python_sptk_rapt.R | `pitch` | `pitch` | other | — | No change needed |
| `nonopt_swipe` | ssff_python_sptk_swipe.R | `pitch` | `pitch` | other | — | No change needed |
| `pitch` | ssff_c_assp_mhspitch.R | `pitch[Hz]` | `pitch[Hz]` | other | Hz | No change needed |
| `pitch_mhs` | ssff_c_assp_mhspitch.R | `pitch[Hz]` | `pitch[Hz]` | other | Hz | No change needed |
| `rfcana` | ssff_c_assp_lp_analysis.R | `gain[dB]` | `gain[dB]` | other | dB | No change needed |
| `rfcana` | ssff_c_assp_lp_analysis.R | `RFC` | `RFC` | other | — | No change needed |
| `trk_cepstrum` | ssff_c_assp_cepstrum.R | `C[dB]` | `C[dB]` | other | dB | No change needed |
| `trk_covarep_iaif` | covarep_iaif.R | `glottal_derivative` | `glottal_derivative` | other | — | No change needed |
| `trk_covarep_iaif` | covarep_iaif.R | `glottal_flow` | `glottal_flow` | other | — | No change needed |
| `trk_covarep_srh` | covarep_srh.R | `SRH` | `SRH` | other | — | No change needed |
| `trk_covarep_srh` | covarep_srh.R | `VUV` | `VUV` | other | — | No change needed |
| `trk_cssSpectrum` | ssff_c_assp_cssSpectrum.R | `CSS[dB]` | `CSS[dB]` | other | dB | No change needed |
| `trk_dftSpectrum` | ssff_c_assp_dftSpectrum.R | `DFT[dB]` | `DFT[dB]` | other | dB | No change needed |
| `trk_excite` | ssff_python_excite.R | `excitation` | `excitation` | other | — | No change needed |
| `trk_forest` | ssff_c_assp_forest.R | `B[Hz]` | `B[Hz]` | other | Hz | No change needed |
| `trk_forest` | ssff_c_assp_forest.R | `F[Hz]` | `F[Hz]` | other | Hz | No change needed |
| `trk_formantp` | ssff_python_pm_pformantb.R | `lv` | `lv` | other | — | No change needed |
| `trk_formantpathp` | ssff_python_pm_pformantpathb.R | `lv` | `lv` | other | — | No change needed |
| `trk_intensityp` | ssff_python_pm_pintensity.R | `intensity` | `intensity` | other | — | No change needed |
| `trk_lpsSpectrum` | ssff_c_assp_lpsSpectrum.R | `LPS[dB]` | `LPS[dB]` | other | dB | No change needed |
| `trk_mhspitch` | ssff_c_assp_mhspitch.R | `pitch[Hz]` | `pitch[Hz]` | other | Hz | No change needed |
| `trk_pitchmark` | ssff_cpp_estk_pitchmark.R | `pitchmarks` | `pitchmarks` | other | — | No change needed |
| `trk_pitchp` | ssff_python_pm_ppitch.R | `ac` | `ac` | other | — | No change needed |
| `trk_pitchp` | ssff_python_pm_ppitch.R | `cc` | `cc` | other | — | No change needed |
| `trk_pyin` | ssff_python_pyin.R | `vprob` | `vprob` | other | — | No change needed |
| `trk_seenc` | ssff_python_seenc.R | `trk_seenc` | `trk_seenc` | other | — | No change needed |
| `trk_snackp` | ssff_python_snack_pitch.R | `rms` | `rms` | other | — | No change needed |
| `trk_spectral_momentsp` | ssff_python_pm_pspectral_moments.R | `cog` | `cog` | other | — | No change needed |
| `trk_spectral_momentsp` | ssff_python_pm_pspectral_moments.R | `kurtosis` | `kurtosis` | other | — | No change needed |
| `trk_spectral_momentsp` | ssff_python_pm_pspectral_moments.R | `sd` | `sd` | other | — | No change needed |
| `trk_spectral_momentsp` | ssff_python_pm_pspectral_moments.R | `skewness` | `skewness` | other | — | No change needed |
| `reaper_pm` | ssff_python_reaper_pm.R | `pm` | `PM[s]` | pitch_mark | s | add unit |
| `arfana` | ssff_c_assp_lp_analysis.R | `RMS[dB]` | `RMS[dB]` | spectral | dB | No change needed |
| `larana` | ssff_c_assp_lp_analysis.R | `RMS[dB]` | `RMS[dB]` | spectral | dB | No change needed |
| `lpcana` | ssff_c_assp_lp_analysis.R | `RMS[dB]` | `RMS[dB]` | spectral | dB | No change needed |
| `mfcc` | ssff_python_torch_mfcc.R | `DYNAMIC` | `DYNAMIC` | spectral | — | No change needed |
| `rfcana` | ssff_c_assp_lp_analysis.R | `RMS[dB]` | `RMS[dB]` | spectral | dB | No change needed |
| `trk_acfana` | ssff_c_assp_acfana.R | `ACF` | `ACF` | spectral | — | No change needed |
| `trk_mfcc` | ssff_cpp_sptk_mfcc.R | `DYNAMIC` | `DYNAMIC` | spectral | — | No change needed |
| `trk_rmsana` | ssff_c_assp_rmsana.R | `RMS[dB]` | `RMS[dB]` | spectral | dB | No change needed |
| `trk_snackf` | ssff_python_snack_formant.R | `DYNAMIC` | `DYNAMIC` | spectral | — | No change needed |
| `trk_zcrana` | ssff_c_assp_zcrana.R | `ZCR[Hz]` | `ZCR[Hz]` | spectral | Hz | No change needed |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Jitter (ddp)` | `Jitter_DDP[%]` | voice_quality | % | add unit, remove parentheses |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Jitter (local, absolute)` | `Jitter_local_abs[us]` | voice_quality | us | add unit, remove parentheses |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Jitter (local)` | `Jitter_local[%]` | voice_quality | % | add unit, remove parentheses |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Jitter (ppq5)` | `Jitter_PPQ5[%]` | voice_quality | % | add unit, remove parentheses |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Jitter (rap)` | `Jitter_RAP[%]` | voice_quality | % | add unit, remove parentheses |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Shimmer (apq11)` | `Shimmer_APQ11[%]` | voice_quality | % | add unit, remove parentheses |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Shimmer (apq3)` | `Shimmer_APQ3[%]` | voice_quality | % | add unit, remove parentheses |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Shimmer (apq5)` | `Shimmer_APQ5[%]` | voice_quality | % | add unit, remove parentheses |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Shimmer (dda)` | `Shimmer_DDA[%]` | voice_quality | % | add unit, remove parentheses |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Shimmer (local, dB)` | `Shimmer_local[dB]` | voice_quality | dB | add unit, remove parentheses |
| `lst_voice_reportp` | list_python_pm_pvoice_report.R | `Shimmer (local)` | `Shimmer_local[%]` | voice_quality | % | add unit, remove parentheses |

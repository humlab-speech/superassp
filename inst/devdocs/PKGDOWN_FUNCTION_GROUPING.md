# superassp v1.2.1 — Exported Function Reference

**Updated**: 2026-02-28
**Exports**: 93 (from NAMESPACE) + 16 S3method registrations

---

## 1. Pitch & F0 Tracking (20 functions)

| Function | Engine | Speed | Description |
|----------|--------|-------|-------------|
| `trk_pitch_rapt()` | C++ SPTK | Fastest | RAPT dynamic-programming pitch tracker |
| `trk_pitch_swipe()` | C++ SPTK | Fastest | SWIPE sawtooth-inspired estimator (noise-robust) |
| `trk_pitch_dio()` | C++ WORLD | Fastest | DIO synthesis-quality F0 |
| `trk_pitch_harvest()` | C++ WORLD | Fastest | Harvest high-accuracy F0 |
| `trk_pitch_reaper()` | C++ SPTK | Fastest | REAPER F0 + epoch detection |
| `trk_pitch_yin()` | Python | Moderate | Classic YIN autocorrelation |
| `trk_pitch_pyin()` | Python | Moderate | Probabilistic YIN with HMM Viterbi |
| `trk_yaapt()` | Python | Moderate | NCC + dynamic programming |
| `trk_pitch_pda()` | C++ ESTK | Fast | Edinburgh Speech Tools PDA |
| `trk_tandem()` | C++ | Fast | Tandem-STRAIGHT pitch |
| `trk_swiftf0()` | Python DL | Fast | CNN real-time pitch (90-130 ms) |
| `trk_pitch_crepe()` | Python DL | Moderate | Deep CNN on raw waveform |
| `trk_sacc()` | Python | Moderate | Subband autocorrelation (noise-robust) |
| `trk_pitch_cc()` | pladdrr | Fast | Praat pitch (cross-correlation) |
| `trk_pitch_ac()` | pladdrr | Fast | Praat pitch (autocorrelation) |
| `trk_pitch_shs()` | pladdrr | Fast | Praat pitch (subharmonic summation) |
| `trk_pitch_spinet()` | pladdrr | Fast | Praat pitch (SPINET) |
| `trk_dv_f0()` | Python | Moderate | DisVoice F0 tracking |
| `trk_pitch_snack()` | Python | Moderate | Snack Toolkit pitch |
| `trk_covarep_srh()` | Python | Moderate | COVAREP Spectral Relative Harmonics |
| `trk_srh()` | Python | Moderate | VAT SRH variant |
| `trk_straight_f0()` | Python | Moderate | Legacy STRAIGHT F0 extraction |
| `trk_egg_f0()` | Python | Moderate | EGG-based F0 |

---

## 2. Formant Analysis (7 functions)

| Function | Engine | Description |
|----------|--------|-------------|
| `trk_formant_forest()` | C ASSP | Linear prediction (autocorrelation + SLA) |
| `trk_formant_burg()` | pladdrr | Praat Burg method with HMM tracking |
| `trk_deepformants()` | Python DL | PyTorch RNN (2x real-time) |
| `trk_formants_tvwlp()` | Python | Time-Varying Weighted LP (GCI-based) |
| `trk_formant_snack()` | Python | Snack Toolkit formants |
| `trk_dv_formants()` | Python | DisVoice formant tracking |
| `lst_deepformants()` | Python DL | Summary statistics of deep-learning formant tracks |

---

## 3. Spectral Analysis (8 functions)

| Function | Engine | Description |
|----------|--------|-------------|
| `trk_dft_spectrum()` | C ASSP | Discrete Fourier Transform (narrow-band) |
| `trk_css_spectrum()` | C ASSP | Cepstrally smoothed spectrum |
| `trk_lps_spectrum()` | C ASSP | Linear prediction smoothed spectrum |
| `trk_cepstrum()` | C ASSP | Short-term cepstral coefficients |
| `trk_seenc()` | Python | Spectral envelope (WORLD CheapTrick) |
| `trk_straight_spec()` | Python | STRAIGHT spectral analysis |
| `trk_spectral_moments()` | pladdrr | Centroid, spread, skewness, kurtosis |
| `trk_mfcc()` | C++ SPTK | Mel-frequency cepstral coefficients |

---

## 4. Energy & Amplitude (4 functions)

| Function | Engine | Description |
|----------|--------|-------------|
| `trk_rms()` | C ASSP | Root Mean Square amplitude |
| `trk_zcr()` | C ASSP | Zero-crossing rate |
| `trk_acf()` | C ASSP | Autocorrelation function |
| `trk_intensity()` | pladdrr | Praat perceived loudness |

---

## 5. Voice Quality & Aperiodicity (15 functions)

### Track functions (6)

| Function | Engine | Description |
|----------|--------|-------------|
| `trk_d4c()` | C++ WORLD | Band aperiodicity (high-quality) |
| `trk_brouhaha()` | Python DL | VAD + SNR + C50 (50-100x optimized) |
| `trk_creak_union()` | Python DL | Creaky voice detection (AM + CD NN) |
| `trk_cpps()` | pladdrr | Cepstral Peak Prominence Smoothed |
| `trk_vuv()` | pladdrr | Voiced/unvoiced detection |
| `trk_praatsauce()` | pladdrr | 36 VoiceSauce-style measures |

### Summary functions (9)

| Function | Engine | Measures | Description |
|----------|--------|----------|-------------|
| `lst_vat()` | Python | 132 | Voice Analysis Toolbox (jitter, shimmer, HNR, MFCC, wavelet) |
| `lst_voice_sauce()` | Python | 40+ | VoiceSauce parameters (F0, harmonics, CPP, spectral) |
| `lst_covarep_vq()` | Python | ~20 | COVAREP voice quality |
| `lst_vq()` | pladdrr | 36 | Voice quality summary (jitter, shimmer, HNR batched) |
| `lst_pharyngeal()` | pladdrr | 68 | Pharyngeal voice quality (most comprehensive) |
| `lst_voice_report()` | pladdrr | ~15 | Praat voice report |
| `lst_voice_tremor()` | pladdrr | ~10 | Voice tremor analysis |
| `lst_dsi()` | pladdrr | 1 | Dysphonia Severity Index |
| `lst_avqi()` | pladdrr | 1 | Acoustic Voice Quality Index |

---

## 6. Prosody & Intonation (2 functions)

| Function | Engine | Measures | Description |
|----------|--------|----------|-------------|
| `lst_dysprosody()` | Python | 193 | MOMEL-INTSINT, spectral tilt, statistics |
| `lst_voxit()` | Python | 11 | Speaking rate, pauses, rhythmic complexity, pitch dynamics |

---

## 7. Source-Filter Decomposition (3 functions)

| Function | Engine | Description |
|----------|--------|-------------|
| `trk_gfmiaif()` | Python | GFM-IAIF glottal flow model separation |
| `trk_covarep_iaif()` | Python | COVAREP iterative adaptive inverse filtering |
| `trk_excite()` | Python | Excitation signal extraction |

---

## 8. OpenSMILE Feature Sets (4 functions)

| Function | Features | Description |
|----------|----------|-------------|
| `lst_GeMAPS()` | 62 | Minimalistic standard (pitch, jitter, formants, loudness, HNR) |
| `lst_eGeMAPS()` | 88 | Extended (+ bandwidths, MFCCs, spectral flux) |
| `lst_emobase()` | ~1K | Emotional voice features |
| `lst_ComParE_2016()` | ~6K | Computational Paralinguistics challenge set |

All default to native C++ (`use_cpp = TRUE`); Python fallback available.

---

## 9. Phonological Classification (2 functions)

| Function | Engine | Description |
|----------|--------|-------------|
| `trk_phonet()` | Python DL | Phonological posterior probabilities (track) |
| `lst_phonet()` | Python DL | Phonological posterior probabilities (summary) |

---

## 10. Epoch Detection (2 functions)

| Function | Engine | Description |
|----------|--------|-------------|
| `trk_pitchmark_estk()` | C++ ESTK | Glottal closure instants (laryngograph) |
| `trk_pitchmark_reaper()` | C++ SPTK | REAPER pitch marks |

---

## 11. Unit Conversion (14 functions)

All prefixed `ucnv_*`. Psychoacoustic scale conversions:

| Function | Conversion |
|----------|------------|
| `ucnv_hz_to_bark()` / `ucnv_bark_to_hz()` | Hz <-> Bark |
| `ucnv_hz_to_erb()` / `ucnv_erb_to_hz()` | Hz <-> ERB |
| `ucnv_hz_to_mel()` / `ucnv_mel_to_hz()` | Hz <-> Mel |
| `ucnv_hz_to_semitone()` / `ucnv_semitone_to_hz()` | Hz <-> Semitone |
| `ucnv_db_and_hz_to_phon()` / `ucnv_phon_and_hz_to_db()` | dB+Hz <-> Phon |
| `ucnv_db_and_hz_to_sone()` / `ucnv_sone_and_hz_to_db()` | dB+Hz <-> Sone |
| `ucnv_phon_to_sone()` / `ucnv_sone_to_phon()` | Phon <-> Sone |

---

## 12. I/O (6 functions)

| Function | Description |
|----------|-------------|
| `read_audio()` | Read any audio/video into AVAudio |
| `read_ssff()` | Read SSFF (Simple Signal File Format) |
| `read_jstf()` | Read JSON Track Format (JSTF) |
| `read_track()` | Unified reader (auto-detects SSFF or JSTF) |
| `write_ssff()` | Write SSFF files |
| `write_jstf()` | Write JSTF files |

---

## 13. Classes & Generics (7 exports)

| Export | Type | Description |
|--------|------|-------------|
| `AVAudio` | S7 class | In-memory audio with preprocessing |
| `is.AsspDataObj()` | Predicate | Type check for AsspDataObj |
| `dur()` | S3 generic | Duration in seconds |
| `numRecs()` | S3 generic | Number of records/frames |
| `rate()` | S3 generic | Sample rate |
| `startTime()` | S3 generic | Start time |
| `tracks()` | S3 generic | Track names |

---

## 14. Deprecated (1 function)

| Function | Replacement |
|----------|-------------|
| `trk_egg_f0() (eggstract)` | `trk_egg_f0()` |

---

## Performance Tiers

| Tier | Engine | Speed | Examples |
|------|--------|-------|---------|
| 1 | C++ (SPTK/WORLD/ESTK) | Fastest | `trk_pitch_rapt`, `trk_pitch_dio`, `trk_mfcc`, `trk_d4c` |
| 2 | C (ASSP) | Fast | `trk_formant_forest`, `trk_dft_spectrum`, `trk_rms` |
| 3 | pladdrr (R/C++) | Fast | `trk_pitch_cc`, `trk_formant_burg`, `trk_cpps`, `lst_vq` |
| 4 | Python DL | Moderate | `trk_swiftf0`, `trk_pitch_crepe`, `trk_brouhaha` |
| 5 | Python classical | Moderate | `trk_pitch_yin`, `trk_pitch_pyin`, `trk_gfmiaif` |

---

## Export Summary

| Category | Count |
|----------|-------|
| `trk_*` | 49 |
| `lst_*` | 17 |
| `ucnv_*` | 14 |
| `read_*` | 4 |
| `write_*` | 2 |
| Classes & generics | 7 |
| **Total exports** | **93** |
| S3method registrations | 16 |

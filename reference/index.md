# Package index

## Pitch & F0 Tracking

Fundamental frequency estimation — 18 algorithms spanning C++ (fastest),
classical signal processing, and deep learning.

- [`trk_pitch_rapt()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_rapt.md)
  : Track fundamental frequency using RAPT (Robust Algorithm for Pitch
  Tracking)
- [`trk_pitch_swipe()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_swipe.md)
  : Track fundamental frequency using SWIPE (Sawtooth Waveform Inspired
  Pitch Estimator)
- [`trk_pitch_dio()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_dio.md)
  : DIO Pitch Tracking (C++ implementation)
- [`trk_pitch_harvest()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_harvest.md)
  : Harvest Pitch Tracking (C++ implementation)
- [`trk_pitch_reaper()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_reaper.md)
  : Track fundamental frequency using REAPER (Robust Epoch And Pitch
  EstimatoR)
- [`trk_pitch_yin()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_yin.md)
  : Track fundamental frequency using the YIN algorithm
- [`trk_pitch_pyin()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_pyin.md)
  : Track fundamental frequency using probabilistic YIN (pYIN)
- [`trk_pitch_pda()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_pda.md)
  : Track fundamental frequency using the ESTk PDA algorithm
- [`trk_tandem()`](https://humlab-speech.github.io/superassp/reference/trk_tandem.md)
  : Track pitch and voiced speech using the TANDEM-STRAIGHT algorithm
- [`trk_pitch_swiftf0()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_swiftf0.md)
  : Track fundamental frequency using SwiftF0 (ONNX)
- [`trk_pitch_crepe()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_crepe.md)
  : Track fundamental frequency and periodicity using CREPE (ONNX)
- [`trk_pitch_cc()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_cc.md)
  : Pitch tracking via Praat's cross-correlation method
- [`trk_pitch_ac()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_ac.md)
  : Pitch tracking via Praat's autocorrelation method
- [`trk_pitch_shs()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_shs.md)
  : Pitch tracking via Praat's subharmonic summation (SHS) method
- [`trk_pitch_spinet()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_spinet.md)
  : Pitch tracking via Praat's SPINET method
- [`trk_pitch_srh()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_srh.md)
  : Track fundamental frequency using the Summation of Residual
  Harmonics (SRH)
- [`trk_pitch_ksv()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_ksv.md)
  : Track fundamental frequency using the KSV periodicity detector
- [`trk_pitch_mhs()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_mhs.md)
  : Track pitch using the Modified Harmonic Sieve algorithm
- [`trk_pitch_snack()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_snack.md)
  : Track fundamental frequency using the Snack/ESPS dp_f0 algorithm

## Formant Analysis

Formant frequency and bandwidth tracking.

- [`trk_formant_forest()`](https://humlab-speech.github.io/superassp/reference/trk_formant_forest.md)
  : Track formant frequencies and bandwidths (FOREST)
- [`trk_formant_deepformants()`](https://humlab-speech.github.io/superassp/reference/trk_formant_deepformants.md)
  : Track formant frequencies using DeepFormants (ONNX)
- [`trk_formant_tvwlp()`](https://humlab-speech.github.io/superassp/reference/trk_formant_tvwlp.md)
  : Track formants using Time-Varying Weighted Linear Prediction (TVWLP)
- [`trk_formant_burg()`](https://humlab-speech.github.io/superassp/reference/trk_formant_burg.md)
  : Formant frequencies and bandwidths via Praat's Burg method
- [`trk_formant_cgdzp()`](https://humlab-speech.github.io/superassp/reference/trk_formant_cgdzp.md)
  : Track formants using Chirp Group Delay Zero-Phase (CGDZP) analysis
- [`trk_formant_snack()`](https://humlab-speech.github.io/superassp/reference/trk_formant_snack.md)
  : Track formants and bandwidths using the Snack/ESPS LPC tracker
- [`trk_formant_formantnet()`](https://humlab-speech.github.io/superassp/reference/trk_formant_formantnet.md)
  : Track formant frequencies and bandwidths using FormantNet (ONNX)

## Spectral Analysis

Spectrum estimation, cepstral analysis, and spectral moments.

- [`trk_dft_spectrum()`](https://humlab-speech.github.io/superassp/reference/trk_dft_spectrum.md)
  : Track short-term DFT power spectrum
- [`trk_css_spectrum()`](https://humlab-speech.github.io/superassp/reference/trk_css_spectrum.md)
  : Track cepstrally-smoothed spectrum
- [`trk_lps_spectrum()`](https://humlab-speech.github.io/superassp/reference/trk_lps_spectrum.md)
  : Track LP-smoothed spectrum
- [`trk_cepstrum()`](https://humlab-speech.github.io/superassp/reference/trk_cepstrum.md)
  : Track short-term cepstral coefficients
- [`trk_spectral_moments()`](https://humlab-speech.github.io/superassp/reference/trk_spectral_moments.md)
  : Spectral moments (CoG, SD, skewness, kurtosis)
- [`trk_mfcc()`](https://humlab-speech.github.io/superassp/reference/trk_mfcc.md)
  : Extract Mel-Frequency Cepstral Coefficients (MFCCs) via SPTK
- [`trk_cheap_trick()`](https://humlab-speech.github.io/superassp/reference/trk_cheap_trick.md)
  : CheapTrick Spectral Envelope Estimation (WORLD vocoder, C++
  implementation)

## Energy & Amplitude

Signal energy, zero-crossings, autocorrelation, and intensity.

- [`trk_rms()`](https://humlab-speech.github.io/superassp/reference/trk_rms.md)
  : Track short-term RMS amplitude
- [`trk_zcr()`](https://humlab-speech.github.io/superassp/reference/trk_zcr.md)
  : Track short-term zero-crossing rate
- [`trk_acf()`](https://humlab-speech.github.io/superassp/reference/trk_acf.md)
  : Track short-term autocorrelation function
- [`trk_intensity()`](https://humlab-speech.github.io/superassp/reference/trk_intensity.md)
  : Sound intensity contour

## Voice Quality & Aperiodicity

Comprehensive voice quality assessment — from single-track measures to
132-parameter dysphonia toolboxes.

- [`trk_d4c()`](https://humlab-speech.github.io/superassp/reference/trk_d4c.md)
  : Estimate band aperiodicity using the D4C algorithm (WORLD vocoder)
- [`trk_cpps()`](https://humlab-speech.github.io/superassp/reference/trk_cpps.md)
  : Cepstral Peak Prominence Smoothed (CPPS)
- [`trk_vuv()`](https://humlab-speech.github.io/superassp/reference/trk_vuv.md)
  : Voiced/unvoiced segmentation via two-pass adaptive pitch detection
- [`trk_praatsauce()`](https://humlab-speech.github.io/superassp/reference/trk_praatsauce.md)
  : Comprehensive voice quality feature set via PraatSauce
- [`lst_covarep_vq()`](https://humlab-speech.github.io/superassp/reference/lst_covarep_vq.md)
  : Extract voice quality parameters (NAQ, QOQ, H1-H2, HRF, PSP) per
  utterance
- [`lst_vq()`](https://humlab-speech.github.io/superassp/reference/lst_vq.md)
  : Voice Quality Measurements using pladdrr
- [`lst_pharyngeal()`](https://humlab-speech.github.io/superassp/reference/lst_pharyngeal.md)
  : Pharyngeal Voice Quality Analysis
- [`lst_voice_report()`](https://humlab-speech.github.io/superassp/reference/lst_voice_report.md)
  : Voice Report Analysis (pladdrr)
- [`lst_voice_tremor()`](https://humlab-speech.github.io/superassp/reference/lst_voice_tremor.md)
  : Vocal Tremor Analysis Using pladdrr
- [`lst_dsi()`](https://humlab-speech.github.io/superassp/reference/lst_dsi.md)
  : Dysphonia Severity Index (DSI) Analysis (pladdrr)
- [`lst_avqi()`](https://humlab-speech.github.io/superassp/reference/lst_avqi.md)
  : Acoustic Voice Quality Index (AVQI) using pladdrr
- [`trk_covarep_creak()`](https://humlab-speech.github.io/superassp/reference/trk_covarep_creak.md)
  : Detect creaky voice (vocal fry) per frame
- [`trk_covarep_env_te()`](https://humlab-speech.github.io/superassp/reference/trk_covarep_env_te.md)
  : Estimate spectral envelope using the True Envelope (Teager energy)
  method
- [`trk_covarep_vad_drugman()`](https://humlab-speech.github.io/superassp/reference/trk_covarep_vad_drugman.md)
  : Detect voiced frames using Drugman's multi-branch VAD
- [`trk_covarep_vq_gci()`](https://humlab-speech.github.io/superassp/reference/trk_covarep_vq_gci.md)
  : Track GCI-anchored voice quality measures as a time series
- [`trk_creak_vat()`](https://humlab-speech.github.io/superassp/reference/trk_creak_vat.md)
  : Detect creaky voice using the Kane-Drugman VAT creak detector
- [`trk_gci_vat()`](https://humlab-speech.github.io/superassp/reference/trk_gci_vat.md)
  : Detect glottal closure instants (GCIs) using SE-VQ via voiceanalysis
- [`trk_iaif_vat()`](https://humlab-speech.github.io/superassp/reference/trk_iaif_vat.md)
  : Estimate glottal flow via IAIF using voiceanalysis
- [`trk_mdq_vat()`](https://humlab-speech.github.io/superassp/reference/trk_mdq_vat.md)
  : Track Maxima Dispersion Quotient (MDQ) for breathy/tense voice
  discrimination
- [`trk_peakslope_vat()`](https://humlab-speech.github.io/superassp/reference/trk_peakslope_vat.md)
  : Peak slope via voiceanalysis Daless wavelet bank
- [`trk_peakslope()`](https://humlab-speech.github.io/superassp/reference/trk_peakslope.md)
  : Track spectral tilt using D'Alessandro PeakSlope (Morlet wavelet)
- [`trk_pitch_vat()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_vat.md)
  : Track fundamental frequency using SRH via the voiceanalysis package
- [`trk_hmpd()`](https://humlab-speech.github.io/superassp/reference/trk_hmpd.md)
  : Extract harmonic model phase distortion features (HMPD): AE, PDM,
  PDD
- [`lst_lf_vat_synthesis()`](https://humlab-speech.github.io/superassp/reference/lst_lf_vat_synthesis.md)
  : Synthesise an LF model glottal pulse via voiceanalysis
- [`lst_vq_vat()`](https://humlab-speech.github.io/superassp/reference/lst_vq_vat.md)
  : Per-GCI voice-quality summary via voiceanalysis
- [`lst_polarity()`](https://humlab-speech.github.io/superassp/reference/lst_polarity.md)
  : Signal polarity detection (RESKEW algorithm)

## Prosody & Intonation

Prosodic features, rhythm, articulation complexity, and pitch modelling.

- [`lst_dysprosody()`](https://humlab-speech.github.io/superassp/reference/lst_dysprosody.md)
  : Extract Dysprosody Prosodic Features
- [`lst_voxit()`](https://humlab-speech.github.io/superassp/reference/lst_voxit.md)
  : Extract Voxit prosodic complexity features from audio files
- [`lst_vowel_space()`](https://humlab-speech.github.io/superassp/reference/lst_vowel_space.md)
  : Vowel space analysis (F1×F2 area ratio)
- [`momel()`](https://humlab-speech.github.io/superassp/reference/momel.md)
  : Run MOMEL algorithm on F0 values
- [`intsint()`](https://humlab-speech.github.io/superassp/reference/intsint.md)
  : Run INTSINT algorithm on MOMEL targets
- [`prosody_measures()`](https://humlab-speech.github.io/superassp/reference/prosody_measures.md)
  : Compute prosodic measures from audio file or Sound object
- [`voxit-analysis-stats`](https://humlab-speech.github.io/superassp/reference/voxit-analysis-stats.md)
  : Voxit analysis statistics
- [`voxit-dsp-utils`](https://humlab-speech.github.io/superassp/reference/voxit-dsp-utils.md)
  : DSP utility functions from Voxit

## Source-Filter Decomposition

Glottal source and vocal tract separation.

- [`trk_gfmiaif()`](https://humlab-speech.github.io/superassp/reference/trk_gfmiaif.md)
  : Decompose speech into vocal tract, glottis, and lip radiation LP
  filters (GFM-IAIF)
- [`trk_covarep_iaif()`](https://humlab-speech.github.io/superassp/reference/trk_covarep_iaif.md)
  : Extract glottal flow waveform using Iterative Adaptive Inverse
  Filtering (IAIF)

## OpenSMILE Feature Sets

Standardized acoustic feature extraction via OpenSMILE C++.

- [`lst_GeMAPS()`](https://humlab-speech.github.io/superassp/reference/lst_GeMAPS.md)
  : Compute the GeMAPS openSMILE feature set (C++ Implementation)
- [`lst_eGeMAPS()`](https://humlab-speech.github.io/superassp/reference/lst_eGeMAPS.md)
  : Compute the eGeMAPS openSMILE feature set
- [`lst_emobase()`](https://humlab-speech.github.io/superassp/reference/lst_emobase.md)
  : Compute the emobase openSMILE feature set
- [`lst_ComParE_2016()`](https://humlab-speech.github.io/superassp/reference/lst_ComParE_2016.md)
  : Compute the ComParE 2016 openSMILE feature set

## Epoch Detection

Glottal closure instants and pitch marks.

- [`trk_pitchmark_estk()`](https://humlab-speech.github.io/superassp/reference/trk_pitchmark_estk.md)
  : Detect glottal closure instants in laryngograph signals using ESTk
  pitchmark
- [`trk_pitchmark_reaper()`](https://humlab-speech.github.io/superassp/reference/trk_pitchmark_reaper.md)
  : Detect glottal closure instants using REAPER (pitch marks)
- [`lst_covarep_gci_sedreams()`](https://humlab-speech.github.io/superassp/reference/lst_covarep_gci_sedreams.md)
  : SEDREAMS Glottal Closure Instant Detection

## Psychoacoustics

Equal-loudness contours and loudness unit conversions (ISO 226, ISO
532).

- [`iso226_phon`](https://humlab-speech.github.io/superassp/reference/iso226_phon.md)
  : ISO 226:2023 Phon (Loudness Level) Conversions
- [`iso532-sone`](https://humlab-speech.github.io/superassp/reference/iso532-sone.md)
  : ISO 532 Sone (Loudness) Conversions

## Unit Conversion

Psychoacoustic scale conversions (Hz, Bark, ERB, Mel, semitone, phon,
sone).

- [`ucnv_bark_to_hz()`](https://humlab-speech.github.io/superassp/reference/ucnv_bark_to_hz.md)
  : Convert Bark Scale to Frequency
- [`ucnv_db_and_hz_to_phon()`](https://humlab-speech.github.io/superassp/reference/ucnv_db_and_hz_to_phon.md)
  : Convert Sound Pressure Level and Frequency to Loudness Level (Phon)
- [`ucnv_db_and_hz_to_sone()`](https://humlab-speech.github.io/superassp/reference/ucnv_db_and_hz_to_sone.md)
  : Convert dB and Hz Directly to Sone
- [`ucnv_erb_to_hz()`](https://humlab-speech.github.io/superassp/reference/ucnv_erb_to_hz.md)
  : Convert ERB-rate Scale to Frequency
- [`ucnv_hz_to_bark()`](https://humlab-speech.github.io/superassp/reference/ucnv_hz_to_bark.md)
  : Convert Frequency to Bark Scale
- [`ucnv_hz_to_erb()`](https://humlab-speech.github.io/superassp/reference/ucnv_hz_to_erb.md)
  : Convert Frequency to ERB-rate Scale
- [`ucnv_hz_to_mel()`](https://humlab-speech.github.io/superassp/reference/ucnv_hz_to_mel.md)
  : Convert Frequency to Mel Scale
- [`ucnv_hz_to_semitone()`](https://humlab-speech.github.io/superassp/reference/ucnv_hz_to_semitone.md)
  : Convert Frequency to Semitones
- [`ucnv_mel_to_hz()`](https://humlab-speech.github.io/superassp/reference/ucnv_mel_to_hz.md)
  : Convert Mel Scale to Frequency
- [`ucnv_phon_and_hz_to_db()`](https://humlab-speech.github.io/superassp/reference/ucnv_phon_and_hz_to_db.md)
  : Convert Loudness Level (Phon) and Frequency to Sound Pressure Level
- [`ucnv_phon_to_sone()`](https://humlab-speech.github.io/superassp/reference/ucnv_phon_to_sone.md)
  : Convert Phon to Sone
- [`ucnv_semitone_to_hz()`](https://humlab-speech.github.io/superassp/reference/ucnv_semitone_to_hz.md)
  : Convert Semitones to Frequency
- [`ucnv_sone_and_hz_to_db()`](https://humlab-speech.github.io/superassp/reference/ucnv_sone_and_hz_to_db.md)
  : Convert Sone and Hz to dB
- [`ucnv_sone_to_phon()`](https://humlab-speech.github.io/superassp/reference/ucnv_sone_to_phon.md)
  : Convert Sone to Phon

## I/O — Audio & SSFF

Load audio, read/write SSFF signal files.

- [`read_audio()`](https://humlab-speech.github.io/superassp/reference/read_audio.md)
  : Read an audio file into an AsspDataObj
- [`av_to_asspDataObj()`](https://humlab-speech.github.io/superassp/reference/av_to_asspDataObj.md)
  : Convert audio file to AsspDataObj
- [`avaudio_to_av()`](https://humlab-speech.github.io/superassp/reference/avaudio_to_av.md)
  : Convert AVAudio to av::read_audio_bin Format
- [`avaudio_to_tempfile()`](https://humlab-speech.github.io/superassp/reference/avaudio_to_tempfile.md)
  : Convert AVAudio to Temporary WAV File
- [`read_ssff()`](https://humlab-speech.github.io/superassp/reference/read_ssff.md)
  : Read an SSFF or audio file into an AsspDataObj
- [`write_ssff()`](https://humlab-speech.github.io/superassp/reference/write_ssff.md)
  : Write an AsspDataObj to an SSFF file
- [`prep_recode()`](https://humlab-speech.github.io/superassp/reference/prep_recode.md)
  : Re-encode Media File with Custom Parameters

## I/O — JSON Track Format (JSTF)

Create, read, write, and manipulate JSON Track Format files.

- [`read_track()`](https://humlab-speech.github.io/superassp/reference/read_track.md)
  : Unified Track Reading Interface
- [`write_track()`](https://humlab-speech.github.io/superassp/reference/write_track.md)
  : Write Track to File
- [`create_json_track_obj()`](https://humlab-speech.github.io/superassp/reference/create_json_track_obj.md)
  : Create a JsonTrackObj
- [`append_json_track_slice()`](https://humlab-speech.github.io/superassp/reference/append_json_track_slice.md)
  : Append a slice to JsonTrackObj
- [`merge_json_tracks()`](https://humlab-speech.github.io/superassp/reference/merge_json_tracks.md)
  : Merge multiple JsonTrackObj files
- [`subset_json_track()`](https://humlab-speech.github.io/superassp/reference/subset_json_track.md)
  : Subset JsonTrackObj
- [`validate_json_track()`](https://humlab-speech.github.io/superassp/reference/validate_json_track.md)
  : Validate JsonTrackObj
- [`store_slice()`](https://humlab-speech.github.io/superassp/reference/store_slice.md)
  : Provides the ability to store a multidimensional feature set related
  to a part of a signal.
- [`json_track_core`](https://humlab-speech.github.io/superassp/reference/json_track_core.md)
  : JSON Track Object Core Functions
- [`json_track_methods`](https://humlab-speech.github.io/superassp/reference/json_track_methods.md)
  : JSON Track Conversion Methods
- [`read_jstf()`](https://humlab-speech.github.io/superassp/reference/read_jstf.md)
  : Read JSTF File
- [`write_jstf()`](https://humlab-speech.github.io/superassp/reference/write_jstf.md)
  : Write JSTF Object to File
- [`jstf_io`](https://humlab-speech.github.io/superassp/reference/jstf_io.md)
  : JSTF (JSON Sparse Track Format) I/O Functions

## Classes & Data Structures

Core data classes, S7 generics, and AsspDataObj accessors.

- [`is_avaudio()`](https://humlab-speech.github.io/superassp/reference/is_avaudio.md)
  : Check if Object is AVAudio
- [`as_avaudio()`](https://humlab-speech.github.io/superassp/reference/as_avaudio.md)
  : Convert to AVAudio Object
- [`as.data.frame(`*`<AsspDataObj>`*`)`](https://humlab-speech.github.io/superassp/reference/AsspDataObj.md)
  [`print(`*`<AsspDataObj>`*`)`](https://humlab-speech.github.io/superassp/reference/AsspDataObj.md)
  [`as_tibble(`*`<AsspDataObj>`*`)`](https://humlab-speech.github.io/superassp/reference/AsspDataObj.md)
  [`cut(`*`<AsspDataObj>`*`)`](https://humlab-speech.github.io/superassp/reference/AsspDataObj.md)
  : AsspDataObj — ASSP Data Object
- [`print(`*`<JsonTrackObj>`*`)`](https://humlab-speech.github.io/superassp/reference/JsonTrackObj.md)
  [`as.data.frame(`*`<JsonTrackObj>`*`)`](https://humlab-speech.github.io/superassp/reference/JsonTrackObj.md)
  [`as_tibble(`*`<JsonTrackObj>`*`)`](https://humlab-speech.github.io/superassp/reference/JsonTrackObj.md)
  [`summary(`*`<JsonTrackObj>`*`)`](https://humlab-speech.github.io/superassp/reference/JsonTrackObj.md)
  : JsonTrackObj — JSON Track Format Object
- [`s7-methods`](https://humlab-speech.github.io/superassp/reference/s7-methods.md)
  : S7 Method System for DSP Functions
- [`sample_rate()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md)
  [`n_records()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md)
  [`signal_duration()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md)
  [`start_time()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md)
  [`track_names()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md)
  [`file_path()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md)
  [`track_formats()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md)
  [`dur()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md)
  [`numRecs()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md)
  [`rate()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md)
  [`startTime()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md)
  [`tracks()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md)
  : Accessor methods for AsspDataObj and JsonTrackObj

## ASSP Constants & Types

Enumerated types and constants from the ASSP signal processing library.

- [`AsspFileFormats`](https://humlab-speech.github.io/superassp/reference/AsspFileFormats.md)
  : list of possibly useful file formats for AsspDataObj corresponding
  to the first element of the fileInfo attribute
- [`AsspLpTypes()`](https://humlab-speech.github.io/superassp/reference/AsspLpTypes.md)
  : AsspLpTypes
- [`AsspSpectTypes()`](https://humlab-speech.github.io/superassp/reference/AsspSpectTypes.md)
  : AsspSpectTypes
- [`AsspWindowTypes()`](https://humlab-speech.github.io/superassp/reference/AsspWindowTypes.md)
  : AsspWindowTypes
- [`isAsspLpType()`](https://humlab-speech.github.io/superassp/reference/isAsspLpType.md)
  : isAsspLpType
- [`isAsspSpectType()`](https://humlab-speech.github.io/superassp/reference/isAsspSpectType.md)
  : isAsspSpectType
- [`isAsspWindowType()`](https://humlab-speech.github.io/superassp/reference/isAsspWindowType.md)
  : isAsspWindowType
- [`wrasspOutputInfos`](https://humlab-speech.github.io/superassp/reference/wrasspOutputInfos.md)
  : list of default output extensions, track names and output type for
  each signal processing function in wrassp

## Package Utilities

Introspection, visualisation, and media conversion helpers.

- [`get_track_label()`](https://humlab-speech.github.io/superassp/reference/get_track_label.md)
  : Get track label for plotting
- [`get_track_label_expr()`](https://humlab-speech.github.io/superassp/reference/get_track_label_expr.md)
  : Get track label as expression for plotting
- [`ggtrack()`](https://humlab-speech.github.io/superassp/reference/ggtrack.md)
  : Create ggplot with automatic track labels
- [`differentiate()`](https://humlab-speech.github.io/superassp/reference/differentiate.md)
  : Derivation of SSFF track objects

## Python / pladdrr Integration Helpers

Audio loading and format conversion helpers for pladdrr backends.

- [`av_load_for_pladdrr()`](https://humlab-speech.github.io/superassp/reference/av_load_for_pladdrr.md)
  : Load audio file as pladdrr Sound object
- [`pladdrr_df_to_superassp()`](https://humlab-speech.github.io/superassp/reference/pladdrr_df_to_superassp.md)
  : Convert pladdrr data frame to superassp format
- [`rmsana_memory()`](https://humlab-speech.github.io/superassp/reference/rmsana_memory.md)
  : Perform RMS analysis on AsspDataObj in memory
- [`process_media_file()`](https://humlab-speech.github.io/superassp/reference/process_media_file.md)
  : Process audio from any media file format

## Legacy Functions

Pre-v1.0 functions retained for backward compatibility. Prefer the
`trk_*` equivalents for new code.

- [`trk_pitch_mhs()`](https://humlab-speech.github.io/superassp/reference/trk_pitch_mhs.md)
  : Track pitch using the Modified Harmonic Sieve algorithm
- [`harmonics()`](https://humlab-speech.github.io/superassp/reference/harmonics.md)
  : Compute the harmonic frequency structure from f0 measurements
- [`trk_afdiff()`](https://humlab-speech.github.io/superassp/reference/trk_afdiff.md)
  : Differentiate an audio waveform
- [`trk_affilter()`](https://humlab-speech.github.io/superassp/reference/trk_affilter.md)
  : Apply a digital filter to audio signals
- [`trk_arf()`](https://humlab-speech.github.io/superassp/reference/trk_arf.md)
  : Track LP-derived vocal tract area function coefficients
- [`trk_lar()`](https://humlab-speech.github.io/superassp/reference/trk_lar.md)
  : Track LP-derived log area ratios
- [`trk_lpc()`](https://humlab-speech.github.io/superassp/reference/trk_lpc.md)
  : Track LP filter coefficients
- [`trk_rfc()`](https://humlab-speech.github.io/superassp/reference/trk_rfc.md)
  : Track LP reflection coefficients
- [`useWrasspLogger`](https://humlab-speech.github.io/superassp/reference/useWrasspLogger.md)
  : package variable to force the usage of the logger set to FALSE by
  default
- [`read.AsspDataObj()`](https://humlab-speech.github.io/superassp/reference/read.AsspDataObj.md)
  : read.AsspDataObj from a signal/parameter file
- [`write.AsspDataObj()`](https://humlab-speech.github.io/superassp/reference/write.AsspDataObj.md)
  : write.AsspDataObj to file

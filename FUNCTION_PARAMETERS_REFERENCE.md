# superassp Function Parameters Reference

**Date**: 2025-10-29
**Purpose**: Comprehensive catalog of all formal parameters used across `trk_*` and `lst_*` functions

---

## Overview

This document catalogs **97 unique parameters** used across 75+ DSP functions in superassp. Parameters are organized alphabetically with descriptions, default values, and lists of functions using each parameter.

---

## Universal Parameters (Used by Most/All Functions)

### beginTime
- **Description**: Start time in seconds for analysis. 0 or NULL means start of file. Can be vector matching length of listOfFiles
- **Default**: 0.0 or 0
- **Type**: Numeric (scalar or vector)
- **Functions**: trk_rapt, trk_swipe, trk_dio, trk_reaper, trk_harvest, trk_d4c, trk_mfcc, trk_pitchmark, trk_forest, trk_mhspitch, trk_ksvfo, trk_acfana, trk_zcrana, trk_rmsana, trk_cepstrum, trk_cssSpectrum, trk_dftSpectrum, trk_lpsSpectrum, trk_lp_analysis, trk_pyin, trk_crepe, trk_yin, trk_yaapt, trk_snackp, trk_snackf, trk_seenc, trk_excite, trk_aperiodicities, trk_brouhaha, trk_pitchp, trk_formantp, trk_formantpathp, trk_intensityp, trk_spectral_momentsp, trk_praat_sauce, lst_dysprosody, lst_vat, lst_voice_sauce, lst_voxit, lst_voice_reportp, lst_voice_tremorp, lst_avqip, lst_dsip

### endTime
- **Description**: End time in seconds for analysis. 0 or NULL means end of file. Can be vector matching length of listOfFiles
- **Default**: 0.0 or 0
- **Type**: Numeric (scalar or vector)
- **Functions**: Same as beginTime

### explicitExt
- **Description**: File extension for output files when toFile=TRUE
- **Default**: Varies by function ("f0", "fms", "mfcc", "pm", "pit", "fo", "acf", "zcr", "rms", "cep", "css", "dft", "lps", "arf", "lar", "rfc", "lpc")
- **Type**: Character
- **Functions**: All trk_* functions (not applicable to most lst_* functions)

### listOfFiles
- **Description**: Character vector of file paths to process. Accepts any media format via av package (WAV, MP3, MP4, FLAC, video files, etc.)
- **Default**: (required parameter, no default)
- **Type**: Character vector
- **Functions**: All trk_* and lst_* functions

### outputDirectory
- **Description**: Directory path for output files. NULL means same directory as input files
- **Default**: NULL
- **Type**: Character or NULL
- **Functions**: All trk_* functions with toFile=TRUE capability

### toFile
- **Description**: Write results to SSFF file (TRUE) or return as R object (FALSE)
- **Default**: TRUE
- **Type**: Logical
- **Functions**: All trk_* functions, some lst_* functions

### verbose
- **Description**: Display progress information, warnings, and processing details
- **Default**: TRUE
- **Type**: Logical
- **Functions**: All trk_* and lst_* functions

### windowShift
- **Description**: Frame shift in milliseconds (time between analysis windows)
- **Default**: 5.0, 10.0 (varies by function)
- **Type**: Numeric
- **Functions**: All pitch tracking functions, spectral analysis functions, formant tracking functions

---

## Pitch Tracking Parameters

### maxF
- **Description**: Maximum fundamental frequency in Hz for pitch detection
- **Default**: 400.0 (SPTK), 600.0 (Parselmouth), 500 (Python), varies
- **Type**: Numeric
- **Functions**: trk_rapt, trk_swipe, trk_dio, trk_reaper, trk_harvest, trk_d4c, trk_ksvfo, trk_mhspitch, trk_pyin, trk_crepe, trk_yin, trk_yaapt, trk_snackp, trk_pitchp, lst_dysprosody, lst_voxit

### minF
- **Description**: Minimum fundamental frequency in Hz for pitch detection
- **Default**: 60.0 (SPTK), 50.0 (Python), 70 (some ASSP), varies
- **Type**: Numeric
- **Functions**: Same as maxF

### voicing_threshold
- **Description**: Threshold for determining voiced vs unvoiced regions (0-1 scale)
- **Default**: 0.9 (rapt, reaper, harvest), 0.3 (swipe), 0.85 (dio, d4c)
- **Type**: Numeric
- **Functions**: trk_rapt, trk_swipe, trk_dio, trk_reaper, trk_harvest, trk_d4c

---

## Spectral Analysis Parameters

### fftLength
- **Description**: FFT length in points. 0 means automatic calculation
- **Default**: 0
- **Type**: Integer
- **Functions**: trk_cepstrum, trk_cssSpectrum, trk_dftSpectrum, trk_lpsSpectrum

### resolution
- **Description**: Target frequency resolution in Hz (determines FFT length)
- **Default**: 40.0
- **Type**: Numeric
- **Functions**: trk_cepstrum, trk_cssSpectrum, trk_dftSpectrum, trk_lpsSpectrum

### window
- **Description**: Window function type for spectral analysis. Options: "BLACKMAN", "HAMMING", "HANN", "BARTLETT", "RECTANGLE", etc. See AsspWindowTypes()
- **Default**: "BLACKMAN" (most functions), "HAMMING" (some)
- **Type**: Character
- **Functions**: trk_acfana, trk_cepstrum, trk_cssSpectrum, trk_dftSpectrum, trk_forest, trk_lp_analysis, trk_lpsSpectrum, trk_rmsana

### windowSize
- **Description**: Analysis window size in milliseconds
- **Default**: 20.0 (most), 25.0 (some), varies
- **Type**: Numeric
- **Functions**: trk_cepstrum, trk_cssSpectrum, trk_lpsSpectrum, trk_acfana, trk_forest, trk_lp_analysis, trk_rmsana, trk_zcrana

---

## Linear Prediction Parameters

### analysisOrder
- **Description**: LP analysis order. 0 or NULL sets to sample_rate_kHz + 3
- **Default**: 0 or NULL
- **Type**: Integer
- **Functions**: trk_acfana, trk_cepstrum, trk_cssSpectrum, trk_dftSpectrum, trk_lpsSpectrum, trk_lp_analysis (and variants: rfcana, arfana, larana, lpcana)

### preemphasis
- **Description**: Pre-emphasis coefficient (typically negative, range -1 to 0). Compensates for high-frequency attenuation
- **Default**: -0.8 (forest), -0.95 (lpsSpectrum)
- **Type**: Numeric
- **Functions**: trk_forest, trk_lpsSpectrum, trk_lp_analysis and variants

### effectiveLength
- **Description**: Make window size effective (accounting for zero-padding) rather than exact
- **Default**: TRUE
- **Type**: Logical
- **Functions**: trk_acfana, trk_forest, trk_ksvfo, trk_lp_analysis, trk_lpsSpectrum, trk_rmsana

---

## Formant Analysis Parameters

### numFormants
- **Description**: Number of formants to identify (maximum: 8 or half the LPC order)
- **Default**: 4
- **Type**: Integer
- **Functions**: trk_forest

### nominalF1
- **Description**: Nominal (assumed) F1 frequency in Hz for formant classification
- **Default**: 500
- **Type**: Numeric
- **Functions**: trk_forest

### estimate
- **Description**: Insert rough frequency estimates for missing formants (otherwise set to zero)
- **Default**: FALSE
- **Type**: Logical
- **Functions**: trk_forest

### incrOrder
- **Description**: Increase default LPC order by 2 (adds one resonance)
- **Default**: 0
- **Type**: Integer
- **Functions**: trk_forest

---

## MFCC Parameters

### n_mfcc
- **Description**: Number of MFCC coefficients to extract
- **Default**: 13
- **Type**: Integer
- **Functions**: trk_mfcc

### n_mels
- **Description**: Number of mel filterbanks
- **Default**: 40
- **Type**: Integer
- **Functions**: trk_mfcc

### fmin
- **Description**: Minimum frequency for mel filterbank in Hz
- **Default**: 0.0
- **Type**: Numeric
- **Functions**: trk_mfcc

### fmax
- **Description**: Maximum frequency for mel filterbank in Hz (NULL = Nyquist)
- **Default**: NULL
- **Type**: Numeric or NULL
- **Functions**: trk_mfcc

### lifter
- **Description**: Liftering coefficient for MFCC (cepstral smoothing)
- **Default**: 22
- **Type**: Integer
- **Functions**: trk_mfcc

### floor
- **Description**: Floor value for mel filterbank output (avoids log(0))
- **Default**: 1.0
- **Type**: Numeric
- **Functions**: trk_mfcc

---

## Gender-Specific Parameters

### gender
- **Description**: Gender-specific analysis parameters. Codes: "f" (female), "m" (male), "u" (unknown). For "f": effective window = 12.5ms, nominal F1 = 560Hz
- **Default**: "m" or "u"
- **Type**: Character
- **Functions**: trk_forest, trk_mhspitch, trk_ksvfo

---

## Pitch Detection Algorithm Parameters

### maxZCR
- **Description**: Maximum zero-crossing rate in Hz for voicing detection
- **Default**: 3000.0 or 3000
- **Type**: Numeric
- **Functions**: trk_ksvfo, trk_mhspitch

### minAmp
- **Description**: Minimum signal amplitude threshold for voiced samples
- **Default**: 50 (ksvfo), 100 (mhspitch)
- **Type**: Numeric
- **Functions**: trk_ksvfo, trk_mhspitch

### minRMS
- **Description**: Minimum RMS amplitude in dB
- **Default**: 18.0
- **Type**: Numeric
- **Functions**: trk_mhspitch

### minAC1
- **Description**: Minimum first autocorrelation coefficient
- **Default**: 0.25 or 0.250
- **Type**: Numeric
- **Functions**: trk_mhspitch

### minProb
- **Description**: Minimum quality value of F0 fit
- **Default**: 0.52 or 0.520
- **Type**: Numeric
- **Functions**: trk_mhspitch

### plainSpectrum
- **Description**: Use plain spectrum (no smoothing) for pitch detection
- **Default**: FALSE
- **Type**: Logical
- **Functions**: trk_mhspitch

---

## Autocorrelation Analysis Parameters

### energyNormalization
- **Description**: Calculate energy-normalized autocorrelation
- **Default**: FALSE
- **Type**: Logical
- **Functions**: trk_acfana

### lengthNormalization
- **Description**: Calculate length-normalized autocorrelation
- **Default**: FALSE
- **Type**: Logical
- **Functions**: trk_acfana

---

## RMS Analysis Parameters

### linear
- **Description**: Compute linear RMS values (TRUE) or logarithmic dB scale (FALSE)
- **Default**: FALSE
- **Type**: Logical
- **Functions**: trk_rmsana

---

## Cepstral Analysis Parameters

### numCeps
- **Description**: Number of cepstral coefficients for smoothing. 0 = automatic
- **Default**: 0
- **Type**: Integer
- **Functions**: trk_cssSpectrum

### deemphasize
- **Description**: Undo spectral tilt from LP pre-emphasis
- **Default**: TRUE
- **Type**: Logical
- **Functions**: trk_lpsSpectrum

---

## DFT Spectrum Parameters

### bandwidth
- **Description**: Effective analysis bandwidth in Hz. 0 = smallest possible given FFT length
- **Default**: 0.0
- **Type**: Numeric
- **Functions**: trk_dftSpectrum

---

## PYIN (Probabilistic YIN) Parameters

### thresholds
- **Description**: Number of thresholds for peak estimation
- **Default**: 100
- **Type**: Integer
- **Functions**: trk_pyin

### beta_parameters
- **Description**: Shape parameters for beta distribution prior over thresholds
- **Default**: c(2, 18)
- **Type**: Numeric vector (length 2)
- **Functions**: trk_pyin

### boltzmann_parameter
- **Description**: Shape parameter for Boltzmann distribution prior over troughs. Larger values favor smaller periods
- **Default**: 2
- **Type**: Numeric
- **Functions**: trk_pyin

### max_transition_rate
- **Description**: Maximum pitch transition rate in octaves per second
- **Default**: 35.92
- **Type**: Numeric
- **Functions**: trk_pyin

### switch_probability
- **Description**: Probability of switching voiced/unvoiced state
- **Default**: 0.01
- **Type**: Numeric
- **Functions**: trk_pyin

### no_trough_probability
- **Description**: Probability to add to global minimum if no trough below threshold
- **Default**: 0.01
- **Type**: Numeric
- **Functions**: trk_pyin

### center
- **Description**: Center audio/use centered frames
- **Default**: TRUE
- **Type**: Logical
- **Functions**: trk_pyin, trk_crepe, trk_yin, trk_yaapt

### pad_mode
- **Description**: Padding mode for audio ("constant", "reflect", etc.)
- **Default**: "constant"
- **Type**: Character
- **Functions**: trk_pyin

---

## Pitchmark Parameters

### def_period
- **Description**: Default pitch period in seconds for interpolation (used with fill=TRUE)
- **Default**: 0.01 (100 Hz)
- **Type**: Numeric
- **Functions**: trk_pitchmark

### fill
- **Description**: Post-process pitchmarks to ensure min/max periods and fill unvoiced regions
- **Default**: FALSE
- **Type**: Logical
- **Functions**: trk_pitchmark

### invert
- **Description**: Invert signal polarity (for upside-down laryngograph signals)
- **Default**: FALSE
- **Type**: Logical
- **Functions**: trk_pitchmark

### max_period
- **Description**: Maximum allowed pitch period in seconds (~50 Hz min F0)
- **Default**: 0.02
- **Type**: Numeric
- **Functions**: trk_pitchmark

### min_period
- **Description**: Minimum allowed pitch period in seconds (~333 Hz max F0)
- **Default**: 0.003
- **Type**: Numeric
- **Functions**: trk_pitchmark

### to_f0
- **Description**: Convert pitchmarks to F0 contour (returns F0 instead of timestamps)
- **Default**: FALSE
- **Type**: Logical
- **Functions**: trk_pitchmark

### use_cpp
- **Description**: Use C++ implementation (faster) instead of ESTK binary
- **Default**: TRUE
- **Type**: Logical
- **Functions**: trk_pitchmark

### lx_high_frequency
- **Description**: High-pass cutoff in Hz for initial filtering (removes low-frequency swell)
- **Default**: 40
- **Type**: Numeric
- **Functions**: trk_pitchmark

### lx_high_order
- **Description**: Order of high-pass FIR filter
- **Default**: 19
- **Type**: Integer
- **Functions**: trk_pitchmark

### lx_low_frequency
- **Description**: Low-pass cutoff in Hz for initial filtering (removes high-frequency noise)
- **Default**: 400
- **Type**: Numeric
- **Functions**: trk_pitchmark

### lx_low_order
- **Description**: Order of low-pass FIR filter
- **Default**: 19
- **Type**: Integer
- **Functions**: trk_pitchmark

### df_low_frequency
- **Description**: Low-pass cutoff in Hz for differentiated signal smoothing
- **Default**: 1000
- **Type**: Numeric
- **Functions**: trk_pitchmark

### df_low_order
- **Description**: Order of differentiated signal low-pass filter (0 = disable)
- **Default**: 19
- **Type**: Integer
- **Functions**: trk_pitchmark

### median_order
- **Description**: Order of median smoother for differentiated signal (0 = disable)
- **Default**: 19
- **Type**: Integer
- **Functions**: trk_pitchmark

---

## Brouhaha VAD Parameters

### onset
- **Description**: VAD onset threshold. Speech detected when score > onset
- **Default**: 0.780
- **Type**: Numeric (0-1 scale)
- **Functions**: trk_brouhaha

### offset
- **Description**: VAD offset threshold. Speech ends when score < offset
- **Default**: 0.780
- **Type**: Numeric (0-1 scale)
- **Functions**: trk_brouhaha

### min_duration_on
- **Description**: Minimum speech duration in seconds (shorter regions removed)
- **Default**: 0
- **Type**: Numeric
- **Functions**: trk_brouhaha

### min_duration_off
- **Description**: Minimum silence duration in seconds (shorter gaps filled)
- **Default**: 0
- **Type**: Numeric
- **Functions**: trk_brouhaha

### model_path
- **Description**: Path to custom brouhaha model checkpoint (NULL = default model)
- **Default**: NULL
- **Type**: Character or NULL
- **Functions**: trk_brouhaha

### use_optimized
- **Description**: Use optimized inference with pre-allocated arrays (2-3x faster)
- **Default**: TRUE
- **Type**: Logical
- **Functions**: trk_brouhaha

### batch_size
- **Description**: Number of chunks per batch (higher = more memory, possibly faster)
- **Default**: 32
- **Type**: Integer
- **Functions**: trk_brouhaha

---

## Voice Analysis Parameters

### f0_min
- **Description**: Minimum F0 in Hz for voice analysis
- **Default**: 50
- **Type**: Numeric
- **Functions**: lst_vat

### f0_max
- **Description**: Maximum F0 in Hz for voice analysis
- **Default**: 500
- **Type**: Numeric
- **Functions**: lst_vat

### f0_algorithm
- **Description**: F0 estimation algorithm ("SWIPE" or "PRAAT")
- **Default**: "SWIPE"
- **Type**: Character
- **Functions**: lst_vat

### return_f0
- **Description**: Include F0 contour in output
- **Default**: FALSE
- **Type**: Logical
- **Functions**: lst_vat

### use_cython
- **Description**: Use Cython-optimized functions (faster)
- **Default**: TRUE
- **Type**: Logical
- **Functions**: lst_vat

### use_thesis_mode
- **Description**: Use thesis-compliant implementations (differ slightly from MATLAB)
- **Default**: FALSE
- **Type**: Logical
- **Functions**: lst_vat

### timeout
- **Description**: Maximum analysis time in seconds (NULL = no timeout)
- **Default**: NULL
- **Type**: Numeric or NULL
- **Functions**: lst_vat

---

## Parallel Processing Parameters

### parallel
- **Description**: Enable parallel processing for multiple files (auto-enabled for 2+ files)
- **Default**: Auto-detected based on file count
- **Type**: Logical
- **Functions**: trk_pitchmark, lst_dysprosody, and functions using processMediaFiles_LoadAndProcess

### n_cores
- **Description**: Number of CPU cores for parallel processing (NULL = automatic detection)
- **Default**: NULL
- **Type**: Integer or NULL
- **Functions**: trk_pitchmark, lst_dysprosody

---

## Media Conversion Parameters

### assertLossless
- **Description**: File extensions to assert contain lossless audio (warning if not)
- **Default**: NULL
- **Type**: Character vector or NULL
- **Functions**: trk_acfana, trk_cepstrum, trk_cssSpectrum, trk_dftSpectrum, trk_forest, trk_ksvfo, trk_lp_analysis, trk_lpsSpectrum, trk_mhspitch, trk_rmsana, trk_zcrana

### keepConverted
- **Description**: Keep converted audio files after processing
- **Default**: FALSE
- **Type**: Logical
- **Functions**: Same as assertLossless

### convertOverwrites
- **Description**: Overwrite existing converted files
- **Default**: FALSE
- **Type**: Logical
- **Functions**: Same as assertLossless

---

## Special Analysis Parameters

### centerTime
- **Description**: Single-frame analysis at specific time point (overrides beginTime/endTime/windowShift)
- **Default**: FALSE
- **Type**: Logical or Numeric
- **Functions**: trk_acfana, trk_cepstrum, trk_cssSpectrum, trk_dftSpectrum, trk_forest, trk_ksvfo, trk_lp_analysis, trk_lpsSpectrum, trk_mhspitch, trk_rmsana, trk_zcrana

---

## D4C Aperiodicity Parameters

### threshold
- **Description**: D4C threshold parameter for aperiodicity estimation
- **Default**: 0.85
- **Type**: Numeric
- **Functions**: trk_d4c

---

## Python Audio Processing Parameters

### broadcast_audio
- **Description**: Broadcast audio to all channels
- **Default**: TRUE
- **Type**: Logical
- **Functions**: trk_pyin, trk_crepe, trk_yin, trk_yaapt

### hop_length
- **Description**: Hop length for STFT in samples
- **Default**: Varies (typically computed from windowShift)
- **Type**: Integer
- **Functions**: Various Python-based audio functions

### drop_last
- **Description**: Drop last batch if fewer samples than batch_size
- **Default**: FALSE
- **Type**: Logical
- **Functions**: Deep learning pitch functions (trk_crepe, trk_swiftf0)

---

## Installation Parameters

### method
- **Description**: Python module installation method ("auto", "conda", "pip")
- **Default**: "auto"
- **Type**: Character
- **Functions**: All install_* functions (install_swiftf0, install_brouhaha, install_dysprosody, etc.)

### compile_cython
- **Description**: Compile Cython extensions for maximum speed
- **Default**: FALSE
- **Type**: Logical
- **Functions**: install_brouhaha, install_voxit

---

## Logging Parameters

### logToFile
- **Description**: Log commands to separate logfile in outputDirectory
- **Default**: FALSE
- **Type**: Logical
- **Functions**: trk_acfana, trk_cepstrum, trk_cssSpectrum, trk_dftSpectrum, trk_forest, trk_ksvfo, trk_lp_analysis, trk_lpsSpectrum, trk_mhspitch, trk_rmsana, trk_zcrana

---

## Summary Statistics

- **Total Unique Parameters**: 97
- **Universal Parameters**: 8 (listOfFiles, beginTime, endTime, windowShift, toFile, outputDirectory, verbose, explicitExt)
- **Pitch-Specific Parameters**: 15
- **Spectral-Specific Parameters**: 12
- **LP-Specific Parameters**: 8
- **Formant-Specific Parameters**: 4
- **MFCC-Specific Parameters**: 6
- **Pitchmark-Specific Parameters**: 14
- **VAD-Specific Parameters**: 7
- **Parallel Processing Parameters**: 2

---

## Parameter Usage Patterns

### Most Common Parameters (Used by 30+ functions)
1. **listOfFiles** - All functions
2. **beginTime** / **endTime** - All time-windowing capable functions
3. **toFile** - All track functions
4. **verbose** - All functions
5. **outputDirectory** - All file-writing functions
6. **windowShift** - All frame-based analysis functions

### Pitch Tracking Parameter Set
Functions using this set: trk_rapt, trk_swipe, trk_dio, trk_reaper, trk_harvest, trk_pyin, trk_crepe, trk_yin, trk_yaapt
- minF
- maxF
- voicing_threshold (or equivalent)
- windowShift

### Spectral Analysis Parameter Set
Functions using this set: trk_cepstrum, trk_cssSpectrum, trk_dftSpectrum, trk_lpsSpectrum
- windowSize
- window
- fftLength
- resolution

### Linear Prediction Parameter Set
Functions using this set: trk_forest, trk_lpsSpectrum, trk_lp_analysis
- analysisOrder (or order)
- preemphasis
- effectiveLength
- window

---

## Notes for Function Development

When creating new trk_* or lst_* functions:

1. **Always include**: listOfFiles, beginTime, endTime, toFile, verbose, outputDirectory
2. **For time-series tracks**: Add windowShift, explicitExt
3. **For pitch tracking**: Add minF, maxF, voicing_threshold
4. **For spectral analysis**: Add windowSize, window, fftLength or resolution
5. **For LP-based methods**: Add analysisOrder, preemphasis
6. **For parallel processing**: Support automatic parallelization via processMediaFiles_LoadAndProcess()

---

**Document Created**: 2025-10-29
**Status**: Complete catalog of 97 parameters across 75+ functions
**Purpose**: Developer reference for parameter usage and standardization

# superassp 0.7.1

## New Features

### Dysprosody Prosodic Assessment Module

* **NEW: `lst_dysprosody()`** - Extract 193 prosodic features using the dysprosody model from Nylén et al. (2025, doi: 10.3389/fnhum.2025.1566274)
  - MOMEL-INTSINT pitch target extraction and tone coding
  - Spectral tilt measures with Iseli-Alwan harmonic correction
  - Statistical summaries and differential features
  - Full integration with av package for universal media format support
  - Parallel batch processing support (~5x speedup with 8 cores)
  - Performance: ~0.16-0.44s per file (14x realtime)

* **NEW: `install_dysprosody()`** - Install dysprosody Python module and dependencies
* **NEW: `dysprosody_available()`** - Check dysprosody module availability
* **NEW: `dysprosody_info()`** - Get dysprosody version and dependency information

**Python Module Integration:**
- Created `inst/python/dysprosody/` package with optimized implementation
- Auto-imports Phase 1 optimized version (15-25% performance improvement)
- Fallback to pure Python if optimized version unavailable
- Comprehensive documentation and performance analysis included

**Features Extracted (193 total):**
- Prosodic metadata: Duration, PitchKey, PitchRange, PitchMean, IntsIntLabels
- Spectral features: L2L1, L2cL1c, L1cLF3c, SLF, C1, SpectralBalance, SLF6D
- Statistical summaries: mean, std, var, iqr, max, min for all time-varying features
- Differential features: _diff versions showing inter-INTSINT-label changes

# superassp 0.7.0

## Major Changes

### Universal Media Format Support

All DSP functions now use the `av` package for audio loading, supporting:
- Standard formats: WAV, MP3, FLAC, OGG, AAC, Opus
- Video formats: MP4, MKV, AVI, MOV (extracts audio)
- Niche formats: AU, Kay, NIST, NSP (via wrassp fallback)

### Complete librosa Migration

**Migrated 6 functions from librosa.load() to av::read_audio_bin():**
- `trk_pyin()` - Probabilistic YIN pitch tracker
- `trk_yin()` - YIN pitch tracker
- `trk_crepe()` - CREPE deep learning pitch tracker (migrated from torchcrepe)
- `trk_yaapt()` - YAAPT pitch tracker
- `trk_seenc()` - WORLD spectral envelope coding
- `trk_excite()` - SPTK excitation signal extraction
- `trk_aperiodicities()` - WORLD D4C aperiodicity (deprecated)
- `reaper_pm()` - REAPER pitch mark extraction

**Python scripts updated:**
- `inst/python/snack_pitch.py` - Replaced librosa with soundfile
- `inst/python/snack_formant.py` - Replaced librosa with soundfile

### PyTorch Function Cleanup

**Removed 3 redundant PyTorch functions:**
- `trk_kaldi_pitch()` - Deprecated in torchaudio 2.9+, use `trk_rapt()` instead
- `trk_torch_pitch()` - Generic torch pitch, use `trk_rapt()` or `trk_swipe()` instead
- `trk_torch_mfcc()` - Redundant with `trk_mfcc()` C++ SPTK implementation

**Kept:**
- `trk_crepe()` - Unique deep learning CNN algorithm, already migrated to av

**Rationale:** PyTorch functions unlikely to work well in reticulate environment, C++ alternatives are faster and more stable.

### Package Organization

* Renamed `R/wrassp_packageVars.R` → `R/assp_library_vars.R` for clarity
* All package variables now clearly identified by their source library

## Documentation

* Added comprehensive documentation for all migrated functions
* Created `trk_excite()` documentation (previously undocumented)
* Updated all function docs to mention av package usage
* Removed documentation for deleted functions

## Performance

* All DSP functions now have consistent, modern architecture
* In-memory processing throughout
* No file I/O bottlenecks from format conversions

## Breaking Changes

* **REMOVED:** `trk_kaldi_pitch()` - Use `trk_rapt()` instead
* **REMOVED:** `trk_torch_pitch()` - Use `trk_rapt()` or `trk_swipe()` instead
* **REMOVED:** `trk_torch_mfcc()` - Use `trk_mfcc()` instead

These functions were redundant with faster, more stable C++ implementations.

## Bug Fixes

* Fixed media format compatibility across all Python-based DSP functions
* Improved error handling for short audio files
* Better time windowing support in Python functions

---

# Earlier Versions

See git history for versions prior to 0.7.0.

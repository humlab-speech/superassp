# Python Environment Documentation

This document describes the Python environment used by the superassp R package for Python-based speech analysis functions.

## Python Configuration (via reticulate)

```
Python:      /opt/miniconda3/bin/python3
Version:     3.12.9 (packaged by Anaconda, Inc.)
libpython:   /opt/miniconda3/lib/libpython3.12.dylib
pythonhome:  /opt/miniconda3
```

**Note:** Python version is forced by `RETICULATE_PYTHON` environment variable.

## Core Scientific Libraries

| Package | Version | Purpose |
|---------|---------|---------|
| numpy | 2.2.6 | Numerical computing, array operations |
| scipy | 1.16.1 | Scientific computing, signal processing |
| pandas | 2.3.2 | Data manipulation, CSV/dataframe handling |
| soundfile | 0.13.1 | Audio I/O (WAV, FLAC, etc.) |
| librosa | 0.11.0 | Audio analysis and feature extraction |

## Speech Analysis Libraries

| Package | Version | Purpose | Used By |
|---------|---------|---------|---------|
| praat-parselmouth | 0.4.6 | Praat functionality in Python | `trk_praat_sauce()`, `trk_pitchp()`, `trk_formantp()` |
| pysptk | 1.0.1 | Speech Signal Processing Toolkit | Spectral analysis functions |
| pyworld | 0.3.5 | WORLD vocoder | `trk_aperiodicities()` |
| praatdet-py | 0.1.0 | Praat formant detection | Development/testing |

## Deep Learning Libraries

| Package | Version | Purpose | Used By |
|---------|---------|---------|---------|
| tensorflow | 2.20.0 | Deep learning framework | `trk_deepformants()`, neural network models |
| torch (PyTorch) | 2.5.1 | Deep learning framework | `trk_deepformants()`, `trk_crepe()` |
| torchaudio | 2.5.1 | Audio processing for PyTorch | Audio loading, resampling |
| pytorch-lightning | 2.5.5 | PyTorch training framework | Model training utilities |

## Optional/Development Packages

| Package | Version | Notes |
|---------|---------|-------|
| torch-audiomentations | 0.12.0 | Data augmentation |
| pytorch-metric-learning | 2.9.0 | Metric learning |
| torchmetrics | 1.8.2 | Evaluation metrics |

## Functions Requiring Python

### Pitch Tracking (Python-based)
- `trk_pitchp()` - Praat pitch tracking (parselmouth)
- `trk_crepe()` - CREPE neural network pitch tracking (PyTorch) *[if installed]*
- `trk_pyin()` - Probabilistic YIN (librosa) *[if installed]*
- `trk_straight_f0()` - STRAIGHT F0 extraction (scipy) **[Currently has segfault issue]**

### Formant Analysis (Python-based)
- `trk_praat_sauce()` - Praat formant tracking with voice quality (parselmouth) ✓ Working
- `trk_formantp()` - Praat formant tracking (parselmouth)
- `trk_deepformants()` - Deep learning formant tracking (PyTorch) *[Requires model download]*

### Voice Quality Analysis (Python-based)
- `trk_aperiodicities()` - WORLD vocoder aperiodicity estimation (pyworld)
- `lst_voice_sauce()` - VoiceSauce measurements (parselmouth + scipy)

## Known Issues

### 1. trk_straight_f0 Segfault
**Status:** Known limitation
**Issue:** Segfaults when calling scipy C extensions through reticulate
**Workaround:** Use alternative pitch trackers (`trk_rapt()`, `trk_swipe()`, `trk_reaper()`)

### 2. trk_deepformants Model Requirements
**Status:** Working when models installed
**Issue:** Requires PyTorch models to be downloaded separately
**Solution:** Run `install_deepformants()` to download models

## Installation Instructions

### Conda Environment Setup

```bash
conda create -n superassp python=3.12
conda activate superassp

# Core scientific libraries
conda install numpy scipy pandas soundfile librosa

# Speech analysis
pip install praat-parselmouth pysptk pyworld

# Deep learning (optional)
conda install pytorch torchaudio -c pytorch
pip install tensorflow
```

### R Configuration

Set environment variable in `.Renviron`:
```
RETICULATE_PYTHON=/opt/miniconda3/bin/python3
```

Or set in R before loading superassp:
```r
Sys.setenv(RETICULATE_PYTHON = "/opt/miniconda3/bin/python3")
library(superassp)
```

## Verification

Check Python configuration:
```r
library(superassp)
reticulate::py_config()
```

Test Python-based functions:
```r
test_file <- system.file("samples/sustained/a32b.wav", package = "superassp")

# Test Praat formant tracking
result <- trk_praat_sauce(test_file, toFile = FALSE)

# Test Praat pitch tracking
result <- trk_pitchp(test_file, toFile = FALSE)

# Test WORLD aperiodicities
result <- trk_aperiodicities(test_file, toFile = FALSE)
```

## Performance Benchmarks

Based on tests with 4.04s audio file (100 iterations):

| Function | Median Time | Python Library |
|----------|-------------|----------------|
| `trk_praat_sauce` | 1107 ms | parselmouth |
| `trk_formantp` | 893 ms | parselmouth |
| `trk_pitchp` | ~200-300 ms | parselmouth |
| `trk_aperiodicities` | ~500-700 ms | pyworld |

**Note:** Python-based functions are generally slower than C++ implementations but provide additional functionality and compatibility with Praat.

## Troubleshooting

### Issue: "Python module not found"
**Solution:** Install the required package in the conda environment:
```bash
pip install <package-name>
```

### Issue: "reticulate cannot find Python"
**Solution:** Set `RETICULATE_PYTHON` environment variable:
```r
Sys.setenv(RETICULATE_PYTHON = "/opt/miniconda3/bin/python3")
```

### Issue: "ImportError: No module named..."
**Solution:** Verify package installation:
```bash
/opt/miniconda3/bin/python3 -m pip list | grep <package>
```

## Last Updated

2025-11-02 - Verified during benchmark testing and bug fix session

## See Also

- [BUGFIX_SESSION_2025-11-02.md](BUGFIX_SESSION_2025-11-02.md) - Recent bug fixes
- [requirements.txt](requirements.txt) - Python package requirements
- [setup.py](setup.py) - Package setup configuration

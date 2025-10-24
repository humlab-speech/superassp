# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

This repository implements MOMEL (MOdelling MELody) and INTSINT (INternational Transcription System for INTonation) algorithms for prosodic analysis of speech. It contains:

- **Python implementations** using Parselmouth (Praat bindings) for modern prosodic feature extraction
- **Legacy binaries** for the original C implementation of MOMEL
- **Perl script** for INTSINT coding

The primary goal is computing prosodic measures from audio files as described in doi: 10.3389/fnhum.2025.1566274.

## Key Files

### Primary Implementation Files

- **dysprosody_pure.py** - ⭐ **RECOMMENDED** Pure Python implementation of `prosody_measures()` with no external dependencies on compiled binaries or Perl. This is the optimized version that should be used for new projects. Features:
  - Complete Python implementation using momel_intsint.py
  - Faithful to the published dysprosody implementation
  - Platform-independent (no compiled binaries needed)
  - Returns identical feature sets as original implementation
  - ~0.2-0.6s processing time per file

- **dysprosody.py** - Original implementation for prosodic feature extraction. Contains the `prosody_measures()` function which:
  - Extracts F0 using automatic min/max estimation
  - Runs MOMEL and INTSINT via subprocess calls to C binaries and Perl
  - Computes spectral tilt measures (C1, SLF, spectral balance, etc.)
  - Returns comprehensive prosodic statistics as pandas Series
  - ⚠️ Requires working momel binaries and Perl installation

- **momel_intsint.py** - Pure Python reimplementation of MOMEL-INTSINT algorithms using Parselmouth. Contains:
  - `momel()` - Pitch target extraction using quadratic spline modeling
  - `intsint()` - Optimal tone label assignment
  - `process_momel_intsint()` - Complete pipeline

- **demo_pure.py** - Demonstration script showing how to use dysprosody_pure.py for single files or batch processing

### Legacy Components

- **momel.c / momel.h** - Original C implementation of MOMEL algorithm
- **momel_linux**, **momel_osx_intel**, **momel_osx_ppc**, **momel_win.exe** - Compiled binaries for different platforms
- **intsint.pl** - Perl script for INTSINT coding (requires Perl installation)

## Architecture

### Two Implementation Approaches

1. **Subprocess-based** (dysprosody.py):
   - Calls compiled MOMEL binaries via subprocess
   - Pipes pitch values as stdin, receives targets as stdout
   - Calls Perl script for INTSINT coding
   - Platform detection: macOS uses `momel_osx_intel`, Linux uses `momel_linux`

2. **Pure Python** (momel_intsint.py):
   - Complete reimplementation without external dependencies
   - Uses Parselmouth for all Praat operations
   - More maintainable, no subprocess overhead
   - Implements same algorithms: glitch elimination, quadratic regression, target clustering

### Data Flow in dysprosody.py

```
Audio file → automatic_min_max_fo() → Pitch object
                                      ↓
                      momel() → subprocess call → MOMEL targets (PitchTier)
                                      ↓
           code_with_intsint() → Perl subprocess → TextGrid with 3 tiers
                                      ↓
                      prosody_measures() → Statistical features (pandas Series)
```

### Key Algorithms

**Automatic F0 Range Estimation** (dysprosody.py:93-100):
- Two-pass pitch extraction
- First pass: use default range
- Compute min_fo = floor(Q25 * 0.75/10) * 10
- Compute max_fo = ceil(min_fo * 2^pitch_span / 10) * 10
- Second pass with optimized range

**MOMEL Target Extraction** (momel_intsint.py:469-509):
1. Eliminate glitches (pitch jumps > RAPP_GLITCH)
2. Find candidates via quadratic regression (`cible()`)
3. Cluster and reduce targets (`reduc()`)
4. Add boundary targets (`borne()`)

**INTSINT Coding** (momel_intsint.py:604-659):
- Optimizes over range (MIN_RANGE to MAX_RANGE) and key (mean ± MEAN_SHIFT)
- Assigns tone labels: M (mid), T (top), B (bottom), H (higher), L (lower), U (up), D (down), S (same)
- Minimizes sum-squared error between observed and estimated targets

**Spectral Tilt Measures** (dysprosody.py:297-366):
- Extracts harmonics L1-L4 with formant correction (Iseli-Alwan algorithm)
- Computes bandwidth using Hawks-Miller estimation
- Returns: L2L1, L2cL1c, L1cLF3c, L1LF3, SLF, C1, Spectral Balance, SLF6D coefficients

## Usage

### ⭐ Recommended: Pure Python implementation (dysprosody_pure.py)

```python
from dysprosody_pure import prosody_measures
import pandas as pd
import glob
import os

# Single file analysis
features = prosody_measures("audio.wav")
# Returns pandas.Series with 193 prosodic features

# Access specific features
print(f"Duration: {features['Duration']:.2f}s")
print(f"Pitch Mean: {features['PitchMean']:.0f} Hz")
print(f"INTSINT Labels: {features['IntsIntLabels']:.0f}")

# Batch processing - recommended approach
files = glob.glob("**/*.wav", recursive=True)
results = {}

for wav_file in files:
    result = prosody_measures(wav_file)
    if result is not None:  # Skip files < 1 second
        basename = os.path.basename(wav_file).removesuffix('.wav')
        results[basename] = result

# Create DataFrame with files as rows, features as columns
df = pd.DataFrame(results).T
df.to_csv("prosody_results.csv")
```

### Demo script

```bash
# Analyze a single file
python demo_pure.py audio.wav

# Batch analyze all WAV files in current directory
python demo_pure.py
```

### Original implementation (dysprosody.py) - requires binaries/Perl

```python
from dysprosody import prosody_measures

# Single file (requires momel binary and intsint.pl)
features = prosody_measures("/path/to/audio.wav")
```

### Low-level MOMEL-INTSINT only (momel_intsint.py)

```python
from momel_intsint import process_momel_intsint

intsint_targets, range_oct, key = process_momel_intsint("audio.wav")

for target in intsint_targets:
    print(f"{target.time:.3f} {target.tone} {target.target:.0f} {target.estimate:.0f}")
```

## Dependencies

### For dysprosody_pure.py (recommended)

- **parselmouth** - Python bindings for Praat
- **numpy** - Numerical operations
- **pandas** - Data structures and analysis
- **scipy** - Statistical functions

### Additional requirements for dysprosody.py (original)

- **Perl** - Required for intsint.pl
- **momel binary** - Platform-specific compiled binary (momel_osx_intel, momel_linux, etc.)

## Platform Compatibility

### dysprosody_pure.py (recommended)
- **All platforms**: Fully platform-independent, works on macOS, Linux, Windows, ARM, x86_64

### dysprosody.py (original)
- **macOS**: Uses `momel_osx_intel` (x86_64 binary)
- **Linux**: Uses `momel_linux` (i386 binary, may require 32-bit libs)
- **Windows**: `momel_win.exe` available but not in subprocess detection logic
- **ARM/M1/M2 Macs**: May have compatibility issues with x86_64 binary

## Important Notes

### General
- Files < 1 second duration are skipped by `prosody_measures()`
- Returns pandas.Series with 193 features including:
  - Prosodic metadata: Duration, PitchKey, PitchRange, PitchMean, IntsIntLabels, UniqueIntsInt
  - Spectral features: L2L1, L2cL1c, L1cLF3c, SLF, C1, Spectral Balance, SLF6D coefficients
  - Statistical summaries: mean, std, variation, iqr, max, min for all time-varying features
  - Differential features: "_diff" versions showing inter-INTSINT-label changes

### dysprosody_pure.py (recommended)
- **No external dependencies** on compiled binaries or Perl
- **Faithful reproduction** of published implementation output
- **Parameter matching**: Uses identical MOMEL/INTSINT parameters as original (window_length=30ms, reduced_window_length=20ms, minimal_distance=5 frames)
- **Processing time**: ~0.2-0.6s per file depending on duration

### dysprosody.py (original)
- All MOMEL binaries and intsint.pl must be in current working directory
- MOMEL outputs time in milliseconds, INTSINT expects milliseconds but outputs seconds
- Uses hardcoded paths (`./momel_osx_intel`, `./momel_linux_intel`, `intsint.pl`)
- Platform detection at dysprosody.py:107-110 checks `platform.system() == "Darwin"` for macOS
- May fail with "division by zero" error in Perl script for certain audio files

## Citations

When using this code, cite:
- Original paper: doi: 10.3389/fnhum.2025.1566274
- Hirst, D., & Espesser, R. (1993). Automatic Modelling Of Fundamental Frequency Using A Quadratic Spline Function
- Hirst, D. (2019). INTSINT: a new algorithm using the OMe scale
- See dysprosody.py:38-47 for complete citation list

## License

Distributed under CC BY 4.0 license (https://creativecommons.org/licenses/by/4.0/)

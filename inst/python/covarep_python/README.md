# COVAREP Python Implementation

Python reimplementation of the COVAREP (Cooperative Voice Analysis Repository for Speech Technologies) toolkit.

## Project Status

**Phase:** 1 - Foundation (Pilot)  
**Started:** October 17, 2025  
**Current Focus:** Core infrastructure and voicebox compatibility

## Features

### Implemented ✅
- Project structure
- Basic infrastructure

### In Progress 🚧
- Voicebox compatibility layer
- F0 tracking (pitch_srh)
- IAIF (Iterative Adaptive Inverse Filtering)

### Planned ⏳
- Glottal source analysis
- Voice quality parameters
- Envelope estimation
- Feature extraction pipeline
- Vocoder and sinusoidal modeling

## Installation

```bash
cd covarep_python
pip install -r requirements.txt
```

### Development Installation

```bash
pip install -e .
```

## Quick Start

```python
from covarep import F0Tracker
import soundfile as sf

# Load audio
audio, fs = sf.read("speech.wav")

# Extract F0
f0_tracker = F0Tracker(method='srh')
f0, vuv = f0_tracker.estimate(audio, fs)
```

## Project Structure

```
covarep_python/
├── covarep/              # Main package
│   ├── glottal/         # Glottal source analysis
│   ├── envelope/        # Spectral envelope estimation
│   ├── features/        # Feature extraction
│   ├── vocoder/         # Vocoder algorithms
│   ├── sinusoidal/      # Sinusoidal modeling
│   ├── f0/              # F0 estimation
│   └── voicebox/        # Voicebox compatibility layer
├── tests/               # Unit tests
├── examples/            # Usage examples
└── docs/                # Documentation
```

## Original COVAREP

This is a Python reimplementation of:
- **Original:** https://github.com/covarep/covarep
- **Citation:** G. Degottex, J. Kane, T. Drugman, T. Raitio and S. Scherer, "COVAREP - A collaborative voice analysis repository for speech technologies", In Proc. IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), Florence, Italy 2014.

## License

LGPL (following original COVAREP licensing)

## Authors

Python Implementation: 2025
Original MATLAB: Gilles Degottex, John Kane, Thomas Drugman, Tuomo Raitio, Stefan Scherer, et al.

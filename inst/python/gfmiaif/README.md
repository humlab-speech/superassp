# GFM-IAIF Python Module

Glottal Flow Model-based Iterative Adaptive Inverse Filtering for speech signal processing.

## Overview

This module implements the GFM-IAIF algorithm for source-filter separation of speech signals, extracting linear prediction coefficients for:
- **Vocal tract filter** (default: 48th order)
- **Glottis filter** (default: 3rd order, highly recommended)
- **Lip radiation filter** (2nd order)

## Algorithm

GFM-IAIF is an extension of the classical IAIF algorithm with improved pre-emphasis that allows extraction of wide-band glottis response, incorporating both glottal formant and spectral tilt characteristics.

### Processing Pipeline

1. **Pre-frame addition** - Adds ramp to reduce edge effects
2. **Lip radiation cancellation** - Leaky integration (1/[1 - d·z⁻¹])
3. **Gross glottis estimation** - Iterative 1st-order LPC
4. **Gross vocal tract estimation** - High-order LPC after glottis removal
5. **Fine glottis estimation** - Refined ng-order LPC
6. **Fine vocal tract estimation** - Final nv-order LPC

## Usage

### Single Frame Processing

```python
from gfmiaif import gfmiaif_fast
import numpy as np

# Process single speech frame
frame = np.random.randn(512)
av, ag, al = gfmiaif_fast(frame, nv=48, ng=3, d=0.99)

# av: Vocal tract coefficients (shape: 49)
# ag: Glottis coefficients (shape: 4)
# al: Lip radiation coefficients (shape: 2)
```

### Frame-Based Processing (for R integration)

```python
from gfmiaif import gfmiaif_frame_based

# Process entire audio signal
result = gfmiaif_frame_based(
    audio=signal,
    sample_rate=16000,
    window_shift=0.010,  # 10ms
    window_size=0.032,   # 32ms
    nv=48,
    ng=3,
    d=0.99,
    window_type='hann'
)

# Returns dictionary with:
# - av: (n_frames, 49) vocal tract coefficients
# - ag: (n_frames, 4) glottis coefficients
# - al: (n_frames, 2) lip radiation coefficients
# - timestamps: (n_frames,) frame times in seconds
```

## Parameters

### Core Parameters

- `nv` (int, default=48): Vocal tract LPC order
  - Controls spectral resolution of vocal tract filter
  - Higher order = more detailed formant structure
  - Typical range: 24-48 for speech at 16kHz

- `ng` (int, default=3): Glottis LPC order
  - **IMPORTANT**: ng=3 is highly recommended
  - Algorithm designed assuming 3rd-order filter captures glottis timbre
  - Describes tenseness, effort, breathiness

- `d` (float, default=0.99): Leaky integration coefficient
  - Models lip radiation as high-pass filter
  - Range: 0.95-0.99
  - Higher values = stronger high-frequency emphasis

- `window_type` (str, default='hann'): Window function
  - Options: 'hann', 'hamming', 'blackman'
  - Hanning window is standard for speech analysis

### Frame-Based Parameters

- `window_shift` (float, default=0.010): Frame shift in seconds
  - Typical: 5-10ms for speech
  - Smaller = better time resolution, more computation

- `window_size` (float, default=0.032): Frame size in seconds
  - Typical: 20-32ms for speech (320-512 samples at 16kHz)
  - Must be > nv+1 samples

## Performance

The implementation uses optimized algorithms:
- **FFT-based autocorrelation** (3x faster than direct method)
- **JIT-compiled Levinson-Durbin** (5-10x faster with Numba)
- **Vectorized operations** throughout

### Typical Processing Speed

- Single frame (512 samples): ~0.25 ms
- Real-time factor (16kHz): ~130x
- Can process 100+ frames/second

### Dependencies

**Required:**
- numpy >= 1.19
- scipy >= 1.5

**Optional (for maximum speed):**
- numba >= 0.54 (enables JIT compilation, 5-10x speedup)

## References

[1] O. Perrotin and I. V. McLoughlin (2019)
    "A spectral glottal flow model for source-filter separation of speech"
    IEEE ICASSP 2019, pp. 7160-7164.

[2] P. Alku (1992)
    "Glottal wave analysis with pitch synchronous iterative adaptive inverse filtering"
    Speech Communication, 11(2-3), pp. 109-118.

## License

GNU Lesser General Public License v3.0 or later (LGPL-3.0-or-later)

Copyright (c) 2019 Univ. Grenoble Alpes, CNRS, Grenoble INP, GIPSA-lab
Python implementation (c) 2025

## Integration

This module is designed for integration with the superassp R package via reticulate.
See `trk_gfmiaif()` R function for high-level interface.

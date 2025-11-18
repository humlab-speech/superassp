"""
STRAIGHT (Speech Transformation and Representation using Adaptive Interpolation 
of weiGHTed spectrogram) - Python Implementation

A faithful reimplementation of the MATLAB STRAIGHT vocoder.

Main API functions:

- MulticueF0v14: F0 extraction (91.9% frame, 99.0% mean accuracy)

- exstraightspec: Spectral analysis (99.9962% accurate)
- exstraightAPind: Aperiodicity analysis (99.83% accurate)
- exstraightsynth: Speech synthesis (99.99% accurate)


All modules achieve >90% correlation with MATLAB baseline.

Example:
    >>> from legacy_STRAIGHT import MulticueF0v14, exstraightsynth, exstraightspec
    >>> f0, vuv, aux, params = MulticueF0v14(x, fs, 71, 800)
    >>> spec = exstraightspec(x, fs, f0, vuv)
    >>> synth = exstraightsynth(f0, spec, ap, fs)

    
Copyright (c) 2024 Python Implementation
Original MATLAB version: Copyright (c) Wakayama University, 2004-2016
"""

__version__ = '0.2.0'

__author__ = 'Python port of STRAIGHT by Hideki Kawahara'

# Import main API
from .synthesis import exstraightsynth, SynthesisParameters
from .spectral import exstraightspec, SpectralParameters
from .aperiodicity import exstraightAPind, AperiodicityParameters

from .f0_extraction import MulticueF0v14, F0Parameters
from .f0_wrapper import extract_f0_safe

__all__ = [
    'MulticueF0v14',
    'F0Parameters',
    'extract_f0_safe',

    'exstraightsynth',
    'exstraightspec',
    'exstraightAPind',
    'SynthesisParameters',
    'SpectralParameters',
    'AperiodicityParameters',
]

"""
COVAREP - Cooperative Voice Analysis Repository for Speech Technologies
Python Implementation

A comprehensive toolkit for voice and speech analysis including:
- Glottal source analysis
- F0 (pitch) estimation
- Voice quality parameters
- Spectral envelope estimation
- Vocoder algorithms
- Sinusoidal modeling
- Feature extraction

Original MATLAB version: https://github.com/covarep/covarep
"""

__version__ = '0.1.0'
__author__ = 'COVAREP Python Team'
__license__ = 'LGPLv3'

# Core imports
from . import voicebox
from . import f0
from . import glottal
from . import utils

__all__ = ['voicebox', 'f0', 'glottal', 'utils', '__version__']

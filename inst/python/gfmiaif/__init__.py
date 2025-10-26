"""
GFM-IAIF: Glottal Flow Model-based Iterative Adaptive Inverse Filtering

This module provides optimized implementations of the GFM-IAIF algorithm
for speech signal processing, extracting linear prediction coefficients
for vocal tract, glottis, and lip radiation filters.

References
----------
[1] O. Perrotin and I. V. McLoughlin (2019)
    "A spectral glottal flow model for source-filter separation of speech",
    IEEE ICASSP 2019, pp. 7160-7164.

Copyright (c) 2019 Univ. Grenoble Alpes, CNRS, Grenoble INP, GIPSA-lab
Python implementation (c) 2025

License: GNU Lesser General Public License v3.0 or later
"""

from .gfmiaif import gfmiaif_fast, gfmiaif_frame_based
from .lpc import lpc_fast

__all__ = ['gfmiaif_fast', 'gfmiaif_frame_based', 'lpc_fast']
__version__ = '1.0.0'

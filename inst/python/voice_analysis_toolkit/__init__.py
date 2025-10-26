"""
Voice Analysis Toolkit - Python Implementation

This is a Python conversion of the MATLAB Voice Analysis Toolkit developed by
John Kane at the Phonetics and Speech Laboratory in Trinity College Dublin.

The toolkit contains a range of functions for glottal source and voice quality
analysis, including:
- SE_VQ: Glottal closure instant (GCI) detection
- IAIF: Iterative adaptive inverse filtering
- Creak detection: Detection and classification of creaky voice
- MDQ: Maxima dispersion quotient for voice quality analysis
- peakSlope: Peak slope measure extraction

Original MATLAB code: (C) 2009-2013 John Kane, Trinity College Dublin
Python conversion: 2025

License: GNU GPL v2 for non-commercial use
"""

__version__ = "1.0.0"

from . import general
from . import se_vq
from . import creak
from . import utils

__all__ = ['general', 'se_vq', 'creak', 'utils']

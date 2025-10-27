"""
Union Method Creak Detection - Optimized for R Integration
Python package for automatic creaky voice detection using the Union Method.

Optimized for use with superassp/protoscribe R packages:
- Accepts NumPy arrays from av package (in-memory processing)
- Returns NumPy arrays for efficient transfer to R
- No disk I/O within Python (handled by R)
- Vectorized operations for speed
"""

__version__ = "1.0.0"
__all__ = ['CreakDetector', 'CreakDetectorExtended', 'extract_creak_features', 
           'detect_creak_union', 'detect_creak_union_extended']

from .detector import CreakDetector, CreakDetectorExtended
from .features import get_all_features as extract_creak_features
from .api import detect_creak_union, detect_creak_union_extended

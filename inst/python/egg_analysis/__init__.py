"""
EGG Analysis for superassp

Python DSP module for electroglottographic (EGG) signal analysis.
Adapted from egg_python package for integration with superassp.

This module provides f0 and open quotient (Oq) extraction from EGG signals,
returning results at equal intervals suitable for SSFF track format.
"""

from .egg_f0 import analyze_egg_f0

__all__ = ['analyze_egg_f0']
__version__ = '1.0.0'

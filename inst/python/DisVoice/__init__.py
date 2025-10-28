"""
DisVoice - Optimized speech feature extraction for superassp

This module provides in-memory, optimized implementations of speech analysis
functions using Parselmouth (Praat Python interface).

Key optimizations:
- 1.5-2.5x faster than file-based approaches
- No temporary file I/O
- Direct numpy array processing
- Compatible with av package audio loading

Main modules:
- praat_functions: Pitch, formant, and voicing extraction
"""

from . import praat_functions

__version__ = "0.1.8"
__all__ = ["praat_functions"]

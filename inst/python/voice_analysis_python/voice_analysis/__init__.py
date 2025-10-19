"""
Voice Analysis Toolbox - Python Implementation

A faithful reimplementation of the MATLAB Voice Analysis Toolbox by Athanasios Tsanas
for computing 132 dysphonia measures from sustained vowel recordings.

Original MATLAB code: Copyright (c) Athanasios Tsanas, 2014
Python implementation: 2025

License: GPL-3.0

Citation Requirements:
    Tsanas, A., Little, M., McSharry, P., & Ramig, L. (2011).
    Nonlinear speech analysis algorithms mapped to a standard metric achieve
    clinically useful quantification of average Parkinson's disease symptom severity.
    Journal of the Royal Society Interface, 8(59), 842-855.
"""

__version__ = "0.1.0"
__author__ = "Voice Analysis Team"

from .core import VoiceAnalyzer, analyze_voice, analyze_voice_file
from .core_parallel import (
    VoiceAnalyzerParallel, 
    analyze_voice_parallel, 
    analyze_voice_file_parallel,
    analyze_batch_parallel
)

__all__ = [
    "VoiceAnalyzer", 
    "analyze_voice", 
    "analyze_voice_file",
    "VoiceAnalyzerParallel",
    "analyze_voice_parallel",
    "analyze_voice_file_parallel",
    "analyze_batch_parallel"
]

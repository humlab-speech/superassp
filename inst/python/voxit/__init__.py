"""
Voxit: Voice and articulation complexity measures

This module provides acoustic summaries for speech analysis, computing
prosodic and rhythmic complexity measures from audio and word alignments.

Features:
- Speaking rate (WPM)
- Pause analysis (counts, durations, rates)
- Rhythmic complexity of pauses
- Pitch statistics (range, entropy, speed, acceleration)

References:
- Original MATLAB implementation: https://github.com/[voxit-url]
- Python reimplementation with optimizations
"""

__version__ = "1.0.0"
__author__ = "Voxit Contributors"

from .voxit_core import compute_voxit_features
from .voxit_optimized import compute_voxit_features_optimized

# Try to use optimized version if numba is available
try:
    from .voxit_numba import compute_voxit_features_numba
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

# Try to use Cython version if compiled
try:
    from .voxit_cython import compute_voxit_features_cython
    HAS_CYTHON = True
except ImportError:
    HAS_CYTHON = False

def compute_features(gentle_data, pitch_data, start_time=None, end_time=None,
                    use_numba=False, use_cython=False):
    """
    Compute Voxit features with automatic optimization selection.
    
    Parameters:
    -----------
    gentle_data : list of dict
        Word alignment data with 'word', 'case', 'start', 'end' keys
    pitch_data : list of dict  
        Pitch track data with 'time', 'frequency' keys
    start_time : float, optional
        Start time in seconds for analysis window
    end_time : float, optional
        End time in seconds for analysis window
    use_numba : bool
        Use numba-optimized version if available
    use_cython : bool
        Use cython-optimized version if available
        
    Returns:
    --------
    dict : Dictionary of computed features
    """
    if use_cython and HAS_CYTHON:
        return compute_voxit_features_cython(gentle_data, pitch_data, 
                                             start_time, end_time)
    elif use_numba and HAS_NUMBA:
        return compute_voxit_features_numba(gentle_data, pitch_data,
                                            start_time, end_time)
    else:
        return compute_voxit_features_optimized(gentle_data, pitch_data,
                                                start_time, end_time)

__all__ = [
    'compute_features',
    'compute_voxit_features',
    'compute_voxit_features_optimized',
    'HAS_NUMBA',
    'HAS_CYTHON'
]

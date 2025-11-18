"""
Voxit: Voice and articulation complexity measures

This module provides acoustic summaries for speech analysis, computing
prosodic and rhythmic complexity measures from audio and word alignments.

Features:
- Speaking rate (WPM)
- Pause analysis (counts, durations, rates)
- Rhythmic complexity of pauses
- Pitch statistics (range, entropy, speed, acceleration)


Performance Note:
-----------------
After comprehensive benchmarking, the pure Python implementation provides
the best performance for typical use cases (20-200 words, 2-25ms processing).
Numba JIT compilation introduces overhead that makes it 2-5x slower.
See PERFORMANCE_ANALYSIS.md for details.

References:
- Original MATLAB implementation: https://github.com/NSKarolyn/voxit
- Python reimplementation with validation

"""

__version__ = "1.0.0"
__author__ = "Voxit Contributors"

# Primary implementation - best performance
from .voxit_core import compute_voxit_features

# Optimized implementation (uses NumPy efficiently)
from .voxit_optimized import compute_voxit_features_optimized

# Numba version (available but NOT recommended - see PERFORMANCE_ANALYSIS.md)

try:
    from .voxit_numba import compute_voxit_features_numba
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False


def compute_features(gentle_data, pitch_data, start_time=None, end_time=None,
                    use_numba=False):
    """
    Compute Voxit features.
    
    Note: use_numba=True is NOT recommended. Benchmarking shows the
    pure Python version is 2-5x faster for typical use cases.

    
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

        Use numba-optimized version (NOT RECOMMENDED - slower than Python)

        
    Returns:
    --------
    dict : Dictionary of computed features
    """

    if use_numba and HAS_NUMBA:
        import warnings
        warnings.warn(
            "Numba version is 2-5x slower than pure Python. "
            "See PERFORMANCE_ANALYSIS.md for details.",
            PerformanceWarning
        )
        return compute_voxit_features_numba(gentle_data, pitch_data,
                                            start_time, end_time)
    else:
        # Use optimized version (best performance)

        return compute_voxit_features_optimized(gentle_data, pitch_data,
                                                start_time, end_time)

__all__ = [
    'compute_features',
    'compute_voxit_features',
    'compute_voxit_features_optimized',
    'HAS_NUMBA'

]

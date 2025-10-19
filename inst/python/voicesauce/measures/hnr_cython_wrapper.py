"""
Wrapper for Cython-optimized HNR with fallback to standard implementation
"""

import numpy as np
from typing import List, Dict

# Try to import Cython-compiled module
try:
    from .hnr_cython import get_hnr_cython
    CYTHON_AVAILABLE = True
except ImportError:
    CYTHON_AVAILABLE = False
    get_hnr_cython = None


def get_hnr_optimized(audio: np.ndarray, fs: int, f0_values: np.ndarray,
                     frame_shift: float = 1.0, n_periods: int = 5,
                     freq_bands: List[float] = None,
                     force_cython: bool = False) -> Dict[str, np.ndarray]:
    """
    Calculate HNR using best available implementation

    Tries implementations in order:
    1. Cython (if compiled and available) - 2-3x faster
    2. Scipy.fft optimized (fallback) - 1.2x faster
    3. Standard numpy (final fallback)

    Args:
        audio: Full audio signal
        fs: Sampling rate
        f0_values: F0 array
        frame_shift: Frame shift in milliseconds
        n_periods: Number of pitch periods
        freq_bands: Frequency bands for HNR (default [500, 1500, 2500, 3500])
        force_cython: If True, raises error if Cython not available

    Returns:
        Dictionary with HNR05, HNR15, HNR25, HNR35 arrays
    """
    if freq_bands is None:
        freq_bands = [500.0, 1500.0, 2500.0, 3500.0]

    # Try Cython implementation first
    if CYTHON_AVAILABLE:
        return get_hnr_cython(audio, fs, f0_values, frame_shift, n_periods, freq_bands)

    if force_cython:
        raise ImportError(
            "Cython HNR module not available. "
            "Compile with: python setup_hnr_cython.py build_ext --inplace"
        )

    # Fallback to scipy.fft optimized version
    try:
        from .hnr_optimized import get_hnr_optimized as get_hnr_scipy
        return get_hnr_scipy(audio, fs, f0_values, frame_shift, n_periods, freq_bands)
    except ImportError:
        pass

    # Final fallback to standard implementation
    from .hnr import get_hnr
    return get_hnr(audio, fs, f0_values, frame_shift, n_periods, freq_bands)


def is_cython_available() -> bool:
    """Check if Cython HNR module is available"""
    return CYTHON_AVAILABLE


def get_hnr_implementation() -> str:
    """
    Get the name of the HNR implementation being used

    Returns:
        'cython', 'scipy', or 'standard'
    """
    if CYTHON_AVAILABLE:
        return 'cython'

    try:
        from .hnr_optimized import NUMBA_AVAILABLE
        return 'scipy'
    except:
        pass

    return 'standard'

"""
Fast pitch tracking with automatic Cython/Numba/Python fallback.
"""

import numpy as np
from scipy import signal

# Try to import accelerated versions
try:
    from .pitch_cython import srh_estimate_pitch_cython, median_filter_f0_cython
    HAS_CYTHON = True
except ImportError:
    HAS_CYTHON = False

from .pitch import srh_pitch_tracking as srh_pitch_tracking_python


def srh_pitch_tracking_fast(wave, fs, f0_min, f0_max, use_cython=True):
    """
    Fast SRH pitch tracking with automatic optimization.

    Uses Cython if available, otherwise falls back to Python.

    Parameters
    ----------
    wave : ndarray
        Speech signal
    fs : int
        Sampling frequency (Hz)
    f0_min : float
        Minimum F0 (Hz)
    f0_max : float
        Maximum F0 (Hz)
    use_cython : bool
        Use Cython if available (default: True)

    Returns
    -------
    f0 : ndarray
        F0 contour (Hz)
    vuv : ndarray
        Voice/unvoiced decision (1=voiced, 0=unvoiced)
    srh_val : ndarray
        SRH correlation values
    """
    wave = np.asarray(wave, dtype=np.float64)

    if HAS_CYTHON and use_cython:
        # Use Cython-accelerated version
        f0, vuv, srh_val = srh_estimate_pitch_cython(
            wave, fs, f0_min, f0_max, K=5, hop_length=160)

        # Median filtering
        f0 = median_filter_f0_cython(f0, vuv, window_size=5)

        return f0, vuv, srh_val
    else:
        # Fall back to Python implementation
        return srh_pitch_tracking_python(wave, fs, f0_min, f0_max)


# Print status on import
if HAS_CYTHON:
    print("✓ Cython available - using optimized pitch tracking")
else:
    print("⚠ Cython not found - using Python fallback (slower)")
    print("  Build with: cd ftrack_python && python setup.py build_ext --inplace")

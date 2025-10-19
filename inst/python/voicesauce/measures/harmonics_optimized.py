"""
Optimized harmonic amplitude extraction using numba JIT compilation
"""

import numpy as np
from typing import Tuple

try:
    from numba import jit, prange
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    # Fallback decorators
    def jit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator
    prange = range


@jit(nopython=True, cache=True, fastmath=True)
def _compute_harmonic_amplitude_fast(segment: np.ndarray, f0: float, fs: int, 
                                      search_range: float = 0.1) -> Tuple[float, float]:
    """
    Fast harmonic amplitude computation using DFT
    Optimized with numba JIT compilation
    
    Args:
        segment: Audio segment
        f0: Target frequency (Hz)
        fs: Sampling rate
        search_range: Search range as fraction of f0
    
    Returns:
        amplitude: Harmonic amplitude in dB
        frequency: Actual frequency found
    """
    n = len(segment)
    
    # Search range
    f_min = max(f0 * (1 - search_range), 1.0)
    f_max = f0 * (1 + search_range)
    
    # Grid search with fine resolution
    n_points = 50
    freq_grid = np.linspace(f_min, f_max, n_points)
    
    max_amp = -np.inf
    best_freq = f0
    
    for freq in freq_grid:
        # Compute DFT at this frequency
        omega = 2 * np.pi * freq / fs
        cos_sum = 0.0
        sin_sum = 0.0
        
        for i in range(n):
            angle = omega * i
            cos_sum += segment[i] * np.cos(angle)
            sin_sum += segment[i] * np.sin(angle)
        
        # Magnitude
        magnitude = np.sqrt(cos_sum**2 + sin_sum**2)
        
        if magnitude > max_amp:
            max_amp = magnitude
            best_freq = freq
    
    # Convert to dB
    amplitude_db = 20 * np.log10(max_amp + 1e-10)
    
    return amplitude_db, best_freq


@jit(nopython=True, parallel=True, cache=True, fastmath=True)
def get_harmonics_fast(audio: np.ndarray, fs: int, f0_values: np.ndarray,
                       sample_shift: int, n_periods: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Fast harmonic extraction using numba parallel processing
    
    Args:
        audio: Full audio signal
        fs: Sampling rate
        f0_values: F0 array
        sample_shift: Samples per frame shift
        n_periods: Number of pitch periods
    
    Returns:
        h1, h2, h4: Harmonic amplitudes
    """
    n_frames = len(f0_values)
    h1 = np.full(n_frames, np.nan)
    h2 = np.full(n_frames, np.nan)
    h4 = np.full(n_frames, np.nan)
    
    # Parallel loop over frames
    for k in prange(n_frames):
        f0 = f0_values[k]
        
        if np.isnan(f0) or f0 <= 0:
            continue
        
        # Calculate segment indices
        ks = k * sample_shift
        n0 = fs / f0
        
        start = int(ks - n_periods/2 * n0)
        end = int(ks + n_periods/2 * n0)
        
        # Check bounds
        if start < 0 or end >= len(audio):
            continue
        
        # Extract segment
        segment = audio[start:end]
        
        if len(segment) == 0:
            continue
        
        # Extract harmonics
        h1[k], _ = _compute_harmonic_amplitude_fast(segment, f0, fs)
        h2[k], _ = _compute_harmonic_amplitude_fast(segment, 2*f0, fs)
        h4[k], _ = _compute_harmonic_amplitude_fast(segment, 4*f0, fs)
    
    return h1, h2, h4


def get_harmonics_optimized(audio: np.ndarray, fs: int, f0_values: np.ndarray,
                            frame_shift: float = 1.0, n_periods: int = 3) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Optimized wrapper for harmonic extraction
    
    Args:
        audio: Full audio signal
        fs: Sampling rate
        f0_values: F0 array
        frame_shift: Frame shift in milliseconds
        n_periods: Number of pitch periods
    
    Returns:
        h1, h2, h4: Harmonic amplitudes
    """
    sample_shift = int(fs * frame_shift / 1000)
    
    if NUMBA_AVAILABLE:
        return get_harmonics_fast(audio, fs, f0_values, sample_shift, n_periods)
    else:
        # Fall back to original implementation
        from .harmonics import get_harmonics
        return get_harmonics(audio, fs, f0_values, frame_shift, n_periods)

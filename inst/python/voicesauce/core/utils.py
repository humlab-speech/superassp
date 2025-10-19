"""
Utility functions for VoiceSauce
"""

import numpy as np
from typing import Optional, Tuple


def align_time_series(target_times: np.ndarray, source_times: np.ndarray,
                     source_values: np.ndarray, 
                     frame_precision: int = 1) -> np.ndarray:
    """
    Align time series data by finding closest time matches
    
    Args:
        target_times: Target time points (ms)
        source_times: Source time points (ms)
        source_values: Source values to interpolate
        frame_precision: Tolerance in frames
    
    Returns:
        aligned_values: Values aligned to target times
    """
    aligned = np.full(len(target_times), np.nan)
    
    for i, t_time in enumerate(target_times):
        # Find closest source time
        diffs = np.abs(source_times - t_time)
        min_idx = np.argmin(diffs)
        
        # Check if within tolerance
        if diffs[min_idx] <= frame_precision:
            aligned[i] = source_values[min_idx]
    
    return aligned


def db_to_linear(db: float) -> float:
    """Convert dB to linear scale"""
    return 10 ** (db / 20)


def linear_to_db(linear: float) -> float:
    """Convert linear to dB scale"""
    if linear <= 0:
        return -np.inf
    return 20 * np.log10(linear)


def hz_to_bark(hz: float) -> float:
    """Convert Hz to Bark scale"""
    return 26.81 * hz / (1960 + hz) - 0.53


def bark_to_hz(bark: float) -> float:
    """Convert Bark to Hz"""
    return 1960 * (bark + 0.53) / (26.28 - bark)


def estimate_bandwidth_hawks_miller(formant: float, formant_number: int) -> float:
    """
    Estimate formant bandwidth using Hawks & Miller formula
    
    Args:
        formant: Formant frequency (Hz)
        formant_number: Formant number (1, 2, 3, ...)
    
    Returns:
        Estimated bandwidth (Hz)
    """
    # Hawks & Miller (1995) formula
    if formant_number == 1:
        return 0.03 * formant + 50
    elif formant_number == 2:
        return 0.03 * formant + 70
    elif formant_number == 3:
        return 0.03 * formant + 110
    elif formant_number == 4:
        return 0.03 * formant + 150
    else:
        return 0.03 * formant + 200


def interpolate_missing_values(values: np.ndarray, method: str = 'linear') -> np.ndarray:
    """
    Interpolate missing (NaN) values
    
    Args:
        values: Input array with NaN values
        method: Interpolation method ('linear', 'nearest', 'zero', 'slinear', 'quadratic', 'cubic')
    
    Returns:
        Interpolated array
    """
    from scipy.interpolate import interp1d
    
    mask = ~np.isnan(values)
    if not mask.any():
        return values
    
    indices = np.arange(len(values))
    valid_indices = indices[mask]
    valid_values = values[mask]
    
    if len(valid_indices) < 2:
        return values
    
    try:
        f = interp1d(valid_indices, valid_values, kind=method, 
                    bounds_error=False, fill_value=np.nan)
        return f(indices)
    except:
        return values


def smooth_values(values: np.ndarray, window_size: int = 5) -> np.ndarray:
    """
    Smooth values using moving average
    
    Args:
        values: Input array
        window_size: Size of smoothing window
    
    Returns:
        Smoothed array
    """
    if window_size < 2:
        return values
    
    # Use convolution for moving average
    kernel = np.ones(window_size) / window_size
    
    # Handle NaN values
    mask = ~np.isnan(values)
    smoothed = values.copy()
    
    if mask.sum() > window_size:
        valid_vals = values[mask]
        smoothed_valid = np.convolve(valid_vals, kernel, mode='same')
        smoothed[mask] = smoothed_valid
    
    return smoothed


def find_peaks_simple(x: np.ndarray, min_distance: int = 1) -> np.ndarray:
    """
    Simple peak finding
    
    Args:
        x: Input signal
        min_distance: Minimum distance between peaks
    
    Returns:
        Peak indices
    """
    from scipy.signal import find_peaks
    peaks, _ = find_peaks(x, distance=min_distance)
    return peaks

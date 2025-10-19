"""
Formant amplitude extraction (A1, A2, A3)
Direct port of func_GetA1A2A3.m
"""

import numpy as np
from typing import Tuple


def get_magnitude_at_frequency(segment: np.ndarray, freq: float, fs: int, fft_len: int = 8192) -> Tuple[float, float]:
    """
    Get maximal spectral magnitude around a frequency
    Port of MATLAB ana_GetMagnitudeMax
    
    Args:
        segment: Audio segment
        freq: Target frequency (Hz)
        fs: Sampling rate
        fft_len: FFT length
    
    Returns:
        magnitude: Maximum magnitude in dB
        freq_max: Frequency of maximum
    """
    if np.isnan(freq) or freq <= 0:
        return np.nan, np.nan
    
    # FFT
    spectrum = np.fft.fft(segment, fft_len)
    magnitude_db = 20 * np.log10(np.abs(spectrum[:fft_len//2]) + 1e-10)
    
    # Frequency resolution
    freq_step = fs / fft_len
    
    # Search range: ±10% around target frequency
    low_freq = max(freq - 0.1 * freq, 0)
    high_freq = min(freq + 0.1 * freq, fs/2 - freq_step)
    
    # Convert to bin indices
    low_bin = int(low_freq / freq_step)
    high_bin = int(high_freq / freq_step) + 1
    
    if low_bin >= len(magnitude_db) or high_bin >= len(magnitude_db):
        return np.nan, np.nan
    
    # Find maximum in range
    search_region = magnitude_db[low_bin:high_bin]
    if len(search_region) == 0:
        return np.nan, np.nan
    
    max_idx = np.argmax(search_region)
    max_magnitude = search_region[max_idx]
    max_freq = (low_bin + max_idx) * freq_step
    
    return max_magnitude, max_freq


def get_formant_amplitudes(audio: np.ndarray, fs: int, 
                          f0_values: np.ndarray,
                          f1_values: np.ndarray,
                          f2_values: np.ndarray,
                          f3_values: np.ndarray,
                          frame_shift: float = 1.0,
                          n_periods: int = 3) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate A1, A2, A3 for entire signal
    Port of func_GetA1A2A3.m
    
    Args:
        audio: Full audio signal
        fs: Sampling rate
        f0_values: F0 array
        f1_values, f2_values, f3_values: Formant arrays
        frame_shift: Frame shift in milliseconds
        n_periods: Number of pitch periods to extract
    
    Returns:
        a1, a2, a3: Formant amplitudes (dB) for each frame
    """
    sample_shift = int(fs * frame_shift / 1000)
    
    a1 = np.full(len(f0_values), np.nan)
    a2 = np.full(len(f0_values), np.nan)
    a3 = np.full(len(f0_values), np.nan)
    
    for k in range(len(f0_values)):
        f0 = f0_values[k]
        f1 = f1_values[k] if k < len(f1_values) else np.nan
        f2 = f2_values[k] if k < len(f2_values) else np.nan
        f3 = f3_values[k] if k < len(f3_values) else np.nan
        
        if np.isnan(f0) or f0 <= 0:
            continue
        
        if np.isnan(f1) or np.isnan(f2) or np.isnan(f3):
            continue
        
        # Calculate sample indices for segment
        ks = int(k * sample_shift)
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
        
        # Extract formant amplitudes
        try:
            a1[k], _ = get_magnitude_at_frequency(segment, f1, fs)
            a2[k], _ = get_magnitude_at_frequency(segment, f2, fs)
            a3[k], _ = get_magnitude_at_frequency(segment, f3, fs)
        except:
            continue
    
    return a1, a2, a3


def compute_harmonic_formant_differences(h1: np.ndarray, 
                                        a1: np.ndarray, 
                                        a2: np.ndarray, 
                                        a3: np.ndarray) -> dict:
    """
    Compute differences between H1 and formant amplitudes
    
    Args:
        h1: H1 amplitudes
        a1, a2, a3: Formant amplitudes
    
    Returns:
        Dictionary with H1-A1, H1-A2, H1-A3 differences
    """
    return {
        'H1_A1': h1 - a1,
        'H1_A2': h1 - a2,
        'H1_A3': h1 - a3
    }

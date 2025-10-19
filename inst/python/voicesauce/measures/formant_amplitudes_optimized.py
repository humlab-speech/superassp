"""
Optimized formant amplitude extraction using vectorized operations
"""

import numpy as np
from scipy.fft import fft
from typing import Tuple


def get_formant_amplitudes_optimized(audio: np.ndarray, fs: int,
                                     f0_values: np.ndarray,
                                     f1_values: np.ndarray,
                                     f2_values: np.ndarray,
                                     f3_values: np.ndarray,
                                     frame_shift: float = 1.0,
                                     n_periods: int = 3,
                                     fft_len: int = 8192) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Optimized formant amplitude extraction using scipy.fft
    
    Uses scipy.fft (faster than numpy.fft) and vectorized operations
    
    Args:
        audio: Full audio signal
        fs: Sampling rate
        f0_values: F0 array
        f1_values, f2_values, f3_values: Formant arrays
        frame_shift: Frame shift in milliseconds
        n_periods: Number of pitch periods
        fft_len: FFT length
    
    Returns:
        a1, a2, a3: Formant amplitudes
    """
    sample_shift = int(fs * frame_shift / 1000)
    freq_step = fs / fft_len
    
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
        
        ks = int(k * sample_shift)
        n0 = fs / f0
        
        start = int(ks - n_periods/2 * n0)
        end = int(ks + n_periods/2 * n0)
        
        if start < 0 or end >= len(audio):
            continue
        
        segment = audio[start:end]
        if len(segment) == 0:
            continue
        
        # Compute spectrum using scipy.fft (faster)
        spectrum = fft(segment, fft_len)
        magnitude_db = 20 * np.log10(np.abs(spectrum[:fft_len//2]) + 1e-10)
        
        # Vectorized extraction for A1
        low_bin = int((f1 - 0.1 * f1) / freq_step)
        high_bin = int((f1 + 0.1 * f1) / freq_step) + 1
        
        if 0 <= low_bin < len(magnitude_db) and low_bin < high_bin <= len(magnitude_db):
            a1[k] = np.max(magnitude_db[low_bin:high_bin])
        
        # Vectorized extraction for A2
        low_bin = int((f2 - 0.1 * f2) / freq_step)
        high_bin = int((f2 + 0.1 * f2) / freq_step) + 1
        
        if 0 <= low_bin < len(magnitude_db) and low_bin < high_bin <= len(magnitude_db):
            a2[k] = np.max(magnitude_db[low_bin:high_bin])
        
        # Vectorized extraction for A3
        low_bin = int((f3 - 0.1 * f3) / freq_step)
        high_bin = int((f3 + 0.1 * f3) / freq_step) + 1
        
        if 0 <= low_bin < len(magnitude_db) and low_bin < high_bin <= len(magnitude_db):
            a3[k] = np.max(magnitude_db[low_bin:high_bin])
    
    return a1, a2, a3

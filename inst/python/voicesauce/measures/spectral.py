"""
Spectral measures (2K, 5K, H4-2K)
Direct port of func_Get2K.m, func_Get5K.m, func_Get2K5K.m, func_GetH42K.m
"""

import numpy as np
from typing import Tuple


def get_spectral_measure_at_freq(audio: np.ndarray, fs: int, f0_values: np.ndarray,
                                 target_freq: float,
                                 frame_shift: float = 1.0, 
                                 n_periods: int = 3,
                                 fft_len: int = 8192) -> np.ndarray:
    """
    Get spectral energy at a target frequency
    Generic function for 2K, 5K, etc.
    
    Args:
        audio: Full audio signal
        fs: Sampling rate
        f0_values: F0 array
        target_freq: Target frequency (Hz)
        frame_shift: Frame shift in milliseconds
        n_periods: Number of pitch periods
        fft_len: FFT length
    
    Returns:
        Spectral values (dB) for each frame
    """
    sample_shift = int(fs * frame_shift / 1000)
    spectral = np.full(len(f0_values), np.nan)
    
    for k, f0 in enumerate(f0_values):
        if np.isnan(f0) or f0 <= 0:
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
        
        # FFT
        spectrum = np.fft.fft(segment, fft_len)
        magnitude_db = 20 * np.log10(np.abs(spectrum[:fft_len//2]) + 1e-10)
        
        # Get value at target frequency
        freq_step = fs / fft_len
        target_bin = int(target_freq / freq_step)
        
        if target_bin < len(magnitude_db):
            spectral[k] = magnitude_db[target_bin]
    
    return spectral


def get_2k(audio: np.ndarray, fs: int, f0_values: np.ndarray,
          frame_shift: float = 1.0, n_periods: int = 3) -> np.ndarray:
    """Get spectral energy at 2kHz"""
    return get_spectral_measure_at_freq(audio, fs, f0_values, 2000.0, 
                                       frame_shift, n_periods)


def get_5k(audio: np.ndarray, fs: int, f0_values: np.ndarray,
          frame_shift: float = 1.0, n_periods: int = 3) -> np.ndarray:
    """Get spectral energy at 5kHz"""
    return get_spectral_measure_at_freq(audio, fs, f0_values, 5000.0, 
                                       frame_shift, n_periods)


def get_2k5k(two_k: np.ndarray, five_k: np.ndarray) -> np.ndarray:
    """Calculate 2K-5K spectral tilt"""
    return two_k - five_k


def get_h42k(h4: np.ndarray, two_k: np.ndarray) -> np.ndarray:
    """Calculate H4-2K"""
    return h4 - two_k

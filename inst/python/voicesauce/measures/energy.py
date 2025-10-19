"""
Energy calculation
Direct port of func_GetEnergy.m
"""

import numpy as np


def get_energy(audio: np.ndarray, fs: int, f0_values: np.ndarray,
              frame_shift: float = 1.0, n_periods: int = 5) -> np.ndarray:
    """
    Calculate energy normalized by F0
    Port of func_GetEnergy.m
    
    Args:
        audio: Full audio signal
        fs: Sampling rate
        f0_values: F0 array
        frame_shift: Frame shift in milliseconds
        n_periods: Number of pitch periods to extract
    
    Returns:
        energy: Energy values for each frame
    """
    sample_shift = int(fs * frame_shift / 1000)
    energy = np.full(len(f0_values), np.nan)
    
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
        if len(segment) > 0:
            energy[k] = np.sum(segment ** 2)
    
    return energy

"""
Harmonic amplitude extraction (H1, H2, H4)
Direct port of func_GetH1_H2_H4.m
"""

import numpy as np
from scipy.optimize import minimize_scalar
from typing import Tuple


def extract_harmonic_amplitude(segment: np.ndarray, f0: float, fs: int) -> Tuple[float, float]:
    """
    Extract harmonic amplitude using frequency-domain optimization
    Port of MATLAB func_GetHarmonics
    
    Args:
        segment: Audio segment
        f0: Target frequency (Hz)
        fs: Sampling rate
    
    Returns:
        amplitude: Harmonic amplitude in dB
        frequency: Actual frequency found
    """
    def objective(f):
        """Objective function to minimize (negative amplitude)"""
        n = np.arange(len(segment))
        v = np.exp(-1j * 2 * np.pi * f * n / fs)
        amplitude_db = 20 * np.log10(np.abs(np.dot(segment, v)) + 1e-10)
        return -amplitude_db
    
    # Search range: ±10% of f0
    df_range = 0.1 * f0
    f_min = max(f0 - df_range, 1.0)
    f_max = f0 + df_range
    
    try:
        # Bounded optimization
        result = minimize_scalar(
            objective,
            bounds=(f_min, f_max),
            method='bounded'
        )
        
        amplitude = -result.fun  # Convert back to positive
        frequency = result.x
        
        return amplitude, frequency
    except:
        return np.nan, np.nan


def get_harmonics(audio: np.ndarray, fs: int, f0_values: np.ndarray,
                 frame_shift: float = 1.0, n_periods: int = 3) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate H1, H2, H4 for entire signal
    Port of func_GetH1_H2_H4.m
    
    Args:
        audio: Full audio signal
        fs: Sampling rate
        f0_values: F0 array (one value per frame)
        frame_shift: Frame shift in milliseconds
        n_periods: Number of pitch periods to extract
    
    Returns:
        h1: H1 values (dB) for each frame
        h2: H2 values (dB) for each frame
        h4: H4 values (dB) for each frame
    """
    sample_shift = int(fs * frame_shift / 1000)
    
    h1 = np.full(len(f0_values), np.nan)
    h2 = np.full(len(f0_values), np.nan)
    h4 = np.full(len(f0_values), np.nan)
    
    for k, f0 in enumerate(f0_values):
        if np.isnan(f0) or f0 <= 0:
            continue
        
        # Calculate sample indices for segment
        ks = int(k * sample_shift)
        n0 = fs / f0  # Samples per period
        
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
        try:
            h1[k], _ = extract_harmonic_amplitude(segment, f0, fs)
            h2[k], _ = extract_harmonic_amplitude(segment, 2*f0, fs)
            h4[k], _ = extract_harmonic_amplitude(segment, 4*f0, fs)
        except:
            continue
    
    return h1, h2, h4


def compute_harmonic_differences(h1: np.ndarray, h2: np.ndarray, h4: np.ndarray) -> dict:
    """
    Compute harmonic amplitude differences
    
    Args:
        h1, h2, h4: Harmonic amplitudes
    
    Returns:
        Dictionary with H1-H2 and H2-H4 differences
    """
    return {
        'H1_H2': h1 - h2,
        'H2_H4': h2 - h4
    }

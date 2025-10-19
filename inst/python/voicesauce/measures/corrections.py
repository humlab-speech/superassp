"""
Iseli-Alwan formant corrections
Direct port of func_correct_iseli_z.m
"""

import numpy as np
from typing import Union


def correct_iseli_alwan(f: Union[float, np.ndarray], 
                       fx: Union[float, np.ndarray], 
                       bx: Union[float, np.ndarray], 
                       fs: int) -> Union[float, np.ndarray]:
    """
    Iseli-Alwan formant correction
    Direct port of func_correct_iseli_z.m
    
    Inverse-filtering of vocal tract resonance (formant). 
    Correction is in dB, and dependent on fs.
    
    Args:
        f: Frequency to correct (e.g., F0, 2*F0)
        fx: Formant frequency
        bx: Formant bandwidth
        fs: Sampling rate
    
    Returns:
        Correction value in dB
    
    References:
        Iseli, M., & Alwan, A. (2004). An improved correction formula for 
        the estimation of harmonic magnitudes and its application to open 
        quotient estimation. ICASSP.
    """
    r = np.exp(-np.pi * bx / fs)
    omega_x = 2 * np.pi * fx / fs
    omega = 2 * np.pi * f / fs
    
    a = r**2 + 1 - 2*r*np.cos(omega_x + omega)
    b = r**2 + 1 - 2*r*np.cos(omega_x - omega)
    num = r**2 + 1 - 2*r*np.cos(omega_x)
    
    corr = -10*(np.log10(a) + np.log10(b)) + 20*np.log10(num)
    
    return corr


def apply_corrections(h1: np.ndarray, h2: np.ndarray, h4: np.ndarray,
                     a1: np.ndarray, a2: np.ndarray, a3: np.ndarray,
                     f0: np.ndarray, 
                     f1: np.ndarray, f2: np.ndarray, f3: np.ndarray,
                     b1: np.ndarray, b2: np.ndarray, b3: np.ndarray,
                     fs: int) -> dict:
    """
    Apply Iseli-Alwan corrections to all measures
    
    Args:
        h1, h2, h4: Harmonic amplitudes
        a1, a2, a3: Formant amplitudes
        f0: Fundamental frequency
        f1, f2, f3: Formant frequencies
        b1, b2, b3: Formant bandwidths
        fs: Sampling rate
    
    Returns:
        Dictionary with corrected measures (suffix 'c')
    """
    # Initialize corrected values
    h1c = h1.copy()
    h2c = h2.copy()
    h4c = h4.copy()
    a1c = a1.copy()
    a2c = a2.copy()
    a3c = a3.copy()
    
    # Apply corrections frame by frame
    for k in range(len(f0)):
        if np.isnan(f0[k]) or f0[k] <= 0:
            continue
        
        if np.isnan(f1[k]) or np.isnan(f2[k]) or np.isnan(f3[k]):
            continue
        
        if np.isnan(b1[k]) or np.isnan(b2[k]) or np.isnan(b3[k]):
            continue
        
        # Correct H1 for F1, F2, F3
        if not np.isnan(h1[k]):
            h1c[k] = h1[k]
            h1c[k] -= correct_iseli_alwan(f0[k], f1[k], b1[k], fs)
            h1c[k] -= correct_iseli_alwan(f0[k], f2[k], b2[k], fs)
            h1c[k] -= correct_iseli_alwan(f0[k], f3[k], b3[k], fs)
        
        # Correct H2 for F1, F2, F3
        if not np.isnan(h2[k]):
            h2c[k] = h2[k]
            h2c[k] -= correct_iseli_alwan(2*f0[k], f1[k], b1[k], fs)
            h2c[k] -= correct_iseli_alwan(2*f0[k], f2[k], b2[k], fs)
            h2c[k] -= correct_iseli_alwan(2*f0[k], f3[k], b3[k], fs)
        
        # Correct H4 for F1, F2, F3
        if not np.isnan(h4[k]):
            h4c[k] = h4[k]
            h4c[k] -= correct_iseli_alwan(4*f0[k], f1[k], b1[k], fs)
            h4c[k] -= correct_iseli_alwan(4*f0[k], f2[k], b2[k], fs)
            h4c[k] -= correct_iseli_alwan(4*f0[k], f3[k], b3[k], fs)
        
        # Correct A1 (only for F1)
        if not np.isnan(a1[k]):
            a1c[k] = a1[k]
            a1c[k] -= correct_iseli_alwan(f1[k], f1[k], b1[k], fs)
        
        # Correct A2 (only for F2)
        if not np.isnan(a2[k]):
            a2c[k] = a2[k]
            a2c[k] -= correct_iseli_alwan(f2[k], f2[k], b2[k], fs)
        
        # Correct A3 (only for F3)
        if not np.isnan(a3[k]):
            a3c[k] = a3[k]
            a3c[k] -= correct_iseli_alwan(f3[k], f3[k], b3[k], fs)
    
    # Calculate corrected differences
    return {
        'H1c': h1c,
        'H2c': h2c,
        'H4c': h4c,
        'A1c': a1c,
        'A2c': a2c,
        'A3c': a3c,
        'H1c_H2c': h1c - h2c,
        'H2c_H4c': h2c - h4c,
        'H1c_A1c': h1c - a1c,
        'H1c_A2c': h1c - a2c,
        'H1c_A3c': h1c - a3c
    }

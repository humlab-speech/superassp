"""
Teager-Kaiser Energy Operator (TKEO)

Ported from TKEO function in voice_analysis_redux.m (lines 598-611)
"""

import numpy as np


def compute_tkeo(x):
    """
    Compute Teager-Kaiser Energy Operator
    
    The TKEO is defined as:
    E[n] = x[n]^2 - x[n-1] * x[n+1]
    
    Parameters:
    -----------
    x : ndarray
        Input signal
        
    Returns:
    --------
    energy : ndarray
        TKEO values (same length as input)
    """
    x = np.asarray(x).ravel()
    energy = np.zeros_like(x, dtype=np.float64)
    
    # Boundary conditions
    energy[0] = x[0]**2
    
    # Main TKEO formula
    for n in range(1, len(x) - 1):
        energy[n] = x[n]**2 - x[n-1] * x[n+1]
    
    # Boundary conditions
    energy[-1] = x[-1]**2
    
    return energy


def compute_tkeo_vectorized(x):
    """
    Vectorized version of TKEO (faster for large arrays)
    
    Parameters:
    -----------
    x : ndarray
        Input signal
        
    Returns:
    --------
    energy : ndarray
        TKEO values
    """
    x = np.asarray(x).ravel()
    energy = np.zeros_like(x, dtype=np.float64)
    
    energy[0] = x[0]**2
    energy[1:-1] = x[1:-1]**2 - x[:-2] * x[2:]
    energy[-1] = x[-1]**2
    
    return energy

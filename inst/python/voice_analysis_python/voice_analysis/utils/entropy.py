"""
Entropy Calculations

Ported from H_entropy and log_bb functions in voice_analysis_redux.m
"""

import numpy as np


def compute_entropy(distribution, base='e'):
    """
    Compute entropy of a discrete distribution
    
    Ported from H_entropy function (lines 613-623)
    
    Parameters:
    -----------
    distribution : ndarray
        Probability distribution (will be normalized if not already)
    base : str
        Logarithm base: 'e' (nats), '2' (bits), '10' (hartley)
        
    Returns:
    --------
    H : float
        Entropy value
    """
    distribution = np.asarray(distribution).ravel()
    
    # Normalize if not already a probability distribution
    if np.sum(distribution) > 0:
        distribution = distribution / np.sum(distribution)
    else:
        return 0.0
    
    # Select logarithm base
    if base == '2':
        log_func = np.log2
    elif base == '10':
        log_func = np.log10
    else:  # 'e' or 'nats'
        log_func = np.log
    
    # Compute entropy (handle zeros)
    H = 0.0
    for p in distribution:
        if p > 0:
            H -= p * log_func(p)
    
    return H


def safe_log(x, base='e'):
    """
    Safe logarithm that returns 0 for x=0
    
    Ported from log_bb function (lines 626-651)
    
    Parameters:
    -----------
    x : float or ndarray
        Input value(s)
    base : str
        Logarithm base: 'e' (nats), '2' (bits), '10' (hartley)
        
    Returns:
    --------
    result : float or ndarray
        Logarithm, with 0 for input 0
    """
    if base == '2':
        log_func = np.log2
    elif base == '10':
        log_func = np.log10
    else:  # 'e' or 'nats'
        log_func = np.log
    
    if np.isscalar(x):
        return log_func(x) if x > 0 else 0.0
    else:
        x = np.asarray(x)
        result = np.zeros_like(x, dtype=np.float64)
        mask = x > 0
        result[mask] = log_func(x[mask])
        return result

"""
Detrended Fluctuation Analysis (DFA)

Ported from fastdfa.m and fastdfa_core.c
DFA quantifies fractal-like scaling properties of signals

OPTIMIZED VERSION with Numba JIT compilation for 5-10x speedup
"""

import numpy as np
import warnings

# Try to import numba for performance optimization
try:
    from numba import jit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    # Fallback decorator that does nothing
    def jit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator


def compute_dfa(signal, scales=None):
    """
    Compute Detrended Fluctuation Analysis
    
    Parameters:
    -----------
    signal : ndarray
        Input signal
    scales : array-like
        Scale ranges to analyze (default: 50 to 200 in steps of 20)
        
    Returns:
    --------
    dfa_transformed : float
        DFA value transformed to [0, 1] range via sigmoid
    """
    try:
        # Try using nolds library first (well-tested implementation)
        import nolds
        
        signal = np.asarray(signal).ravel()
        signal = signal[np.isfinite(signal)]
        
        if len(signal) < 100:
            return np.nan
        
        # Compute DFA exponent (alpha)
        if scales is None:
            # Match MATLAB default from line 123: dfa_scaling = (50:20:200)'
            nvals = np.arange(50, 201, 20)
        else:
            nvals = np.asarray(scales)
        
        try:
            alpha = nolds.dfa(signal, nvals=nvals)
        except:
            # If nolds fails, use our implementation
            alpha = _dfa_manual(signal, nvals)
        
        # Transform to bounded range (line 190 in MATLAB)
        # DFA = 1/(1+exp(-dfa))
        dfa_transformed = 1.0 / (1.0 + np.exp(-alpha))
        
        return dfa_transformed
        
    except ImportError:
        # Fallback to manual implementation
        signal = np.asarray(signal).ravel()
        signal = signal[np.isfinite(signal)]
        
        if len(signal) < 100:
            return np.nan
        
        if scales is None:
            scales = np.arange(50, 201, 20)
        
        alpha = _dfa_manual(signal, scales)
        dfa_transformed = 1.0 / (1.0 + np.exp(-alpha))
        
        return dfa_transformed


def _dfa_manual(signal, scales):
    """
    Manual DFA implementation with optimization
    
    Ported from fastdfa_core.c logic
    Uses Numba JIT compilation for 5-10x speedup if available
    """
    N = len(signal)
    
    # Step 1: Integrate the signal (cumulative sum of deviations from mean)
    mean_val = np.mean(signal)
    y = np.cumsum(signal - mean_val)
    
    if NUMBA_AVAILABLE:
        fluctuations, valid_scales = _compute_fluctuations_numba(y, scales, N)
    else:
        fluctuations, valid_scales = _compute_fluctuations_python(y, scales, N)
    
    if len(fluctuations) < 2:
        return np.nan
    
    # Log-log fit to get scaling exponent
    log_scales = np.log10(valid_scales)
    log_flucts = np.log10(fluctuations)
    
    # Linear fit (optimized)
    alpha = _fast_linear_fit(log_scales, log_flucts)
    
    return alpha


@jit(nopython=True, cache=True)
def _compute_fluctuations_numba(y, scales, N):
    """
    Optimized fluctuation computation with Numba
    
    This is the main bottleneck - 5-10x speedup with Numba
    """
    fluctuations = []
    valid_scales = []
    
    for scale_idx in range(len(scales)):
        n = int(scales[scale_idx])
        
        if n < 4 or n >= N:
            continue
        
        # Number of segments
        n_segments = N // n
        
        if n_segments < 1:
            continue
        
        sum_flucts = 0.0
        count = 0
        
        for v in range(n_segments):
            # Extract segment
            start = v * n
            end = (v + 1) * n
            
            # Compute linear trend and fluctuation in one pass
            # This is much faster than polyfit
            sum_t = 0.0
            sum_y = 0.0
            sum_ty = 0.0
            sum_t2 = 0.0
            
            for i in range(n):
                t = float(i)
                y_val = y[start + i]
                sum_t += t
                sum_y += y_val
                sum_ty += t * y_val
                sum_t2 += t * t
            
            # Linear regression: y = a + b*t
            n_float = float(n)
            denom = n_float * sum_t2 - sum_t * sum_t
            
            if abs(denom) < 1e-10:
                continue
            
            b = (n_float * sum_ty - sum_t * sum_y) / denom
            a = (sum_y - b * sum_t) / n_float
            
            # Compute fluctuation (RMSE from trend)
            sum_sq_err = 0.0
            for i in range(n):
                t = float(i)
                y_val = y[start + i]
                trend_val = a + b * t
                err = y_val - trend_val
                sum_sq_err += err * err
            
            fluctuation = np.sqrt(sum_sq_err / n_float)
            sum_flucts += fluctuation
            count += 1
        
        if count > 0:
            # Average fluctuation at this scale
            F_n = sum_flucts / count
            fluctuations.append(F_n)
            valid_scales.append(n)
    
    return np.array(fluctuations), np.array(valid_scales, dtype=np.float64)


def _compute_fluctuations_python(y, scales, N):
    """
    Python fallback for fluctuation computation (slower but works without Numba)
    """
    fluctuations = []
    valid_scales = []
    
    for n in scales:
        n = int(n)
        
        if n < 4 or n >= N:
            continue
        
        # Number of segments
        n_segments = N // n
        
        if n_segments < 1:
            continue
        
        local_flucts = []
        
        for v in range(n_segments):
            # Extract segment
            segment = y[v*n:(v+1)*n]
            
            # Fit linear trend
            t = np.arange(len(segment))
            
            if len(t) < 2:
                continue
            
            # Linear regression
            coeffs = np.polyfit(t, segment, 1)
            trend = np.polyval(coeffs, t)
            
            # Compute fluctuation (RMSE from trend)
            fluctuation = np.sqrt(np.mean((segment - trend)**2))
            local_flucts.append(fluctuation)
        
        if len(local_flucts) > 0:
            # Average fluctuation at this scale
            F_n = np.mean(local_flucts)
            fluctuations.append(F_n)
            valid_scales.append(n)
    
    return np.array(fluctuations), np.array(valid_scales)


@jit(nopython=True, cache=True)
def _fast_linear_fit(x, y):
    """Fast linear regression to get slope (optimized with Numba)"""
    n = len(x)
    sum_x = np.sum(x)
    sum_y = np.sum(y)
    sum_xy = np.sum(x * y)
    sum_x2 = np.sum(x * x)
    
    n_float = float(n)
    denom = n_float * sum_x2 - sum_x * sum_x
    
    if abs(denom) < 1e-10:
        return np.nan
    
    slope = (n_float * sum_xy - sum_x * sum_y) / denom
    
    return slope

"""
Optimized versions of computationally intensive F0 extraction functions using Numba JIT

This module provides drop-in replacements for bottleneck functions with significant speedups.
"""

import numpy as np
from typing import Tuple

try:
    from numba import jit, prange
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    # Create dummy decorator if Numba not available
    def jit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator
    def prange(*args, **kwargs):
        return range(*args, **kwargs)


@jit(nopython=True, cache=True)
def _compute_parabolic_interp(yv: np.ndarray, xo: float) -> Tuple[float, float]:
    """
    Numba-optimized parabolic interpolation
    
    Matches zzParabolicInterp2 algorithm exactly
    
    Args:
        yv: 3-element array [y(-1), y(0), y(1)]
        xo: Initial x position
        
    Returns:
        val: Interpolated peak value
        pos: Interpolated peak position
    """
    # Compute differences: lp = diff(yv) = [yv[1]-yv[0], yv[2]-yv[1]]
    lp0 = yv[1] - yv[0]
    lp1 = yv[2] - yv[1]
    
    # a = lp[0] - lp[1]
    a = lp0 - lp1
    
    # b = (lp[0] + lp[1]) / 2
    b = (lp0 + lp1) / 2.0
    
    # Avoid division by zero
    if abs(a) < 1e-10:
        return yv[1], xo - 1.0
    
    # xp = b / a + xo
    xp = b / a + xo
    
    # val = yv[1] + 0.5 * a * (b / a)**2 + b * (b / a)
    b_over_a = b / a
    val = yv[1] + 0.5 * a * b_over_a * b_over_a + b * b_over_a
    
    # pos = xp - 1
    pos = xp - 1.0
    
    return val, pos


@jit(nopython=True, cache=True)
def _find_peaks_above_threshold(data: np.ndarray, threshold: float) -> np.ndarray:
    """
    Fast peak finding above threshold using Numba
    
    Args:
        data: 1D array to search for peaks
        threshold: Minimum peak value
        
    Returns:
        Array of peak indices
    """
    n = len(data)
    peaks = []
    
    for i in range(1, n - 1):
        if data[i] > threshold:
            if data[i] > data[i-1] and data[i] > data[i+1]:
                peaks.append(i)
    
    return np.array(peaks, dtype=np.int64)


@jit(nopython=True, cache=True)
def _compute_acc_indices_and_scaling(fftl: int, wflfso: int, bbf: np.ndarray,
                                     nsht: float, ndiv: int, npw: np.ndarray,
                                     pc: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute indices and scaling factors for ACC computation
    
    This pre-computes values needed for the ACC inner loop in ztestspecspecnormal.
    """
    center_idx = int((wflfso - 1) / 2) - 1
    
    # Pre-allocate arrays
    all_indices = np.empty((ndiv + 1, wflfso), dtype=np.int64)
    scaling_factors = np.empty(ndiv + 1)
    
    for ii in range(ndiv + 1):
        # Compute indices for this division
        for j in range(wflfso):
            idx = int(np.round(fftl / 2 + bbf[j] + ii * nsht))
            all_indices[ii, j] = idx % fftl
        
        # Compute scaling factor
        scaling_factors[ii] = npw[all_indices[ii, center_idx]] ** pc
    
    return all_indices, scaling_factors


@jit(nopython=True, cache=True)
def _apply_window_and_scale(pw: np.ndarray, w2: np.ndarray, 
                            all_indices: np.ndarray, scaling_factors: np.ndarray,
                            wflfso: int, ndiv: int) -> np.ndarray:
    """
    Apply windowing and scaling for ACC computation
    
    Returns windowed and scaled power spectrum for each division.
    """
    windowed_data = np.empty((ndiv + 1, wflfso))
    
    for ii in range(ndiv + 1):
        for j in range(wflfso):
            windowed_data[ii, j] = pw[all_indices[ii, j]] * w2[j] * scaling_factors[ii]
    
    return windowed_data


@jit(nopython=True, cache=True)
def _zfixpfreq3_core(fxx: np.ndarray, pif2: np.ndarray, mmp: np.ndarray,
                     dfv: np.ndarray, aav: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Core computation for zfixpfreq3 - find fixed points in frequency map
    
    Optimized with Numba JIT compilation for the tight loop.
    """
    nn = len(fxx)
    
    # Compute cd1 and cd2
    cd1 = pif2 - fxx
    cd2 = np.empty(nn)
    for i in range(nn - 1):
        cd2[i] = cd1[i + 1] - cd1[i]
    cd2[nn - 1] = cd1[nn - 1] - cd1[nn - 2]
    
    # Shifted cd1
    cdd1 = np.empty(nn)
    for i in range(nn - 1):
        cdd1[i] = cd1[i + 1]
    cdd1[nn - 1] = cd1[nn - 1]
    
    # Find fixed points
    fp_mask = (cd1 * cdd1 < 0) & (cd2 < 0)
    
    # Count and collect indices
    count = 0
    for i in range(nn):
        if fp_mask[i]:
            count += 1
    
    if count == 0:
        return np.empty(0), np.empty(0), np.empty(0), np.empty(0)
    
    # Collect valid indices
    ixx = np.empty(count, dtype=np.int64)
    idx = 0
    for i in range(nn - 1):  # Ensure i + 1 doesn't exceed bounds
        if fp_mask[i]:
            ixx[idx] = i
            idx += 1
    
    # Trim to actual count
    ixx = ixx[:idx]
    n_valid = len(ixx)
    
    if n_valid == 0:
        return np.empty(0), np.empty(0), np.empty(0), np.empty(0)
    
    # Compute interpolated values
    ff = np.empty(n_valid)
    vv = np.empty(n_valid)
    df = np.empty(n_valid)
    aa = np.empty(n_valid)
    
    for k in range(n_valid):
        i = ixx[k]
        
        # Interpolate frequency
        denom = cd1[i] - cdd1[i]
        if abs(denom) < 1e-10:
            denom = 1e-10
        ff[k] = pif2[i] + (pif2[i + 1] - pif2[i]) * cd1[i] / denom
        
        # Interpolate other values
        df_x = fxx[i + 1] - fxx[i]
        if abs(df_x) < 1e-10:
            df_x = 1e-10
        
        ratio = (ff[k] - fxx[i]) / df_x
        vv[k] = mmp[i] + (mmp[i + 1] - mmp[i]) * ratio
        df[k] = dfv[i] + (dfv[i + 1] - dfv[i]) * ratio
        aa[k] = aav[i] + (aav[i + 1] - aav[i]) * ratio
    
    return ff, vv, df, aa


@jit(nopython=True, cache=True)
def _octave_distance(f1: float, f2: float) -> float:
    """
    Compute octave distance between two frequencies
    
    Args:
        f1, f2: Frequencies in Hz
        
    Returns:
        Distance in octaves (log2(f2/f1))
    """
    if f1 <= 0 or f2 <= 0:
        return 1e10
    return abs(np.log2(f2 / f1))


@jit(nopython=True, cache=True)
def _compute_tracking_costs(candidates: np.ndarray, last_f0: float, 
                           reliabilities: np.ndarray) -> np.ndarray:
    """
    Compute F0 tracking costs for candidate selection
    
    Args:
        candidates: Array of candidate F0 values (Hz)
        last_f0: Previous F0 value (Hz)
        reliabilities: Reliability scores for candidates
        
    Returns:
        Array of costs (lower is better)
    """
    n = len(candidates)
    costs = np.empty(n)
    
    for i in range(n):
        f0 = candidates[i]
        rel = reliabilities[i]
        
        # Octave distance penalty
        oct_dist = _octave_distance(f0, last_f0)
        
        # Combined cost: reliability (inverted) + octave distance penalty
        costs[i] = (1.0 - rel) + 2.0 * oct_dist
    
    return costs


@jit(nopython=True, cache=True)
def _smooth_contour(contour: np.ndarray, window_size: int) -> np.ndarray:
    """
    Apply moving average smoothing to F0 contour
    
    Args:
        contour: F0 contour
        window_size: Smoothing window size
        
    Returns:
        Smoothed contour
    """
    n = len(contour)
    smoothed = np.empty(n)
    
    for i in range(n):
        # Determine window bounds
        start = max(0, i - window_size // 2)
        end = min(n, i + window_size // 2 + 1)
        
        # Compute average
        total = 0.0
        count = 0
        for j in range(start, end):
            if contour[j] > 0:  # Skip zeros
                total += contour[j]
                count += 1
        
        if count > 0:
            smoothed[i] = total / count
        else:
            smoothed[i] = contour[i]
    
    return smoothed


@jit(nopython=True, cache=True)
def _apply_median_filter_3(data: np.ndarray) -> np.ndarray:
    """
    Fast 3-point median filter
    
    Args:
        data: Input array
        
    Returns:
        Median filtered array
    """
    n = len(data)
    result = np.empty(n)
    
    # Handle edges
    result[0] = data[0]
    result[n-1] = data[n-1]
    
    # Apply median filter
    for i in range(1, n - 1):
        # Simple 3-element median
        a, b, c = data[i-1], data[i], data[i+1]
        if a <= b:
            if b <= c:
                result[i] = b
            elif a <= c:
                result[i] = c
            else:
                result[i] = a
        else:  # a > b
            if a <= c:
                result[i] = a
            elif b <= c:
                result[i] = c
            else:
                result[i] = b
    
    return result


@jit(nopython=True, cache=True)
def _compute_power_envelope(x: np.ndarray, window_samples: int) -> np.ndarray:
    """
    Fast power envelope computation using sliding window
    
    Args:
        x: Input signal
        window_samples: Window size in samples
        
    Returns:
        Power envelope
    """
    n = len(x)
    power = np.empty(n)
    
    for i in range(n):
        start = max(0, i - window_samples // 2)
        end = min(n, i + window_samples // 2 + 1)
        
        # Compute energy in window
        energy = 0.0
        for j in range(start, end):
            energy += x[j] * x[j]
        
        power[i] = energy / (end - start)
    
    return power


# Wrapper functions that integrate Numba-optimized code with NumPy operations

def parabolic_interp_opt(yv: np.ndarray, xo: float) -> Tuple[float, float]:
    """
    Optimized parabolic interpolation (wrapper)
    """
    return _compute_parabolic_interp(yv, xo)


def find_peaks_opt(data: np.ndarray, threshold: float) -> np.ndarray:
    """
    Optimized peak finding (wrapper)
    """
    return _find_peaks_above_threshold(data, threshold)


def compute_tracking_costs_opt(candidates: np.ndarray, last_f0: float,
                                reliabilities: np.ndarray) -> np.ndarray:
    """
    Optimized tracking cost computation (wrapper)
    """
    return _compute_tracking_costs(candidates, last_f0, reliabilities)


def zfixpfreq3_opt(fxx: np.ndarray, pif2: np.ndarray, mmp: np.ndarray,
                   dfv: np.ndarray, pm: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Optimized version of zfixpfreq3 - find fixed points in frequency map
    """
    aav = np.abs(pm)
    return _zfixpfreq3_core(fxx, pif2, mmp, dfv, aav)


def compute_acc_indices_opt(fftl: int, wflfso: int, bbf: np.ndarray,
                            nsht: float, ndiv: int, npw: np.ndarray,
                            pc: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Pre-compute indices and scaling for ztestspecspecnormal optimization
    """
    return _compute_acc_indices_and_scaling(fftl, wflfso, bbf, nsht, ndiv, npw, pc)


def apply_window_and_scale_opt(pw: np.ndarray, w2: np.ndarray,
                                all_indices: np.ndarray, scaling_factors: np.ndarray,
                                wflfso: int, ndiv: int) -> np.ndarray:
    """
    Apply windowing and scaling for ztestspecspecnormal optimization
    """
    return _apply_window_and_scale(pw, w2, all_indices, scaling_factors, wflfso, ndiv)


# Performance analysis utilities

def benchmark_optimization():
    """
    Benchmark optimized functions against originals
    """
    import time
    
    print("=== Numba Optimization Benchmarks ===\n")
    
    # Test parabolic interpolation
    yv = np.array([0.5, 1.0, 0.8])
    n_iter = 100000
    
    start = time.time()
    for _ in range(n_iter):
        result = _compute_parabolic_interp(yv, 10.0)
    elapsed = time.time() - start
    print(f"Parabolic interpolation: {n_iter} iterations in {elapsed:.3f}s")
    print(f"  -> {n_iter/elapsed:.0f} ops/sec\n")
    
    # Test peak finding
    data = np.random.randn(1000)
    threshold = 0.5
    n_iter = 10000
    
    start = time.time()
    for _ in range(n_iter):
        peaks = _find_peaks_above_threshold(data, threshold)
    elapsed = time.time() - start
    print(f"Peak finding: {n_iter} iterations in {elapsed:.3f}s")
    print(f"  -> {n_iter/elapsed:.0f} ops/sec\n")
    
    # Test tracking costs
    candidates = np.array([100.0, 150.0, 200.0])
    reliabilities = np.array([0.8, 0.6, 0.4])
    last_f0 = 110.0
    n_iter = 100000
    
    start = time.time()
    for _ in range(n_iter):
        costs = _compute_tracking_costs(candidates, last_f0, reliabilities)
    elapsed = time.time() - start
    print(f"Tracking costs: {n_iter} iterations in {elapsed:.3f}s")
    print(f"  -> {n_iter/elapsed:.0f} ops/sec\n")


if __name__ == "__main__":
    benchmark_optimization()

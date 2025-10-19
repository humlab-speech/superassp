"""
Optimized CPP calculation using vectorization and caching
"""

import numpy as np
from scipy.signal import find_peaks
from scipy.fft import fft, ifft

try:
    from numba import jit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    def jit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator


def get_cpp_optimized(audio: np.ndarray, fs: int, f0_values: np.ndarray,
                      frame_shift: float = 1.0, n_periods: int = 5) -> np.ndarray:
    """
    Optimized CPP calculation with pre-windowing and vectorization
    
    Uses scipy.fft for better performance and window caching
    
    Args:
        audio: Full audio signal
        fs: Sampling rate
        f0_values: F0 array
        frame_shift: Frame shift in milliseconds
        n_periods: Number of pitch periods
    
    Returns:
        cpp: CPP values for each frame
    """
    sample_shift = int(fs * frame_shift / 1000)
    cpp = np.full(len(f0_values), np.nan)
    
    # Pre-compute Hamming windows for common lengths
    window_cache = {}
    
    # Quefrency cutoff
    n_ms = int(fs / 1000)
    
    for k in range(len(f0_values)):
        f0 = f0_values[k]
        
        if np.isnan(f0) or f0 <= 0:
            continue
        
        ks = int(k * sample_shift)
        n0 = fs / f0
        
        start = int(ks - n_periods/2 * n0)
        end = int(ks + n_periods/2 * n0)
        
        if start < 0 or end >= len(audio):
            continue
        
        segment = audio[start:end]
        if len(segment) < 2:
            continue
        
        # Get or create window (caching)
        seg_len = len(segment)
        if seg_len not in window_cache:
            window_cache[seg_len] = np.hamming(seg_len)
        
        windowed = segment * window_cache[seg_len]
        
        # Cepstrum using scipy.fft (faster than numpy.fft)
        spectrum = fft(windowed)
        log_spectrum = np.log(np.abs(spectrum) + 1e-10)
        cepstrum = ifft(log_spectrum)
        
        # Convert to dB
        cepstrum_db = 10 * np.log10(np.abs(cepstrum)**2 + 1e-10)
        
        # Use first half
        half_len = len(cepstrum_db) // 2
        cepstrum_db = cepstrum_db[:half_len]
        
        if n_ms >= half_len:
            continue
        
        search_region = cepstrum_db[n_ms:]
        
        if len(search_region) < 2:
            continue
        
        # Find peaks (vectorized)
        peaks, _ = find_peaks(search_region, distance=int(n0))
        
        if len(peaks) == 0:
            continue
        
        # Find peak closest to expected quefrency
        expected_idx = int(n0)
        if expected_idx >= len(search_region):
            expected_idx = len(search_region) // 2
        
        closest_peak_idx = np.argmin(np.abs(peaks - expected_idx))
        peak_idx = peaks[closest_peak_idx]
        peak_val = search_region[peak_idx]
        
        # Fit linear baseline (vectorized)
        x = np.arange(len(search_region))
        coeffs = np.polyfit(x, search_region, 1)
        baseline_val = np.polyval(coeffs, peak_idx)
        
        cpp[k] = peak_val - baseline_val
    
    return cpp

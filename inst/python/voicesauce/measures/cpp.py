"""
Cepstral Peak Prominence (CPP)
Direct port of func_GetCPP.m (Hillenbrand 1994)
"""

import numpy as np
from scipy.signal import find_peaks


def get_cpp(audio: np.ndarray, fs: int, f0_values: np.ndarray,
           frame_shift: float = 1.0, n_periods: int = 5) -> np.ndarray:
    """
    Calculate Cepstral Peak Prominence
    Port of func_GetCPP.m (Hillenbrand et al., 1994)
    
    Args:
        audio: Full audio signal
        fs: Sampling rate
        f0_values: F0 array
        frame_shift: Frame shift in milliseconds
        n_periods: Number of pitch periods to extract (default 5)
    
    Returns:
        cpp: CPP values for each frame
    """
    sample_shift = int(fs * frame_shift / 1000)
    cpp = np.full(len(f0_values), np.nan)
    
    # Quefrency cutoff: ignore below 1ms
    n_ms = int(fs / 1000)
    
    for k, f0 in enumerate(f0_values):
        if np.isnan(f0) or f0 <= 0:
            continue
        
        # Calculate segment
        ks = int(k * sample_shift)
        n0 = fs / f0
        
        start = int(ks - n_periods/2 * n0)
        end = int(ks + n_periods/2 * n0)
        
        if start < 0 or end >= len(audio):
            continue
        
        # Extract and window
        segment = audio[start:end]
        if len(segment) < 2:
            continue
            
        segment = segment * np.hamming(len(segment))
        
        # Cepstrum: FFT -> log -> IFFT
        spectrum = np.fft.fft(segment)
        log_spectrum = np.log(np.abs(spectrum) + 1e-10)
        cepstrum = np.fft.ifft(log_spectrum)
        
        # Convert to dB
        cepstrum_db = 10 * np.log10(np.abs(cepstrum)**2 + 1e-10)
        
        # Use first half
        half_len = len(cepstrum_db) // 2
        cepstrum_db = cepstrum_db[:half_len]
        
        if n_ms >= len(cepstrum_db):
            continue
        
        # Search for peak in valid region
        search_region = cepstrum_db[n_ms:]
        
        if len(search_region) < 2:
            continue
        
        # Find peaks
        peaks, properties = find_peaks(search_region, distance=int(n0))
        
        if len(peaks) == 0:
            continue
        
        # Find peak closest to expected quefrency (n0)
        expected_idx = int(n0)
        if expected_idx >= len(search_region):
            expected_idx = len(search_region) // 2
            
        closest_peak_idx = np.argmin(np.abs(peaks - expected_idx))
        peak_idx = peaks[closest_peak_idx]
        peak_val = search_region[peak_idx]
        
        # Fit linear baseline
        x = np.arange(len(search_region))
        try:
            coeffs = np.polyfit(x, search_region, 1)
            baseline_val = np.polyval(coeffs, peak_idx)
            
            # CPP = peak height above baseline
            cpp[k] = peak_val - baseline_val
        except:
            continue
    
    return cpp

"""
Harmonic-to-Noise Ratio (HNR)
Direct port of func_GetHNR.m (de Krom 1993)
"""

import numpy as np
from typing import List


def get_hnr(audio: np.ndarray, fs: int, f0_values: np.ndarray,
           frame_shift: float = 1.0, n_periods: int = 5,
           freq_bands: List[float] = None) -> dict:
    """
    Calculate Harmonic-to-Noise Ratio
    Port of func_GetHNR.m (de Krom, 1993)
    
    Args:
        audio: Full audio signal
        fs: Sampling rate
        f0_values: F0 array
        frame_shift: Frame shift in milliseconds
        n_periods: Number of pitch periods to extract (default 5)
        freq_bands: Frequency bands for HNR (default [500, 1500, 2500, 3500])
    
    Returns:
        Dictionary with HNR05, HNR15, HNR25, HNR35
    """
    if freq_bands is None:
        freq_bands = [500, 1500, 2500, 3500]
    
    sample_shift = int(fs * frame_shift / 1000)
    n_bands = len(freq_bands)
    
    hnr_values = {
        'HNR05': np.full(len(f0_values), np.nan),
        'HNR15': np.full(len(f0_values), np.nan),
        'HNR25': np.full(len(f0_values), np.nan),
        'HNR35': np.full(len(f0_values), np.nan)
    }
    
    for k, f0 in enumerate(f0_values):
        if np.isnan(f0) or f0 <= 0:
            continue
        
        # Calculate segment
        ks = int(k * sample_shift)
        n0 = fs / f0
        
        start = int(ks - n_periods/2 * n0)
        end = int(ks + n_periods/2 * n0)
        
        # Make length odd
        if (end - start) % 2 == 0:
            end -= 1
        
        if start < 0 or end >= len(audio):
            continue
        
        # Extract segment
        segment = audio[start:end]
        if len(segment) < 2:
            continue
        
        # Calculate HNR
        try:
            hnr = _calculate_hnr_frame(segment, f0, fs, freq_bands)
            
            # Assign to output
            hnr_keys = ['HNR05', 'HNR15', 'HNR25', 'HNR35']
            for i, key in enumerate(hnr_keys[:len(hnr)]):
                hnr_values[key][k] = hnr[i]
        except:
            continue
    
    return hnr_values


def _calculate_hnr_frame(segment: np.ndarray, f0: float, fs: int, 
                        freq_bands: List[float]) -> np.ndarray:
    """
    Calculate HNR for a single frame
    Based on de Krom (1993) cepstrum-based technique
    
    Args:
        segment: Audio segment
        f0: Fundamental frequency
        fs: Sampling rate
        freq_bands: Frequency bands
    
    Returns:
        HNR values for each frequency band
    """
    n_bins = len(segment)
    n0 = int(fs / f0)
    n0_delta = int(n0 * 0.1)  # Search 10% either side
    
    # Window and FFT
    segment = segment * np.hamming(len(segment))
    spectrum = np.fft.fft(segment, n_bins)
    log_spectrum = np.log10(np.abs(spectrum) + 1e-10)
    cepstrum = np.fft.ifft(log_spectrum)
    
    # Find and lifter out rahmonic peaks
    n_peaks = len(segment) // 2 // n0
    
    for i in range(1, n_peaks + 1):
        center = i * n0
        if center + n0_delta >= len(cepstrum) // 2:
            break
        
        # Find peak in region
        start_idx = max(0, center - n0_delta)
        end_idx = min(len(cepstrum) // 2, center + n0_delta)
        
        region = np.abs(cepstrum[start_idx:end_idx])
        if len(region) == 0:
            continue
            
        peak_idx = start_idx + np.argmax(region)
        
        # Lifter out peak (simple zero-out for now)
        # A more sophisticated approach would use derivative zero-crossings
        lifter_start = max(0, peak_idx - 2)
        lifter_end = min(len(cepstrum), peak_idx + 3)
        cepstrum[lifter_start:lifter_end] = 0
    
    # Reconstruct noise spectrum
    mid_len = len(cepstrum) // 2 + 1
    cepstrum[mid_len:] = np.conj(cepstrum[mid_len-2:0:-1])
    
    noise_spectrum = np.real(np.fft.fft(cepstrum))
    harmonic_spectrum = log_spectrum - noise_spectrum
    
    # Apply baseline corrections
    h_delta = f0 / fs * len(segment)
    for freq in np.arange(h_delta, len(segment)//2, h_delta):
        f_start = int(freq - h_delta)
        f_end = int(freq)
        if f_end >= len(noise_spectrum) // 2:
            break
        region = harmonic_spectrum[f_start:f_end]
        if len(region) > 0:
            b_delta = np.abs(np.min(region))
            noise_spectrum[f_start:f_end] -= b_delta
    
    # Calculate HNR for each frequency band
    harmonic_spectrum = log_spectrum - noise_spectrum
    hnr = np.zeros(len(freq_bands))
    
    for i, freq in enumerate(freq_bands):
        freq_bin = int(freq / fs * len(segment))
        if freq_bin < len(harmonic_spectrum) // 2:
            hnr[i] = 20 * np.mean(harmonic_spectrum[1:freq_bin]) - 20 * np.mean(noise_spectrum[1:freq_bin])
    
    return hnr

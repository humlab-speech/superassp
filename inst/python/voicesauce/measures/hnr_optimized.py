"""
Optimized HNR calculation using scipy.fft and vectorization
Note: numba doesn't support np.fft in nopython mode, so using scipy.fft optimization instead
"""

import numpy as np
from typing import List
from scipy.fft import fft, ifft

NUMBA_AVAILABLE = True  # Set to true since we have scipy optimizations


def _calculate_hnr_frame_fast(segment: np.ndarray, f0: float, fs: int,
                              freq_bands: np.ndarray) -> np.ndarray:
    """
    Fast HNR calculation for a single frame using scipy.fft
    Simplified version that follows the standard implementation

    Args:
        segment: Audio segment
        f0: Fundamental frequency
        fs: Sampling rate
        freq_bands: Frequency bands (as numpy array)

    Returns:
        HNR values for each frequency band
    """
    n_bins = len(segment)
    n0 = int(fs / f0)
    n0_delta = int(n0 * 0.1)

    # Window
    windowed = segment * np.hamming(len(segment))

    # FFT using scipy (faster than numpy)
    spectrum = fft(windowed, n_bins)
    log_spectrum = np.log10(np.abs(spectrum) + 1e-10)
    cepstrum = ifft(log_spectrum)

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

        # Lifter out peak (simple zero-out)
        lifter_start = max(0, peak_idx - 2)
        lifter_end = min(len(cepstrum), peak_idx + 3)
        cepstrum[lifter_start:lifter_end] = 0

    # Reconstruct noise spectrum
    mid_len = len(cepstrum) // 2 + 1
    cepstrum[mid_len:] = np.conj(cepstrum[mid_len-2:0:-1])

    noise_spectrum = np.real(fft(cepstrum))
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
        if freq_bin < len(harmonic_spectrum) // 2 and freq_bin > 1:
            hnr[i] = 20 * np.mean(harmonic_spectrum[1:freq_bin]) - 20 * np.mean(noise_spectrum[1:freq_bin])

    return hnr


def get_hnr_fast(audio: np.ndarray, fs: int, f0_values: np.ndarray,
                sample_shift: int, n_periods: int,
                freq_bands: np.ndarray) -> np.ndarray:
    """
    Fast HNR calculation using scipy.fft optimization

    Args:
        audio: Full audio signal
        fs: Sampling rate
        f0_values: F0 array
        sample_shift: Samples per frame shift
        n_periods: Number of pitch periods
        freq_bands: Frequency bands as numpy array

    Returns:
        HNR values array (n_frames x n_bands)
    """
    n_frames = len(f0_values)
    n_bands = len(freq_bands)
    hnr_array = np.full((n_frames, n_bands), np.nan)

    for k in range(n_frames):
        f0 = f0_values[k]

        if np.isnan(f0) or f0 <= 0:
            continue

        ks = k * sample_shift
        n0 = fs / f0

        start = int(ks - n_periods/2 * n0)
        end = int(ks + n_periods/2 * n0)

        # Make length odd
        if (end - start) % 2 == 0:
            end -= 1

        if start < 0 or end >= len(audio):
            continue

        segment = audio[start:end]
        if len(segment) < 2:
            continue

        # Calculate HNR
        hnr_frame = _calculate_hnr_frame_fast(segment, f0, fs, freq_bands)
        hnr_array[k, :] = hnr_frame

    return hnr_array


def get_hnr_optimized(audio: np.ndarray, fs: int, f0_values: np.ndarray,
                     frame_shift: float = 1.0, n_periods: int = 5,
                     freq_bands: List[float] = None) -> dict:
    """
    Optimized HNR calculation wrapper

    Args:
        audio: Full audio signal
        fs: Sampling rate
        f0_values: F0 array
        frame_shift: Frame shift in milliseconds
        n_periods: Number of pitch periods
        freq_bands: Frequency bands for HNR

    Returns:
        Dictionary with HNR05, HNR15, HNR25, HNR35
    """
    if freq_bands is None:
        freq_bands = [500, 1500, 2500, 3500]

    sample_shift = int(fs * frame_shift / 1000)

    hnr_values = {
        'HNR05': np.full(len(f0_values), np.nan),
        'HNR15': np.full(len(f0_values), np.nan),
        'HNR25': np.full(len(f0_values), np.nan),
        'HNR35': np.full(len(f0_values), np.nan)
    }

    if NUMBA_AVAILABLE:
        freq_bands_arr = np.array(freq_bands, dtype=np.float64)
        hnr_array = get_hnr_fast(audio, fs, f0_values, sample_shift,
                                 n_periods, freq_bands_arr)

        # Assign to output dictionary
        hnr_keys = ['HNR05', 'HNR15', 'HNR25', 'HNR35']
        for i, key in enumerate(hnr_keys[:hnr_array.shape[1]]):
            hnr_values[key] = hnr_array[:, i]
    else:
        # Fallback to standard implementation
        from .hnr import get_hnr
        return get_hnr(audio, fs, f0_values, frame_shift, n_periods, freq_bands)

    return hnr_values

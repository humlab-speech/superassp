"""
PeakSlope measurement using wavelet decomposition.

References
----------
Kane, J., Gobl, C., (2011) "Identifying regions of non-modal phonation using
features of the wavelet transform", Proceedings of Interspeech
"""

import numpy as np
import pywt


def get_peak_slope(s, fs):
    """
    Calculate peakSlope measurement on a fixed frame basis.

    The peakSlope measures the slope of peak amplitudes across wavelet scales,
    useful for identifying non-modal phonation regions.

    Parameters
    ----------
    s : array_like
        Speech signal (samples)
    fs : float
        Sampling frequency (Hz)

    Returns
    -------
    peak_slope : ndarray
        PeakSlope parameter measured every 10 ms

    Notes
    -----
    This is a simplified Python implementation. The original MATLAB version uses
    custom Daless wavelets. This version uses PyWavelets with Daubechies wavelets
    as an approximation.

    The MATLAB version uses daless_MW mother wavelets which are not available
    in standard Python libraries. For faithful reproduction, the custom wavelet
    functions would need to be ported.
    """
    s = np.asarray(s).flatten()

    # Frame parameters
    frame_len_ms = 40  # Frame length to ensure one pulse at f0=25 Hz
    frame_shift_ms = 10  # Frame shift 10 ms
    frame_len = int((frame_len_ms / 1000) * fs)
    frame_shift = int((frame_shift_ms / 1000) * fs)

    n_frames = int(np.floor((len(s) - frame_len) / frame_shift))
    peak_slope = np.zeros(n_frames)

    # Wavelet decomposition using PyWavelets
    # Use multi-level DWT as approximation to custom Daless wavelets
    max_level = min(6, pywt.dwt_max_level(len(s), 'db4'))
    coeffs = pywt.wavedec(s, 'db4', level=max_level)

    # Reconstruct each level
    y = []
    for i in range(len(coeffs)):
        c_temp = [np.zeros_like(c) for c in coeffs]
        c_temp[i] = coeffs[i]
        y_i = pywt.waverec(c_temp, 'db4')
        # Ensure same length as original signal
        if len(y_i) > len(s):
            y_i = y_i[:len(s)]
        elif len(y_i) < len(s):
            y_i = np.pad(y_i, (0, len(s) - len(y_i)))
        y.append(y_i)

    # Measure peakSlope per frame
    start = 0
    finish = start + frame_len

    for m in range(n_frames):
        if finish > len(s):
            break

        maxima = []
        for n in range(len(y)):
            maxima.append(np.max(np.abs(y[n][start:finish])))

        # Reverse order and convert to log scale
        maxima = np.log10(np.array(maxima[::-1]) + 1e-10)
        t = np.arange(len(maxima))

        # Linear regression
        if len(t) > 1:
            p = np.polyfit(t, maxima, 1)
            peak_slope[m] = p[0]

        start += frame_shift
        finish = start + frame_len

    # Remove NaN values
    peak_slope[np.isnan(peak_slope)] = 0

    print("Warning: This is a simplified implementation using standard wavelets.")
    print("For faithful reproduction, Daless wavelets from MATLAB code should be ported.")

    return peak_slope

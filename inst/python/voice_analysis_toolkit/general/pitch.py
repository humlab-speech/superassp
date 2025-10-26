"""
Pitch tracking functions.

This module contains Python implementations of pitch tracking algorithms from
the Voice Analysis Toolkit, particularly the SRH (Summation of Residual Harmonics)
method developed by Thomas Drugman.

References
----------
T. Drugman, A. Alwan, "Joint Robust Voicing Detection and Pitch Estimation
Based on Residual Harmonics", Interspeech11, Firenze, Italy, 2011
"""

import numpy as np
from scipy.signal import resample, windows
from .signal_utils import get_lpc_residual


def srh_estimate_pitch(sig, fs, f0_min, f0_max):
    """
    Estimate pitch using SRH (Summation of Residual Harmonics) criterion.

    Parameters
    ----------
    sig : array_like
        Input signal (LPC residual)
    fs : float
        Sampling frequency (Hz)
    f0_min : float
        Minimum F0 value for search (Hz)
    f0_max : float
        Maximum F0 value for search (Hz)

    Returns
    -------
    f0 : ndarray
        F0 values with 10ms hop size, padded with 5 zeros at start and end
    srh_val : ndarray
        SRH values with 10ms hop size, padded with 5 zeros at start and end
    time : ndarray
        Time indices (sample positions) for each frame
    """
    sig = np.asarray(sig).flatten()

    start = 0
    stop = int(np.round(100 / 1000 * fs))  # 100ms window
    shift = int(np.round(10 / 1000 * fs))   # 10ms shift

    n_frames = int(np.floor((len(sig) - stop) / shift)) + 1
    srh_tot = np.zeros((n_frames, f0_max))
    f0s = np.zeros(n_frames)
    srh_vals = np.zeros(n_frames)
    time = np.zeros(n_frames)

    # Use symmetric Blackman window to match MATLAB
    black_win = windows.blackman(stop - start, sym=True)
    index = 0

    while stop <= len(sig):
        time[index] = start
        seg = sig[start:stop]
        seg = seg * black_win
        seg = seg - np.mean(seg)

        # FFT
        spec = np.fft.fft(seg, int(fs))
        spec = np.abs(spec[:int(fs // 2)])
        spec = spec / (np.sqrt(np.sum(spec**2)) + 1e-10)

        srhs = np.zeros(f0_max)

        # SRH spectral criterion
        for freq in range(f0_min, f0_max):
            # Sum harmonics at integer multiples
            harmonics_sum = (spec[freq] + spec[2*freq] + spec[3*freq] +
                           spec[4*freq] + spec[5*freq])
            # Subtract inter-harmonics at half-integer multiples
            interharmonics_sum = (spec[int(np.round(1.5*freq))] +
                                spec[int(np.round(2.5*freq))] +
                                spec[int(np.round(3.5*freq))] +
                                spec[int(np.round(4.5*freq))])

            srhs[freq] = harmonics_sum - interharmonics_sum

        srh_tot[index, :] = srhs

        # Find maximum
        max_idx = np.argmax(srhs)
        f0_frame = max_idx

        f0s[index] = f0_frame
        srh_vals[index] = srhs[f0_frame]

        start += shift
        stop += shift
        index += 1

    # Pad with zeros (5 frames at start and end)
    f0 = np.concatenate([np.zeros(5), f0s, np.zeros(5)])
    srh_val = np.concatenate([np.zeros(5), srh_vals, np.zeros(5)])

    return f0, srh_val, time


def srh_pitch_tracking(wave, fs, f0_min=20, f0_max=500):
    """
    Track pitch using SRH (Summation of Residual Harmonics) method.

    This function performs robust F0 and voicing detection based on the
    residual harmonics of the LPC residual signal.

    Parameters
    ----------
    wave : array_like
        Speech signal
    fs : float
        Sampling frequency (Hz)
    f0_min : float, optional
        Minimum F0 value for search (Hz), default: 20
    f0_max : float, optional
        Maximum F0 value for search (Hz), default: 500

    Returns
    -------
    f0 : ndarray
        Fundamental frequency contour with 10ms hop size (Hz)
    vuv_decisions : ndarray
        Voiced/unvoiced decisions (0=unvoiced, 1=voiced)
    srh_val : ndarray
        SRH values (confidence measure)
    time : ndarray
        Time indices for each frame

    References
    ----------
    T. Drugman, A. Alwan, "Joint Robust Voicing Detection and Pitch Estimation
    Based on Residual Harmonics", Interspeech11, Firenze, Italy, 2011

    Notes
    -----
    Original MATLAB code by Thomas Drugman, TCTS Lab, University of Mons, Belgium.
    Python conversion: 2025
    """
    wave = np.asarray(wave).flatten()

    # Resample to 16 kHz if needed
    if fs > 16000:
        wave = resample(wave, int(len(wave) * 16000 / fs))
        fs = 16000

    lpc_order = int(np.round(3 / 4 * fs / 1000))
    n_iter = 2

    # Get LPC residual
    L = int(np.round(25 / 1000 * fs))      # 25ms window
    shift = int(np.round(5 / 1000 * fs))   # 5ms shift
    res, _ = get_lpc_residual(wave, L, shift, lpc_order)

    # Estimate pitch track in 2 iterations
    for iteration in range(n_iter):
        f0, srh_val, time = srh_estimate_pitch(res, fs, f0_min, f0_max)

        # Find valid F0 estimates
        posi_tmp = srh_val > 0.1

        if np.sum(posi_tmp) > 1:
            f0_mean_est = np.median(f0[posi_tmp])

            # Update F0 search range
            f0_min = int(np.round(0.5 * f0_mean_est))
            f0_max = int(np.round(2 * f0_mean_est))

    # Voiced-unvoiced decisions based on SRH value
    vuv_decisions = np.zeros(len(f0), dtype=int)
    threshold = 0.07
    if np.std(srh_val) > 0.05:
        threshold = 0.085

    pos = srh_val > threshold
    vuv_decisions[pos] = 1

    return f0, vuv_decisions, srh_val, time

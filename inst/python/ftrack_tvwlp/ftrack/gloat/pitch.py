"""
SRH (Summation of Residual Harmonics) pitch tracking.

Based on:
T.Drugman, A.Alwan, "Joint Robust Voicing Detection and Pitch Estimation
Based on Residual Harmonics", Interspeech11, Firenze, Italy, 2011
"""

import numpy as np
from scipy import signal
from .utils import get_lpc_residual


def srh_pitch_tracking(wave, fs, f0_min, f0_max):
    """
    SRH pitch tracking algorithm.

    Parameters
    ----------
    wave : ndarray
        Speech signal
    fs : int
        Sampling frequency in Hz
    f0_min : float
        Minimum F0 to search (Hz)
    f0_max : float
        Maximum F0 to search (Hz)

    Returns
    -------
    f0 : ndarray
        F0 values with 10ms hop size
    vuv_decisions : ndarray
        Voiced/unvoiced decisions (1=voiced, 0=unvoiced)
    srh_val : ndarray
        SRH values for each frame
    """
    wave = np.asarray(wave, dtype=np.float64)

    # Resample to 16kHz if necessary
    if fs > 16000:
        n_samples_new = int(len(wave) * 16000 / fs)
        wave = signal.resample(wave, n_samples_new)
        fs = 16000

    # LPC order
    lpc_order = int(np.round(3 / 4 * fs / 1000))
    n_iter = 2

    # Get LPC residual
    res = get_lpc_residual(wave, int(np.round(25 / 1000 * fs)),
                           int(np.round(5 / 1000 * fs)), lpc_order)

    # Estimate pitch in 2 iterations
    for iteration in range(n_iter):
        f0, srh_val = _srh_estimate_pitch(res, fs, f0_min, f0_max)

        # Refine F0 range based on median estimate
        pos_tmp = srh_val > 0.1
        if np.sum(pos_tmp) > 1:
            f0_mean_est = np.median(f0[pos_tmp])
            f0_min = int(np.round(0.5 * f0_mean_est))
            f0_max = int(np.round(2 * f0_mean_est))

    # Voiced-unvoiced decisions
    vuv_decisions = np.zeros(len(f0), dtype=int)
    threshold = 0.07
    if np.std(srh_val) > 0.05:
        threshold = 0.085

    vuv_decisions[srh_val > threshold] = 1

    return f0, vuv_decisions, srh_val


def _srh_estimate_pitch(sig, fs, f0_min, f0_max):
    """
    Internal function to estimate pitch using SRH criterion.

    Parameters
    ----------
    sig : ndarray
        Signal (typically LPC residual)
    fs : int
        Sampling frequency
    f0_min : int
        Minimum F0 (Hz)
    f0_max : int
        Maximum F0 (Hz)

    Returns
    -------
    f0s : ndarray
        F0 estimates for each frame
    srh_val : ndarray
        SRH values for each frame
    """
    start = 0
    stop = int(np.round(100 / 1000 * fs))
    shift = int(np.round(10 / 1000 * fs))

    n_frames = int(np.floor((len(sig) - stop) / shift)) + 1
    srh_tot = np.zeros((n_frames, f0_max))
    f0s = np.zeros(n_frames)
    srh_val = np.zeros(n_frames)

    black_win = signal.windows.blackman(stop - start)

    index = 0
    while stop <= len(sig):
        seg = sig[start:stop].copy()
        seg = seg * black_win
        seg = seg - np.mean(seg)

        # FFT
        spec = np.fft.fft(seg, n=fs)
        spec = np.abs(spec[:fs // 2])
        spec_energy = np.sqrt(np.sum(spec ** 2))
        if spec_energy > 0:
            spec = spec / spec_energy

        # SRH spectral criterion
        srhs = np.zeros(f0_max)

        for freq in range(f0_min, f0_max):
            # Harmonic peaks
            h1 = spec[freq] if freq < len(spec) else 0
            h2 = spec[2 * freq] if 2 * freq < len(spec) else 0
            h3 = spec[3 * freq] if 3 * freq < len(spec) else 0
            h4 = spec[4 * freq] if 4 * freq < len(spec) else 0
            h5 = spec[5 * freq] if 5 * freq < len(spec) else 0

            # Inter-harmonic valleys
            v1 = spec[int(np.round(1.5 * freq))] if int(np.round(1.5 * freq)) < len(spec) else 0
            v2 = spec[int(np.round(2.5 * freq))] if int(np.round(2.5 * freq)) < len(spec) else 0
            v3 = spec[int(np.round(3.5 * freq))] if int(np.round(3.5 * freq)) < len(spec) else 0
            v4 = spec[int(np.round(4.5 * freq))] if int(np.round(4.5 * freq)) < len(spec) else 0

            srhs[freq] = (h1 + h2 + h3 + h4 + h5) - (v1 + v2 + v3 + v4)

        srh_tot[index, :] = srhs

        # Find maximum
        max_idx = np.argmax(srhs)
        f0_frame = max_idx

        f0s[index] = f0_frame
        srh_val[index] = srhs[f0_frame]

        start += shift
        stop += shift
        index += 1

    # Pad with zeros (MATLAB adds 5 zeros at beginning and end)
    f0 = np.concatenate([np.zeros(5), f0s, np.zeros(5)])
    srh_val = np.concatenate([np.zeros(5), srh_val, np.zeros(5)])

    return f0, srh_val

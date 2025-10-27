"""
Signal processing utilities for formant tracking.
"""

import numpy as np
from scipy import signal


def lpc(x, order):
    """
    Linear Predictive Coding using Levinson-Durbin recursion.

    Parameters
    ----------
    x : ndarray
        Input signal
    order : int
        LPC order

    Returns
    -------
    a : ndarray
        LPC coefficients [1, a1, a2, ..., a_order]
    e : float
        Prediction error
    """
    x = np.asarray(x, dtype=np.float64)

    # Compute autocorrelation
    r = np.correlate(x, x, mode='full')
    r = r[len(r) // 2:]  # Keep only non-negative lags
    r = r[:order + 1]

    # Levinson-Durbin recursion
    # Initialize
    a = np.zeros(order + 1, dtype=np.float64)
    a[0] = 1.0
    e = r[0]

    if e == 0:
        # Avoid division by zero
        return a, 0.0

    for i in range(1, order + 1):
        # Compute reflection coefficient
        lambda_i = r[i]
        for j in range(1, i):
            lambda_i -= a[j] * r[i - j]
        lambda_i /= e

        # Update coefficients
        a_new = a.copy()
        a_new[i] = lambda_i
        for j in range(1, i):
            a_new[j] = a[j] - lambda_i * a[i - j]

        a = a_new

        # Update error
        e = e * (1 - lambda_i ** 2)

        if e <= 0:
            break

    return a, e


def get_lpc_residual(wave, L, shift, order):
    """
    Compute LPC residual signal using overlapping windows.

    Parameters
    ----------
    wave : ndarray
        Input speech signal
    L : int
        Window length in samples (e.g., 25ms)
    shift : int
        Window shift in samples (e.g., 5ms)
    order : int
        LPC order

    Returns
    -------
    res : ndarray
        LPC residual signal (same length as wave)
    """
    wave = np.asarray(wave, dtype=np.float64)
    n_samples = len(wave)

    # Hanning window
    hann_win = signal.windows.hann(L + 1)

    res = np.zeros(n_samples, dtype=np.float64)

    start = 0
    stop = start + L

    while stop < n_samples:
        segment = wave[start:stop + 1].copy()
        segment = segment * hann_win

        # Compute LPC coefficients
        A, e = lpc(segment, order)

        # Compute inverse filter output (residual)
        inv = signal.lfilter(A, 1, segment)

        # Normalize residual energy to match segment energy
        seg_energy = np.sum(segment ** 2)
        inv_energy = np.sum(inv ** 2)
        if inv_energy > 0:
            inv = inv * np.sqrt(seg_energy / inv_energy)

        # Accumulate residual (overlap-add)
        res[start:stop + 1] += inv

        start += shift
        stop += shift

    # Normalize
    res_max = np.max(np.abs(res))
    if res_max > 0:
        res = res / res_max

    return res

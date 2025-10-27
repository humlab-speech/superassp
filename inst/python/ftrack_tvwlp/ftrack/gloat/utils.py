"""
Utility functions for GLOAT toolkit.
"""

import numpy as np
from scipy import signal as sp_signal


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
    r = r[len(r) // 2:]
    r = r[:order + 1]

    # Levinson-Durbin recursion
    a = np.zeros(order + 1, dtype=np.float64)
    a[0] = 1.0
    e = r[0]

    if e == 0:
        return a, 0.0

    for i in range(1, order + 1):
        lambda_i = r[i]
        for j in range(1, i):
            lambda_i -= a[j] * r[i - j]
        lambda_i /= e

        a_new = a.copy()
        a_new[i] = lambda_i
        for j in range(1, i):
            a_new[j] = a[j] - lambda_i * a[i - j]

        a = a_new
        e = e * (1 - lambda_i ** 2)

        if e <= 0:
            break

    return a, e


def get_lpc_residual(wave, L, shift, order):
    """
    Compute LPC residual signal.

    Parameters
    ----------
    wave : ndarray
        Input speech signal
    L : int
        Window length in samples
    shift : int
        Window shift in samples
    order : int
        LPC order

    Returns
    -------
    res : ndarray
        LPC residual signal
    """
    wave = np.asarray(wave, dtype=np.float64)
    n_samples = len(wave)

    hann_win = sp_signal.windows.hann(L + 1)
    res = np.zeros(n_samples, dtype=np.float64)

    start = 0
    stop = start + L

    while stop < n_samples:
        segment = wave[start:stop + 1].copy()
        segment = segment * hann_win

        A, e = lpc(segment, order)
        inv = sp_signal.lfilter(A, 1, segment)

        seg_energy = np.sum(segment ** 2)
        inv_energy = np.sum(inv ** 2)
        if inv_energy > 0:
            inv = inv * np.sqrt(seg_energy / inv_energy)

        res[start:stop + 1] += inv

        start += shift
        stop += shift

    res_max = np.max(np.abs(res))
    if res_max > 0:
        res = res / res_max

    return res

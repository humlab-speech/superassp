"""Resonator for creaky voice detection (RCVD) applied to residual signal."""

import numpy as np
from scipy.signal import filtfilt


def rcvd_reson_gci(res, fs, F0mean):
    """
    Apply resonator used in RCVD (creaky voice detection) to LP-residual.

    Parameters
    ----------
    res : array_like
        LP-residual signal
    fs : float
        Sampling frequency (Hz)
    F0mean : float
        Mean fundamental frequency (Hz)

    Returns
    -------
    y : ndarray
        Resonator output, normalized to [-1, 1]

    Notes
    -----
    The resonator is configured with a narrow bandwidth using settings
    from the RCVD creaky voice detection algorithm.
    """
    res = np.asarray(res).flatten()

    # Configure resonator
    phi = 2 * np.pi * 1 * F0mean / fs
    rho = 0.9  # Narrow bandwidth

    # Resonator coefficients
    b = np.array([1, 0, 0])
    a = np.array([1, -2 * rho * np.cos(phi), rho**2])

    # Apply zero-phase filtering (forward and backward)
    rep = filtfilt(b, a, res)

    # Normalize
    max_val = np.max(np.abs(rep))
    if max_val > 0:
        y = rep / max_val
    else:
        y = rep

    return y

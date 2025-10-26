"""
IAIF (Iterative Adaptive Inverse Filtering) implementation.

This module implements the IAIF algorithm for glottal source estimation.

References
----------
P. Alku, "Glottal wave analysis with Pitch Synchronous Iterative Adaptive
Inverse Filtering", Speech Communication, vol. 11, no. 2-3, pp. 109-118, 1992.
"""

import numpy as np
from scipy.signal import firwin, lfilter, hamming
from .signal_utils import calc_residual


def iaif(x, fs, gci=None, p=None):
    """
    Perform Iterative and Adaptive Inverse Filtering (IAIF).

    The IAIF algorithm estimates the glottal flow derivative by iteratively
    separating the vocal tract and glottal source contributions.

    Parameters
    ----------
    x : array_like
        Speech signal (in samples)
    fs : float
        Sampling frequency (Hz)
    gci : array_like, optional
        Glottal closure instants (in samples). If not provided, SE_VQ will be used.
    p : int, optional
        LPC prediction order. If not provided, defaults to round(fs/1000) + 2

    Returns
    -------
    g_iaif : ndarray
        Glottal flow derivative estimate
    ar_lpc : ndarray
        LPC coefficients from final estimation
    e_lpc : ndarray
        Prediction errors from final estimation

    Notes
    -----
    Original MATLAB implementation by John Kane @ The Phonetics and Speech Lab,
    Trinity College Dublin, August 2012.
    Python conversion: 2025

    The algorithm follows these steps:
    1. High-pass filter to remove DC component
    2. Pre-emphasis with LPC order 1
    3. First glottal source estimation with LPC order p
    4. Vocal tract estimation with LPC order 4
    5. Second (refined) glottal source estimation with LPC order p
    """
    x = np.asarray(x).flatten()

    # Set default LPC filter order
    if p is None:
        p = int(np.round(fs / 1000)) + 2

    # Get GCI if not provided
    if gci is None:
        from ..se_vq import se_vq
        gci, _, _, _ = se_vq(x, fs)

    gci = np.asarray(gci, dtype=int)

    # ------------------------------------------------
    # Part 1: High-pass filter to eliminate DC component
    nc = 704
    fc = 50  # 30 Hz in original Alku (1992), 50 Hz in implementation
    nfc = fc / (fs / 2)

    # Design FIR high-pass filter
    b_hp = firwin(nc + 1, nfc, pass_zero=False)

    # Apply filter
    x_hp = lfilter(b_hp, 1, x)

    # Compensate for filter delay
    x_filt = np.concatenate([x_hp[nc//2:], x[-(nc//2):]])

    # ------------------------------------------------
    # Parts 2 & 3: Pre-emphasis (LPC order 1)
    ord_lpc1 = 1
    x_emph, ar_lpc1, e_lpc1 = calc_residual(x_filt, x_filt, ord_lpc1, gci)

    # ------------------------------------------------
    # Parts 4 & 5: First estimation of glottal source derivative
    ord_lpc2 = p
    residual1, ar_lpc2, e_lpc2 = calc_residual(x_filt, x_emph, ord_lpc2, gci)

    # Part 6: Integration (cancelling lip radiation)
    ug1 = residual1

    # ------------------------------------------------
    # Parts 7 & 8: Elimination of source effect from speech spectrum
    ord_lpc3 = 4
    vt_signal, ar_lpc3, e_lpc3 = calc_residual(x_filt, ug1, ord_lpc3, gci)

    # ------------------------------------------------
    # Parts 9 & 10: Second estimation of glottal source signal
    ord_lpc4 = p
    residual2, ar_lpc, e_lpc = calc_residual(x_filt, vt_signal, ord_lpc4, gci)

    # Final glottal flow derivative estimate
    g_iaif = residual2

    return g_iaif, ar_lpc, e_lpc

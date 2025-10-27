"""
Formant extraction from time-varying LP coefficients.
"""

import numpy as np


def tvlptoformants_akitofi(aki, Nx, npeaks, fs):
    """
    Extract formant frequencies from time-varying LP coefficients.

    Parameters
    ----------
    aki : ndarray
        Time-varying LP coefficient matrix of shape [q+1, p]
        where q is polynomial order and p is LP order
    Nx : int
        Number of time samples
    npeaks : int
        Number of formants to extract
    fs : float
        Sampling frequency in Hz

    Returns
    -------
    fi : ndarray
        Formant frequencies of shape [Nx, npeaks] in Hz
    ak : ndarray
        Time-varying LPC filter coefficients of shape [p+1, Nx]
    """
    q_plus_1, p = aki.shape
    q = q_plus_1 - 1

    tn = np.arange(Nx, dtype=np.float64)

    # Reconstruct time-varying LP coefficients
    # a_k(n) = sum_{i=0}^{q} aki[i,k] * n^i
    akn = np.zeros((p, Nx), dtype=np.float64)

    for k in range(p):
        for i in range(q + 1):
            akn[k, :] += aki[i, k] * (tn ** i)

    # Form LPC filter coefficients: [1, -a1(n), -a2(n), ..., -ap(n)]
    ak = np.vstack([np.ones((1, Nx)), -akn])

    # Extract formants at each time instant
    fi = np.zeros((Nx, npeaks), dtype=np.float64)

    for i in range(Nx):
        # Find roots of the LPC polynomial at time i
        roots_i = np.roots(ak[:, i])

        # Convert to frequencies: angle(root) * fs / (2*pi)
        # Only keep positive frequencies
        angles = np.angle(roots_i)
        freqs = angles * fs / (2 * np.pi)
        freqs_pos = freqs[freqs > 0]

        # Sort and take the first npeaks
        freqs_sorted = np.sort(freqs_pos)

        # Fill in formants (up to available or npeaks, whichever is smaller)
        n_available = min(npeaks, len(freqs_sorted))
        if n_available > 0:
            fi[i, :n_available] = freqs_sorted[:n_available]

    return fi, ak

"""
LPC (Linear Predictive Coding) utility functions.

This module contains Python implementations of LPC-related functions from the
Voice Analysis Toolkit, originally implemented in MATLAB by Mike Brookes
(VoiceBox) and used by John Kane.

References:
    - VoiceBox MATLAB toolbox for speech processing
    - http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
"""

import numpy as np
from scipy import signal
from scipy.linalg import toeplitz


def lpcauto(s, p=12, t=None):
    """
    Perform autocorrelation LPC analysis.

    Parameters
    ----------
    s : array_like
        Input signal
    p : int, optional
        LPC order (default: 12)
    t : array_like, optional
        Frame size details: [len, anal, skip] where:
        - len: length of the frame
        - anal: analysis length (default: len)
        - skip: samples to skip at beginning (default: 0)

    Returns
    -------
    ar : ndarray
        AR coefficients with ar[:,0] = 1, shape (nf, p+1)
    e : ndarray
        Energy in the residual, shape (nf,)
    k : ndarray
        First and last sample of analysis interval, shape (nf, 2)

    Notes
    -----
    This is a Python implementation of lpcauto.m from VoiceBox.
    The analysis interval is always multiplied by a Hamming window.
    """
    s = np.asarray(s).flatten()

    if t is None:
        t = [[len(s), len(s), 0]]
    else:
        t = np.atleast_2d(t)

    nf, ng = t.shape
    if ng < 2:
        t = np.column_stack([t, t])
    if ng < 3:
        t = np.column_stack([t, np.zeros((nf, 1))])

    if nf == 1:
        nf = int(np.floor(1 + (len(s) - t[0, 1] - t[0, 2]) / t[0, 0]))
        tr = 0
    else:
        tr = 1

    ar = np.zeros((nf, p + 1))
    ar[:, 0] = 1
    e = np.zeros(nf)
    k = np.zeros((nf, 2), dtype=int)

    t1 = 1
    it = 0
    nw = -1
    r = np.arange(p + 1)

    for jf in range(nf):
        k[jf, 0] = int(np.ceil(t1 + t[it, 2]))
        k[jf, 1] = int(np.ceil(t1 + t[it, 2] + t[it, 1] - 1))
        cs = np.arange(k[jf, 0] - 1, k[jf, 1])  # -1 for 0-indexing
        nc = len(cs)
        pp = min(p, nc)
        dd = s[cs]

        if nc != nw:
            ww = np.hamming(nc)
            nw = nc
            y = np.zeros(nc + p)
            c = np.arange(nc)

        wd = dd * ww  # windowed data
        y[:nc] = wd  # data with p appended zeros

        # Create data matrix
        z = np.zeros((nc, pp + 1))
        for i in range(pp + 1):
            z[:, i] = y[c + r[i]]

        rr = wd @ z
        rm = toeplitz(rr[:pp])
        rk = np.linalg.matrix_rank(rm)

        if rk > 0:
            if rk < pp:
                rm = rm[:rk, :rk]
                ar[jf, 1:rk+1] = -np.linalg.solve(rm, rr[1:rk+1])
            else:
                ar[jf, 1:pp+1] = -np.linalg.solve(rm, rr[1:pp+1])

        e[jf] = rr @ ar[jf, :pp+1]

        t1 += t[it, 0]
        it += tr

    return ar, e, k


def lpcar2rf(ar):
    """
    Convert autoregressive coefficients to reflection coefficients.

    Parameters
    ----------
    ar : ndarray
        Autoregressive coefficients, shape (nf, p+1)

    Returns
    -------
    rf : ndarray
        Reflection coefficients with rf[:,0]=1, shape (nf, p+1)
    """
    ar = np.atleast_2d(ar)
    nf, p1 = ar.shape

    if p1 == 1:
        return np.ones((nf, 1))

    # Normalize if needed
    if np.any(ar[:, 0] != 1):
        ar = ar / ar[:, [0]]

    rf = ar.copy()

    for j in range(p1 - 2, 0, -1):
        k = rf[:, j + 1]
        d = 1 / (1 - k**2)
        for i in range(1, j + 1):
            rf[:, i] = (rf[:, i] - k * rf[:, j + 1 - i]) * d

    return rf


def lpcar2ra(ar):
    """
    Convert AR filter to inverse filter autocorrelation coefficients.

    Parameters
    ----------
    ar : ndarray
        AR coefficients, shape (nf, p+1)

    Returns
    -------
    ra : ndarray
        Inverse filter autocorrelation coefficients, shape (nf, p+1)
    """
    ar = np.atleast_2d(ar)
    nf, p1 = ar.shape
    ra = np.zeros((nf, p1))

    for i in range(p1):
        ra[:, i] = np.sum(ar[:, :p1-i] * ar[:, i:p1], axis=1)

    return ra


def lpcrf2rr(rf, p=None):
    """
    Convert reflection coefficients to autocorrelation coefficients.

    Parameters
    ----------
    rf : ndarray
        Reflection coefficients, shape (nf, n+1)
    p : int, optional
        Number of rr coefficients to calculate (default=n)

    Returns
    -------
    rr : ndarray
        Autocorrelation coefficients, shape (nf, p+1)
    ar : ndarray
        AR filter coefficients, shape (nf, n+1)
    """
    rf = np.atleast_2d(rf)
    nf, p1 = rf.shape
    p0 = p1 - 1

    if p0 == 0:
        return np.ones((nf, 1)), np.ones((nf, 1))

    a = rf[:, 1:2]
    rr = np.column_stack([np.ones(nf), -a[:, 0], np.zeros((nf, p0 - 1))])
    e = (a[:, 0]**2 - 1)

    for n in range(1, p0):
        k = rf[:, n + 1]
        rr[:, n + 1] = k * e - np.sum(rr[:, n:0:-1] * a, axis=1)
        a_new = np.zeros((nf, n + 1))
        a_new[:, :n] = a + k[:, np.newaxis] * a[:, n-1::-1]
        a_new[:, n] = k
        a = a_new
        e = e * (1 - k**2)

    ar = np.column_stack([np.ones(nf), a])
    r0 = 1 / np.sum(rr * ar, axis=1)
    rr = rr * r0[:, np.newaxis]

    if p is not None:
        if p < p0:
            rr = rr[:, :p+1]
        else:
            extra = np.zeros((nf, p - p0))
            rr = np.column_stack([rr, extra])
            af = -ar[:, p1-1:0:-1]
            for i in range(p0, p):
                rr[:, i + 1] = np.sum(af * rr[:, i-p0:i+1], axis=1)

    return rr, ar


def lpcar2rr(ar, p=None):
    """
    Convert autoregressive coefficients to autocorrelation coefficients.

    Parameters
    ----------
    ar : ndarray
        Autoregressive coefficients including 0th coefficient, shape (nf, n+1)
    p : int, optional
        Number of autocorrelation coefficients (default: n)

    Returns
    -------
    rr : ndarray
        Autocorrelation coefficients including 0th order, shape (nf, p+1)
    """
    ar = np.atleast_2d(ar)
    k = ar[:, 0]**(-2)

    if ar.shape[1] == 1:
        return k.reshape(-1, 1)

    rf = lpcar2rf(ar)

    if p is not None:
        rr, _ = lpcrf2rr(rf, p)
        return rr * k[:, np.newaxis]
    else:
        rr, _ = lpcrf2rr(rf)
        return rr * k[:, np.newaxis]


def distitar(ar1, ar2, mode='d'):
    """
    Calculate the Itakura distance between AR coefficients.

    Parameters
    ----------
    ar1, ar2 : ndarray
        AR coefficient sets to compare. Each row contains a set of coefficients.
    mode : str, optional
        'x': Full distance matrix from every row of ar1 to every row of ar2
        'd': Distance between corresponding rows (default)
        'e': Calculate exp(d) instead of d

    Returns
    -------
    d : ndarray
        Itakura distances

    Notes
    -----
    The Itakura distance is gain-independent and measures spectral dissimilarity.

    References
    ----------
    A.H. Gray Jr and J.D. Markel, "Distance measures for speech processing",
    IEEE ASSP-24(5): 380-391, Oct 1976
    """
    ar1 = np.atleast_2d(ar1)
    ar2 = np.atleast_2d(ar2)

    nf1 = ar1.shape[0]
    nf2 = ar2.shape[0]

    m2 = lpcar2ra(ar2)
    m2[:, 0] = 0.5 * m2[:, 0]

    if mode == 'd' or (mode != 'x' and nf1 == nf2):
        nx = min(nf1, nf2)
        d = 2 * np.sum(lpcar2rr(ar1[:nx, :]) * m2[:nx, :], axis=1) * \
            ((ar1[:nx, 0] / ar2[:nx, 0])**2)
    else:
        d = 2 * (lpcar2rr(ar1) @ m2.T) * \
            ((ar1[:, 0:1] / ar2[:, 0:1].T)**2)

    if 'e' not in mode:
        d = np.log(d)

    return d

"""
Numba-optimized utility functions for GLOAT toolkit.
"""

import numpy as np
from numba import jit


@jit(nopython=True)
def lpc_levinson_durbin_numba(r, order):
    """
    Levinson-Durbin recursion for LPC computation (Numba JIT compiled).

    Expected speedup: 10-20x

    Parameters
    ----------
    r : ndarray
        Autocorrelation coefficients [0 to order]
    order : int
        LPC order

    Returns
    -------
    a : ndarray
        LPC coefficients [1, a1, a2, ..., a_order]
    e : float
        Prediction error
    """
    a = np.zeros(order + 1)
    a[0] = 1.0
    e = r[0]

    if e == 0:
        return a, 0.0

    for i in range(1, order + 1):
        # Compute reflection coefficient
        lambda_i = r[i]
        for j in range(1, i):
            lambda_i -= a[j] * r[i - j]

        if e == 0:
            break

        lambda_i /= e

        # Update coefficients
        a_old = a.copy()
        a[i] = lambda_i
        for j in range(1, i):
            a[j] = a_old[j] - lambda_i * a_old[i - j]

        # Update error
        e = e * (1 - lambda_i * lambda_i)

        if e <= 0:
            break

    return a, e


@jit(nopython=True)
def compute_autocorrelation_numba(x, order):
    """
    Compute autocorrelation coefficients (Numba JIT compiled).

    Expected speedup: 10-15x

    Parameters
    ----------
    x : ndarray
        Input signal
    order : int
        Maximum lag

    Returns
    -------
    r : ndarray
        Autocorrelation [r(0), r(1), ..., r(order)]
    """
    n = len(x)
    r = np.zeros(order + 1)

    for lag in range(order + 1):
        sum_val = 0.0
        for i in range(n - lag):
            sum_val += x[i] * x[i + lag]
        r[lag] = sum_val

    return r


@jit(nopython=True)
def lpc_analysis_numba(segment, order):
    """
    Complete LPC analysis on a segment (Numba JIT compiled).

    Parameters
    ----------
    segment : ndarray
        Windowed signal segment
    order : int
        LPC order

    Returns
    -------
    A : ndarray
        LPC filter coefficients
    e : float
        Prediction error
    """
    # Compute autocorrelation
    r = compute_autocorrelation_numba(segment, order)

    # Levinson-Durbin recursion
    A, e = lpc_levinson_durbin_numba(r, order)

    return A, e


@jit(nopython=True)
def apply_lpc_filter_numba(A, segment):
    """
    Apply LPC inverse filter to segment (Numba JIT compiled).

    Parameters
    ----------
    A : ndarray
        LPC coefficients
    segment : ndarray
        Input segment

    Returns
    -------
    residual : ndarray
        Filtered output (residual)
    """
    order = len(A) - 1
    n = len(segment)
    residual = np.zeros(n)

    for i in range(n):
        sum_val = A[0] * segment[i]
        for j in range(1, min(i + 1, order + 1)):
            sum_val += A[j] * segment[i - j]
        residual[i] = sum_val

    return residual

# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
"""
Cython-optimized utility functions for LPC analysis.
"""

import numpy as np
cimport numpy as cnp
cimport cython
from libc.math cimport sqrt

ctypedef cnp.float64_t DTYPE_t


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def compute_autocorrelation_cython(double[:] x, int order):
    """
    Compute autocorrelation coefficients (Cython).

    Expected speedup: 10-15x

    Parameters
    ----------
    x : array
        Input signal
    order : int
        Maximum lag

    Returns
    -------
    r : ndarray
        Autocorrelation [r(0), r(1), ..., r(order)]
    """
    cdef int n = x.shape[0]
    cdef cnp.ndarray[DTYPE_t, ndim=1] r = np.zeros(order + 1, dtype=np.float64)
    cdef double[:] r_view = r

    cdef int lag, i
    cdef double sum_val

    with nogil:
        for lag in range(order + 1):
            sum_val = 0.0
            for i in range(n - lag):
                sum_val += x[i] * x[i + lag]
            r_view[lag] = sum_val

    return r


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def lpc_levinson_durbin_cython(double[:] r, int order):
    """
    Levinson-Durbin recursion for LPC (Cython).

    Expected speedup: 10-20x

    Parameters
    ----------
    r : array
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
    cdef cnp.ndarray[DTYPE_t, ndim=1] a = np.zeros(order + 1, dtype=np.float64)
    cdef cnp.ndarray[DTYPE_t, ndim=1] a_old = np.zeros(order + 1, dtype=np.float64)
    cdef double[:] a_view = a
    cdef double[:] a_old_view = a_old

    a_view[0] = 1.0
    cdef double e = r[0]

    cdef int i, j
    cdef double lambda_i

    if e == 0:
        return a, 0.0

    with nogil:
        for i in range(1, order + 1):
            # Compute reflection coefficient
            lambda_i = r[i]
            for j in range(1, i):
                lambda_i -= a_view[j] * r[i - j]

            if e == 0:
                break

            lambda_i /= e

            # Copy old coefficients
            for j in range(i + 1):
                a_old_view[j] = a_view[j]

            # Update coefficients
            a_view[i] = lambda_i
            for j in range(1, i):
                a_view[j] = a_old_view[j] - lambda_i * a_old_view[i - j]

            # Update error
            e = e * (1.0 - lambda_i * lambda_i)

            if e <= 0:
                break

    return a, e


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def apply_lpc_filter_cython(double[:] A, double[:] segment):
    """
    Apply LPC inverse filter to segment (Cython).

    Parameters
    ----------
    A : array
        LPC coefficients
    segment : array
        Input segment

    Returns
    -------
    residual : ndarray
        Filtered output (residual)
    """
    cdef int order = A.shape[0] - 1
    cdef int n = segment.shape[0]
    cdef cnp.ndarray[DTYPE_t, ndim=1] residual = np.zeros(n, dtype=np.float64)
    cdef double[:] res_view = residual

    cdef int i, j
    cdef double sum_val

    with nogil:
        for i in range(n):
            sum_val = A[0] * segment[i]
            for j in range(1, min(i + 1, order + 1)):
                sum_val += A[j] * segment[i - j]
            res_view[i] = sum_val

    return residual


@cython.boundscheck(False)
@cython.wraparound(False)
def lpc_analysis_cython(double[:] segment, int order):
    """
    Complete LPC analysis on a segment (Cython).

    Parameters
    ----------
    segment : array
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
    cdef cnp.ndarray[DTYPE_t, ndim=1] r = compute_autocorrelation_cython(segment, order)

    # Levinson-Durbin recursion
    A, e = lpc_levinson_durbin_cython(r, order)

    return A, e

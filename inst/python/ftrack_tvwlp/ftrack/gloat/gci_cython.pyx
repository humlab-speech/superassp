# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
"""
Cython-optimized GCI detection functions.

Replaces critical loops in SEDREAMS GCI detection with compiled C code.
"""

import numpy as np
cimport numpy as cnp
cimport cython
from libc.math cimport round, fabs
from libc.stdlib cimport malloc, free

ctypedef cnp.float64_t DTYPE_t


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def compute_mean_based_signal_cython(
    double[:] wave,
    double[:] black_win,
    int half_l
):
    """
    Compute mean-based signal using sliding window (Cython).

    Expected speedup: 200-300x over Python loops.

    Parameters
    ----------
    wave : array
        Input signal
    black_win : array
        Blackman window
    half_l : int
        Half window length

    Returns
    -------
    mean_based_signal : ndarray
        Mean-based signal
    """
    cdef int n_samples = wave.shape[0]
    cdef cnp.ndarray[DTYPE_t, ndim=1] mean_based_signal = np.zeros(n_samples, dtype=np.float64)
    cdef double[:] mbs_view = mean_based_signal

    cdef int win_len = 2 * half_l + 1
    cdef double win_sum = 0.0
    cdef int m, i
    cdef double sum_val

    # Precompute window sum
    with nogil:
        for i in range(win_len):
            win_sum += black_win[i]

        # Sliding window
        for m in range(half_l, n_samples - half_l):
            sum_val = 0.0
            for i in range(win_len):
                sum_val += wave[m - half_l + i] * black_win[i]
            mbs_view[m] = sum_val / win_sum

    return mean_based_signal


@cython.boundscheck(False)
@cython.wraparound(False)
def find_extrema_cython(double[:] signal):
    """
    Find local minima and maxima (Cython).

    Expected speedup: 50-100x over Python loops.

    Parameters
    ----------
    signal : array
        Input signal

    Returns
    -------
    maxima : ndarray
        Indices of local maxima
    minima : ndarray
        Indices of local minima
    """
    cdef int n = signal.shape[0]
    cdef int* maxima_buf = <int*>malloc(n * sizeof(int))
    cdef int* minima_buf = <int*>malloc(n * sizeof(int))
    cdef int max_count = 0
    cdef int min_count = 0
    cdef int m
    cdef cnp.ndarray[cnp.int64_t, ndim=1] maxima
    cdef cnp.ndarray[cnp.int64_t, ndim=1] minima

    try:
        with nogil:
            for m in range(1, n - 1):
                if signal[m] > signal[m - 1] and signal[m] > signal[m + 1]:
                    maxima_buf[max_count] = m
                    max_count += 1
                elif signal[m] < signal[m - 1] and signal[m] < signal[m + 1]:
                    minima_buf[min_count] = m
                    min_count += 1

        # Copy to numpy arrays
        maxima = np.zeros(max_count, dtype=np.int64)
        minima = np.zeros(min_count, dtype=np.int64)

        for m in range(max_count):
            maxima[m] = maxima_buf[m]
        for m in range(min_count):
            minima[m] = minima_buf[m]

        return maxima, minima

    finally:
        free(maxima_buf)
        free(minima_buf)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def compute_gci_from_residual_cython(
    double[:] res,
    cnp.int64_t[:] minis,
    cnp.int64_t[:] maxis,
    double ratio_gci
):
    """
    Detect GCI locations from LPC residual (Cython).

    Expected speedup: 10-20x over Python loops.

    Parameters
    ----------
    res : array
        LPC residual signal
    minis : array
        Minima indices
    maxis : array
        Maxima indices
    ratio_gci : float
        Median ratio of GCI position

    Returns
    -------
    gci : ndarray
        GCI sample locations
    """
    cdef int n_cycles = minis.shape[0]
    cdef int* gci_buf = <int*>malloc(n_cycles * sizeof(int))
    cdef int gci_count = 0

    cdef int k, interv, start, stop, idx, max_idx
    cdef double alpha, max_val
    cdef int res_len = res.shape[0]
    cdef cnp.ndarray[cnp.int64_t, ndim=1] gci

    try:
        with nogil:
            for k in range(n_cycles):
                if k >= maxis.shape[0]:
                    break

                interv = maxis[k] - minis[k]

                alpha = ratio_gci - 0.25
                start = minis[k] + <int>round(alpha * interv)

                alpha = ratio_gci + 0.35
                stop = minis[k] + <int>round(alpha * interv)

                if start < 0:
                    start = 0
                if stop >= res_len:
                    stop = res_len - 1

                if start >= stop:
                    continue

                # Find maximum in window
                max_val = res[start]
                max_idx = start

                for idx in range(start, stop + 1):
                    if res[idx] > max_val:
                        max_val = res[idx]
                        max_idx = idx

                gci_buf[gci_count] = max_idx
                gci_count += 1

        # Copy to numpy array
        gci = np.zeros(gci_count, dtype=np.int64)
        for k in range(gci_count):
            gci[k] = gci_buf[k]

        return gci

    finally:
        free(gci_buf)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def compute_relative_gci_positions_cython(
    double[:] res,
    cnp.int64_t[:] posis,
    cnp.int64_t[:] minis,
    cnp.int64_t[:] maxis
):
    """
    Compute relative GCI positions for median estimation (Cython).

    Parameters
    ----------
    res : array
        Absolute LPC residual
    posis : array
        Candidate peak positions
    minis : array
        Minima from mean-based signal
    maxis : array
        Maxima from mean-based signal

    Returns
    -------
    rel_posis : ndarray
        Relative positions
    """
    cdef int n_posis = posis.shape[0]
    cdef double* rel_buf = <double*>malloc(n_posis * sizeof(double))
    cdef int rel_count = 0

    cdef int k, j, pos, interv
    cdef double min_dist, dist, rel_pos
    cdef int n_minis = minis.shape[0]
    cdef int n_maxis = maxis.shape[0]
    cdef cnp.ndarray[DTYPE_t, ndim=1] rel_posis

    try:
        with nogil:
            for k in range(n_posis):
                # Find closest minimum
                min_dist = fabs(<double>(minis[0] - posis[k]))
                pos = 0

                for j in range(1, n_minis):
                    dist = fabs(<double>(minis[j] - posis[k]))
                    if dist < min_dist:
                        min_dist = dist
                        pos = j

                if pos >= n_maxis:
                    continue

                interv = maxis[pos] - minis[pos]
                if interv > 0:
                    rel_pos = <double>(posis[k] - minis[pos]) / <double>interv
                    rel_buf[rel_count] = rel_pos
                    rel_count += 1

        # Copy to numpy array
        rel_posis = np.zeros(rel_count, dtype=np.float64)
        for k in range(rel_count):
            rel_posis[k] = rel_buf[k]

        return rel_posis

    finally:
        free(rel_buf)

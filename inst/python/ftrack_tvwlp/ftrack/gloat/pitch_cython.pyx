# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True, initializedcheck=False
"""
Cython-optimized SRH pitch tracking.

This replaces the bottleneck _srh_estimate_pitch function with compiled C code.
Expected speedup: 10-50x on pitch estimation loops.
"""

import numpy as np
cimport numpy as cnp
cimport cython
from libc.math cimport round, floor, sqrt, fabs
from libc.stdlib cimport malloc, free

ctypedef cnp.float64_t DTYPE_t
ctypedef cnp.int64_t INT_t


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline double compute_srh_value(
    double[:] signal,
    int start_idx,
    int T_candidate,
    int K,
    int L
) nogil:
    """
    Compute SRH value for a single pitch candidate.

    Parameters
    ----------
    signal : memoryview
        Input signal
    start_idx : int
        Starting index in signal
    T_candidate : int
        Period length to test (samples)
    K : int
        Number of harmonics
    L : int
        Window size

    Returns
    -------
    srh_val : double
        SRH correlation value
    """
    cdef double sum_val = 0.0
    cdef double norm = 0.0
    cdef int h, tau, idx1, idx2
    cdef int n_samples = signal.shape[0]
    cdef double val1, val2

    # Summation of residual harmonics
    for h in range(1, K + 1):
        tau = h * T_candidate

        for idx1 in range(start_idx, min(start_idx + L, n_samples)):
            idx2 = idx1 + tau

            if idx2 >= n_samples:
                break

            val1 = signal[idx1]
            val2 = signal[idx2]
            sum_val += val1 * val2
            norm += val1 * val1

    if norm > 0:
        return sum_val / sqrt(norm)
    else:
        return 0.0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def srh_estimate_pitch_cython(
    double[:] signal,
    int fs,
    double f0_min,
    double f0_max,
    int K=5,
    int hop_length=160
):
    """
    Cython-optimized SRH pitch estimation.

    Parameters
    ----------
    signal : array
        Speech signal
    fs : int
        Sampling rate
    f0_min : float
        Minimum F0 (Hz)
    f0_max : float
        Maximum F0 (Hz)
    K : int
        Number of harmonics (default: 5)
    hop_length : int
        Frame shift (default: 160 samples = 10ms @ 16kHz)

    Returns
    -------
    f0 : ndarray
        F0 contour (Hz)
    vuv : ndarray
        Voice/unvoiced decision (1=voiced, 0=unvoiced)
    srh_val : ndarray
        SRH correlation values
    """
    cdef int n_samples = signal.shape[0]
    cdef int n_frames = (n_samples - 1) // hop_length + 1
    cdef int T_min = <int>round(fs / f0_max)
    cdef int T_max = <int>round(fs / f0_min)
    cdef int L = 2 * T_max  # Window size

    # Output arrays
    cdef cnp.ndarray[DTYPE_t, ndim=1] f0_out = np.zeros(n_frames, dtype=np.float64)
    cdef cnp.ndarray[DTYPE_t, ndim=1] vuv_out = np.zeros(n_frames, dtype=np.float64)
    cdef cnp.ndarray[DTYPE_t, ndim=1] srh_out = np.zeros(n_frames, dtype=np.float64)

    cdef double[:] f0_view = f0_out
    cdef double[:] vuv_view = vuv_out
    cdef double[:] srh_view = srh_out

    # Temporary variables
    cdef int frame_idx, start_idx, T_candidate, T_best
    cdef double srh_current, srh_best
    cdef double threshold = 0.3

    # Process each frame
    with nogil:
        for frame_idx in range(n_frames):
            start_idx = frame_idx * hop_length

            if start_idx + L >= n_samples:
                # Not enough samples
                f0_view[frame_idx] = 0.0
                vuv_view[frame_idx] = 0.0
                srh_view[frame_idx] = 0.0
                continue

            # Search over period range
            srh_best = -999999.0
            T_best = T_min

            for T_candidate in range(T_min, T_max + 1):
                srh_current = compute_srh_value(signal, start_idx, T_candidate, K, L)

                if srh_current > srh_best:
                    srh_best = srh_current
                    T_best = T_candidate

            # Store results
            srh_view[frame_idx] = srh_best

            if srh_best > threshold:
                f0_view[frame_idx] = <double>fs / <double>T_best
                vuv_view[frame_idx] = 1.0
            else:
                f0_view[frame_idx] = 0.0
                vuv_view[frame_idx] = 0.0

    return f0_out, vuv_out, srh_out


@cython.boundscheck(False)
@cython.wraparound(False)
def median_filter_f0_cython(
    double[:] f0,
    double[:] vuv,
    int window_size=5
):
    """
    Cython-optimized median filtering for F0 tracks.

    Parameters
    ----------
    f0 : array
        F0 contour
    vuv : array
        Voice/unvoiced flags
    window_size : int
        Median filter window size (default: 5)

    Returns
    -------
    f0_filtered : ndarray
        Filtered F0 contour
    """
    cdef int n_frames = f0.shape[0]
    cdef cnp.ndarray[DTYPE_t, ndim=1] f0_filtered = np.zeros(n_frames, dtype=np.float64)
    cdef double[:] f0_filt_view = f0_filtered

    cdef int i, j, k, half_win, count
    cdef double* window_buf
    cdef double temp, median_val

    half_win = window_size // 2
    window_buf = <double*>malloc(window_size * sizeof(double))

    try:
        for i in range(n_frames):
            if vuv[i] == 0:
                f0_filt_view[i] = 0.0
                continue

            # Collect voiced frames in window
            count = 0
            for j in range(max(0, i - half_win), min(n_frames, i + half_win + 1)):
                if vuv[j] > 0:
                    window_buf[count] = f0[j]
                    count += 1

            if count == 0:
                f0_filt_view[i] = 0.0
                continue

            # Simple bubble sort for median (small window)
            for j in range(count):
                for k in range(j + 1, count):
                    if window_buf[j] > window_buf[k]:
                        temp = window_buf[j]
                        window_buf[j] = window_buf[k]
                        window_buf[k] = temp

            # Get median
            if count % 2 == 1:
                median_val = window_buf[count // 2]
            else:
                median_val = (window_buf[count // 2 - 1] + window_buf[count // 2]) / 2.0

            f0_filt_view[i] = median_val

    finally:
        free(window_buf)

    return f0_filtered

# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True
# cython: initializedcheck=False, nonecheck=False
"""
Ultra-fast Cython implementations for metric computation.

These implementations provide significant speedups for CPU inference by:
- Using static typing and C-level performance
- Parallel processing with OpenMP (if available)
- Avoiding Python interpreter overhead
- Direct memory access

Expected speedup: 5-20x for CPU operations
"""

import numpy as np
cimport numpy as cnp
cimport cython
from cython.parallel import prange
from libc.math cimport fabs

# Initialize numpy C API
cnp.import_array()


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def compute_stat_scores_cython(cnp.ndarray[cnp.int64_t, ndim=1] preds,
                                cnp.ndarray[cnp.int64_t, ndim=1] target):
    """Ultra-fast Cython implementation for computing TP, FP, TN, FN.

    This is significantly faster than the PyTorch version for CPU tensors.

    Parameters
    ----------
    preds : np.ndarray
        Binary predictions (0 or 1), shape (n,)
    target : np.ndarray
        Binary targets (0 or 1), shape (n,)

    Returns
    -------
    tp, fp, tn, fn : int
        True positive, false positive, true negative, false negative counts
    """
    cdef long long n = preds.shape[0]
    cdef long long i
    cdef long long tp = 0, fp = 0, tn = 0, fn = 0
    cdef long long p, t

    # Parallel loop for large arrays (OpenMP if available)
    # For small arrays, overhead dominates so we use nogil instead
    if n > 10000:
        # Parallel reduction
        for i in prange(n, nogil=True):
            p = preds[i]
            t = target[i]

            if p == 1:
                if t == 1:
                    tp += 1
                else:
                    fp += 1
            else:
                if t == 0:
                    tn += 1
                else:
                    fn += 1
    else:
        # Sequential for small arrays
        with nogil:
            for i in range(n):
                p = preds[i]
                t = target[i]

                if p == 1:
                    if t == 1:
                        tp += 1
                    else:
                        fp += 1
                else:
                    if t == 0:
                        tn += 1
                    else:
                        fn += 1

    return tp, fp, tn, fn


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def compute_stat_scores_2d_cython(cnp.ndarray[cnp.int64_t, ndim=2] preds,
                                   cnp.ndarray[cnp.int64_t, ndim=2] target):
    """Cython implementation for 2D stat scores (multiple thresholds).

    Parameters
    ----------
    preds : np.ndarray
        Binary predictions, shape (num_thresholds, n)
    target : np.ndarray
        Binary targets, shape (num_thresholds, n)

    Returns
    -------
    tp, fp, tn, fn : np.ndarray
        Arrays of shape (num_thresholds,) with counts for each threshold
    """
    cdef long long num_thresholds = preds.shape[0]
    cdef long long n = preds.shape[1]
    cdef long long i, j
    cdef cnp.ndarray[cnp.int64_t, ndim=1] tp_arr = np.zeros(num_thresholds, dtype=np.int64)
    cdef cnp.ndarray[cnp.int64_t, ndim=1] fp_arr = np.zeros(num_thresholds, dtype=np.int64)
    cdef cnp.ndarray[cnp.int64_t, ndim=1] tn_arr = np.zeros(num_thresholds, dtype=np.int64)
    cdef cnp.ndarray[cnp.int64_t, ndim=1] fn_arr = np.zeros(num_thresholds, dtype=np.int64)
    cdef long long p, t

    # Parallel over thresholds
    for i in prange(num_thresholds, nogil=True):
        for j in range(n):
            p = preds[i, j]
            t = target[i, j]

            if p == 1:
                if t == 1:
                    tp_arr[i] += 1
                else:
                    fp_arr[i] += 1
            else:
                if t == 0:
                    tn_arr[i] += 1
                else:
                    fn_arr[i] += 1

    return tp_arr, fp_arr, tn_arr, fn_arr


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def compute_mae_cython(cnp.ndarray[cnp.float32_t, ndim=1] preds,
                       cnp.ndarray[cnp.float32_t, ndim=1] target,
                       cnp.ndarray[cnp.float32_t, ndim=1] weights=None):
    """Ultra-fast Mean Absolute Error computation.

    Parameters
    ----------
    preds : np.ndarray
        Predictions, shape (n,)
    target : np.ndarray
        Targets, shape (n,)
    weights : np.ndarray, optional
        Weights for each sample, shape (n,)

    Returns
    -------
    mae : float
        Mean absolute error
    """
    cdef long long n = preds.shape[0]
    cdef long long i
    cdef double sum_abs_error = 0.0
    cdef double total_weight = 0.0
    cdef float w

    if weights is not None:
        # Weighted MAE
        with nogil:
            for i in range(n):
                w = weights[i]
                sum_abs_error += fabs(preds[i] - target[i]) * w
                total_weight += w
    else:
        # Unweighted MAE
        with nogil:
            for i in range(n):
                sum_abs_error += fabs(preds[i] - target[i])
            total_weight = n

    if total_weight > 0:
        return sum_abs_error / total_weight
    else:
        return 0.0


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def binarize_with_hysteresis_cython(cnp.ndarray[cnp.float32_t, ndim=1] scores,
                                     float onset,
                                     float offset):
    """Ultra-fast binarization with hysteresis thresholding.

    This is used in the pipeline for VAD post-processing.

    Parameters
    ----------
    scores : np.ndarray
        Frame-level scores, shape (num_frames,)
    onset : float
        Onset threshold (start speech)
    offset : float
        Offset threshold (end speech)

    Returns
    -------
    binary : np.ndarray
        Binary array, shape (num_frames,)
    """
    cdef long long num_frames = scores.shape[0]
    cdef cnp.ndarray[cnp.int8_t, ndim=1] binary = np.zeros(num_frames, dtype=np.int8)
    cdef long long i
    cdef bint active = False
    cdef float score

    with nogil:
        for i in range(num_frames):
            score = scores[i]

            if not active:
                # Start new segment if score >= onset
                if score >= onset:
                    active = True
                    binary[i] = 1
            else:
                # Continue segment until score < offset
                if score >= offset:
                    binary[i] = 1
                else:
                    active = False

    return binary


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def remove_short_segments_cython(cnp.ndarray[cnp.int8_t, ndim=1] binary,
                                  int min_duration_frames,
                                  int value):
    """Remove short segments of a specific value.

    Parameters
    ----------
    binary : np.ndarray
        Binary array (0 or 1), shape (num_frames,)
    min_duration_frames : int
        Minimum duration in frames
    value : int
        Value to filter (0 or 1)

    Returns
    -------
    filtered : np.ndarray
        Filtered binary array
    """
    cdef long long num_frames = binary.shape[0]
    cdef cnp.ndarray[cnp.int8_t, ndim=1] result = binary.copy()
    cdef long long i = 0
    cdef long long start, end, duration

    with nogil:
        while i < num_frames:
            if result[i] == value:
                # Found start of segment
                start = i
                while i < num_frames and result[i] == value:
                    i += 1
                end = i

                # Check duration
                duration = end - start
                if duration < min_duration_frames:
                    # Remove short segment
                    for j in range(start, end):
                        result[j] = 1 - value
            else:
                i += 1

    return result


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def apply_binarization_full_cython(cnp.ndarray[cnp.float32_t, ndim=1] scores,
                                    float onset,
                                    float offset,
                                    float min_duration_on,
                                    float min_duration_off,
                                    float frame_duration):
    """Complete binarization pipeline in Cython.

    This combines hysteresis + short segment removal in one optimized function.

    Parameters
    ----------
    scores : np.ndarray
        Frame-level scores
    onset : float
        Onset threshold
    offset : float
        Offset threshold
    min_duration_on : float
        Minimum speech duration in seconds
    min_duration_off : float
        Minimum gap duration in seconds
    frame_duration : float
        Frame duration in seconds

    Returns
    -------
    binary : np.ndarray
        Final binary array
    """
    cdef cnp.ndarray[cnp.int8_t, ndim=1] binary
    cdef int min_frames_on, min_frames_off

    # Apply hysteresis
    binary = binarize_with_hysteresis_cython(scores, onset, offset)

    # Remove short speech segments
    if min_duration_on > 0:
        min_frames_on = int(min_duration_on / frame_duration)
        binary = remove_short_segments_cython(binary, min_frames_on, 1)

    # Fill short gaps
    if min_duration_off > 0:
        min_frames_off = int(min_duration_off / frame_duration)
        binary = remove_short_segments_cython(binary, min_frames_off, 0)

    return binary

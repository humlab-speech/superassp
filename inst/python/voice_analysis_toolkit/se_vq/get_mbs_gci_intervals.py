"""GCI search interval detection from Mean-Based Signal."""

import numpy as np
from scipy.signal import find_peaks


def get_mbs_gci_intervals(mbs, fs, T0mean, F0max=500):
    """
    Detect intervals for searching GCI locations from mean-based signal.

    Parameters
    ----------
    mbs : array_like
        Mean-based signal
    fs : float
        Sampling frequency (Hz)
    T0mean : float or array_like
        Mean period length in samples
    F0max : float, optional
        Maximum F0 value (Hz), default: 500

    Returns
    -------
    interval : ndarray
        N x 2 array where each row contains [start, stop] sample indices
        for GCI search intervals
    """
    mbs = np.asarray(mbs).flatten()

    if np.isscalar(T0mean):
        T0mean_arr = np.full(len(mbs), T0mean)
    else:
        T0mean_arr = np.asarray(T0mean)

    F0max = F0max * 2
    T0max = int(np.round(fs / F0max))

    # Find negative peaks in MBS
    _, idx = find_peaks(-mbs, distance=T0max)
    N = len(idx)

    search_rate = 0.28
    search_left_rate = 0.01
    interval = np.zeros((N, 2), dtype=int)

    # Define search intervals around each negative peak
    for n in range(N):
        if len(T0mean_arr) > 1:
            start = idx[n] - int(np.round(T0mean_arr[idx[n]] * search_left_rate))
            stop = idx[n] + int(np.round(T0mean_arr[idx[n]] * search_rate))
        else:
            start = idx[n] - int(np.round(T0mean_arr[0] * search_left_rate))
            stop = idx[n] + int(np.round(T0mean_arr[0] * search_rate))

        if start < 0:
            start = 0

        # Check bounds
        if stop >= len(mbs) and start < len(mbs):
            stop = len(mbs) - 1
        elif stop >= len(mbs) and start >= len(mbs):
            break

        interval[n, 0] = start
        interval[n, 1] = stop

    # Remove any zero rows
    interval = interval[:n + 1]

    return interval

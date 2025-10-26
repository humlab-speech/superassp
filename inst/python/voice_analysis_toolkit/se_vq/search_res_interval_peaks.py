"""Search for GCI candidates in residual signal within intervals."""

import numpy as np


def search_res_interval_peaks(res, interval, Ncand=5):
    """
    Search for N candidate GCI peaks in residual within each search interval.

    Parameters
    ----------
    res : array_like
        LP-residual signal
    interval : ndarray
        N x 2 array of [start, stop] indices for search intervals
    Ncand : int, optional
        Number of candidate peaks to find per interval, default: 5

    Returns
    -------
    gci : ndarray
        N x Ncand array of candidate GCI locations (sample indices)
    rel_amp : ndarray
        N x Ncand array of relative amplitudes (normalized inverse amplitudes)

    Notes
    -----
    For each interval, the Ncand highest peaks are selected. If an interval
    contains fewer than Ncand samples, only the maximum peak is used.
    """
    res = np.asarray(res).flatten()
    interval = np.asarray(interval, dtype=int)

    N = interval.shape[0]
    gci = np.zeros((N, Ncand))
    rel_amp = np.zeros((N, Ncand))

    for n in range(N):
        start = interval[n, 0]
        stop = interval[n, 1]

        if stop <= start:
            gci[n, :] = np.nan
        elif stop - start < Ncand:
            # Not enough samples, just use the maximum
            seg = res[start:stop + 1]
            idx = np.argmax(seg)
            gci_cur = idx + start

            gci[n, :] = gci_cur
            rel_amp[n, :] = 0
        else:
            # Find top Ncand peaks
            seg = res[start:stop + 1]
            idx_sorted = np.argsort(seg)[::-1]  # Descending order
            amp_sorted = seg[idx_sorted]

            gci_cur = idx_sorted[:Ncand] + start
            gci[n, :] = gci_cur

            # Relative amplitude (normalized inverse)
            max_amp = amp_sorted[0]
            if max_amp > 0:
                rel_amp[n, :] = 1 - (amp_sorted[:Ncand] / max_amp)
            else:
                rel_amp[n, :] = 0

    return gci, rel_amp

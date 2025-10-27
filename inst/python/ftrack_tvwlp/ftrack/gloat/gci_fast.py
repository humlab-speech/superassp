"""
Optimized SEDREAMS GCI detection with automatic Numba fallback.

This module automatically uses Numba JIT if available, otherwise falls back to Python.
"""

import numpy as np
from scipy import signal

# Try to import Numba optimizations
try:
    from .gci_numba import (
        compute_mean_based_signal_numba,
        find_extrema_numba,
        compute_gci_from_residual_numba,
        compute_relative_gci_positions_numba
    )
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

from .utils import get_lpc_residual


def sedreams_gci_detection_fast(wave, fs, f0_mean, use_numba=True):
    """
    Fast SEDREAMS GCI detection with automatic Numba optimization.

    If Numba is available and use_numba=True, uses JIT-compiled functions
    for 50-300x speedup on critical loops. Otherwise falls back to Python.

    Parameters
    ----------
    wave : ndarray
        Speech signal
    fs : int
        Sampling frequency (Hz)
    f0_mean : float
        Average pitch value (Hz)
    use_numba : bool
        Use Numba if available (default: True)

    Returns
    -------
    gci : ndarray
        Glottal Closure Instant locations (samples, 0-indexed)
    mean_based_signal : ndarray
        Mean-based signal
    """
    wave = np.asarray(wave, dtype=np.float64)
    n_samples = len(wave)

    # Get LPC residual
    res = get_lpc_residual(wave, int(np.round(25 / 1000 * fs)),
                           int(np.round(5 / 1000 * fs)),
                           int(np.round(fs / 1000) + 2))

    # Calculate mean-based signal
    t0_mean = int(np.round(fs / f0_mean))
    half_l = int(np.round((1.6 * t0_mean) / 2))
    black_win = signal.windows.blackman(2 * half_l + 1)

    if HAS_NUMBA and use_numba:
        # Fast Numba version (200-300x faster)
        mean_based_signal = compute_mean_based_signal_numba(wave, black_win, half_l)
    else:
        # Python fallback
        mean_based_signal = np.zeros(n_samples, dtype=np.float64)
        for m in range(half_l, n_samples - half_l):
            vec = wave[m - half_l:m + half_l + 1].copy()
            vec = vec * black_win
            mean_based_signal[m] = np.mean(vec)

    # Remove low-frequency contents
    ws = 30 / (fs / 2)
    wp = 50 / (fs / 2)
    rp = 3
    rs = 60

    n, wn = signal.ellipord(wp, ws, rp, rs)
    b, a = signal.ellip(n, rp, rs, wn, btype='high')

    mean_based_signal = signal.filtfilt(b, a, mean_based_signal)
    mbs_max = np.max(np.abs(mean_based_signal))
    if mbs_max > 0:
        mean_based_signal = mean_based_signal / mbs_max

    # Detect minima and maxima
    if HAS_NUMBA and use_numba:
        # Fast Numba version (50-100x faster)
        pot_maxis, pot_minis = find_extrema_numba(mean_based_signal)
    else:
        # Python fallback
        pot_maxis = []
        pot_minis = []
        for m in range(1, len(mean_based_signal) - 1):
            if (mean_based_signal[m] > mean_based_signal[m - 1] and
                    mean_based_signal[m] > mean_based_signal[m + 1]):
                pot_maxis.append(m)
            elif (mean_based_signal[m] < mean_based_signal[m - 1] and
                  mean_based_signal[m] < mean_based_signal[m + 1]):
                pot_minis.append(m)
        pot_maxis = np.array(pot_maxis, dtype=np.int64)
        pot_minis = np.array(pot_minis, dtype=np.int64)

    if len(pot_maxis) == 0 or len(pot_minis) == 0:
        return np.array([], dtype=int), mean_based_signal

    # Ensure first max comes before first min
    while len(pot_maxis) > 0 and len(pot_minis) > 0 and pot_maxis[0] < pot_minis[0]:
        pot_maxis = pot_maxis[1:]

    # Ensure last min comes before last max
    while len(pot_maxis) > 0 and len(pot_minis) > 0 and pot_minis[-1] > pot_maxis[-1]:
        pot_minis = pot_minis[:-1]

    if len(pot_maxis) == 0 or len(pot_minis) == 0:
        return np.array([], dtype=int), mean_based_signal

    minis = pot_minis
    maxis = pot_maxis

    # Determine median position of GCIs
    res = res / np.max(np.abs(res))

    abs_res = np.abs(res)
    diff1 = np.diff(np.concatenate([[0], abs_res]))
    diff2 = np.diff(diff1 > 0)
    posis = np.where(diff2 == -1)[0]
    posis = posis[abs_res[posis] > 0.4]

    if len(posis) == 0:
        posis = np.where(abs_res > 0.4)[0]

    if HAS_NUMBA and use_numba:
        # Fast Numba version
        rel_posis = compute_relative_gci_positions_numba(abs_res, posis, minis, maxis)
    else:
        # Python fallback
        rel_posis = []
        for k in range(len(posis)):
            dists = np.abs(minis - posis[k])
            pos = np.argmin(dists)

            if pos >= len(maxis):
                continue

            interv = maxis[pos] - minis[pos]
            if interv > 0:
                rel_pos = (posis[k] - minis[pos]) / interv
                rel_posis.append(rel_pos)
        rel_posis = np.array(rel_posis)

    if len(rel_posis) == 0:
        ratio_gci = 0.5
    else:
        ratio_gci = np.median(rel_posis)

    # Detect GCIs from residual
    if HAS_NUMBA and use_numba:
        # Fast Numba version
        gci = compute_gci_from_residual_numba(res, minis, maxis, ratio_gci)
    else:
        # Python fallback
        gci_list = []
        for k in range(len(minis)):
            if k >= len(maxis):
                break

            interv = maxis[k] - minis[k]
            alpha = ratio_gci - 0.25
            start = minis[k] + int(np.round(alpha * interv))
            alpha = ratio_gci + 0.35
            stop = minis[k] + int(np.round(alpha * interv))

            start = max(0, start)
            stop = min(len(res) - 1, stop)

            if start >= stop:
                continue

            vec = res[start:stop + 1]
            if len(vec) == 0:
                continue

            max_idx = np.argmax(vec)
            gci_sample = start + max_idx
            gci_list.append(gci_sample)

        gci = np.array(gci_list, dtype=int)

    return gci, mean_based_signal


# Print status on import
if HAS_NUMBA:
    print("✓ Numba JIT available - using optimized GCI detection")
else:
    print("⚠ Numba not found - using Python fallback (slower)")
    print("  Install with: pip install numba")

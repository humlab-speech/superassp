"""
SEDREAMS (Speech Event Detection using the Residual Excitation And a Mean-based Signal)
GCI (Glottal Closure Instant) detection.

Based on:
T. Drugman, M. Thomas, J. Gudnason, P. Naylor, T. Dutoit,
"Detection of Glottal Closure Instants from Speech Signals: a Quantitative Review",
IEEE Transactions on Audio, Speech and Language Processing.

and

T.Drugman, T.Dutoit, "Glottal Closure and Opening Instant Detection from Speech Signals",
Interspeech09, Brighton, U.K, 2009
"""

import numpy as np
from scipy import signal
from .utils import get_lpc_residual


def sedreams_gci_detection(wave, fs, f0_mean):
    """
    SEDREAMS GCI detection algorithm.

    Parameters
    ----------
    wave : ndarray
        Speech signal
    fs : int
        Sampling frequency (Hz)
    f0_mean : float
        Average pitch value (Hz) for the speaker

    Returns
    -------
    gci : ndarray
        Glottal Closure Instant locations (in samples, 0-indexed)
    mean_based_signal : ndarray
        Mean-based signal used by SEDREAMS
    """
    wave = np.asarray(wave, dtype=np.float64)
    n_samples = len(wave)

    # Get LPC residual
    res = get_lpc_residual(wave, int(np.round(25 / 1000 * fs)),
                           int(np.round(5 / 1000 * fs)),
                           int(np.round(fs / 1000) + 2))
    res_orig = res.copy()

    # Calculate mean-based signal
    mean_based_signal = np.zeros(n_samples, dtype=np.float64)
    t0_mean = int(np.round(fs / f0_mean))

    half_l = int(np.round((1.6 * t0_mean) / 2))
    black_win = signal.windows.blackman(2 * half_l + 1)

    for m in range(half_l, n_samples - half_l):
        vec = wave[m - half_l:m + half_l + 1].copy()
        vec = vec * black_win
        mean_based_signal[m] = np.mean(vec)

    # Remove low-frequency contents of mean-based signal using elliptic high-pass filter
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

    # Detect minima and maxima of mean-based signal
    pot_maxis = []
    pot_minis = []

    for m in range(1, len(mean_based_signal) - 1):
        if (mean_based_signal[m] > mean_based_signal[m - 1] and
                mean_based_signal[m] > mean_based_signal[m + 1]):
            pot_maxis.append(m)
        elif (mean_based_signal[m] < mean_based_signal[m - 1] and
              mean_based_signal[m] < mean_based_signal[m + 1]):
            pot_minis.append(m)

    if len(pot_maxis) == 0 or len(pot_minis) == 0:
        # No GCIs detected
        return np.array([], dtype=int), mean_based_signal

    pot_maxis = np.array(pot_maxis)
    pot_minis = np.array(pot_minis)

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

    # Determine median position of GCIs within the cycle
    res = res / np.max(np.abs(res))

    # Find prominent peaks in residual
    # MATLAB: Posis=find(diff(diff([0 abs(res)])>0)==-1);
    # Posis=Posis(abs(res(Posis))>0.4);
    abs_res = np.abs(res)
    diff1 = np.diff(np.concatenate([[0], abs_res]))
    diff2 = np.diff(diff1 > 0)
    posis = np.where(diff2 == -1)[0]
    posis = posis[abs_res[posis] > 0.4]

    if len(posis) == 0:
        # Fallback: use all points above threshold
        posis = np.where(abs_res > 0.4)[0]

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

    if len(rel_posis) == 0:
        ratio_gci = 0.5  # Default
    else:
        ratio_gci = np.median(rel_posis)

    # Detect GCIs from residual using presence intervals
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

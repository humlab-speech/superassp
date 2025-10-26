"""Mean-Based Signal (MBS) calculation for SE_VQ algorithm."""

import numpy as np
from scipy.signal import blackman
from scipy.interpolate import interp1d
from ..general.signal_utils import zero_phase_hp_filt


def get_mbs(x, fs, T0mean):
    """
    Calculate the mean-based signal as used in SEDREAMS/SE_VQ.

    Parameters
    ----------
    x : array_like
        Speech signal
    fs : float
        Sampling frequency (Hz)
    T0mean : float or array_like
        Mean period length in samples (can be array for varying period)

    Returns
    -------
    mbs : ndarray
        Mean-based signal

    Notes
    -----
    The mean-based signal is computed by taking a local mean within windows
    centered at each sample point. The window size is proportional to T0mean.
    """
    x = np.asarray(x).flatten()

    if np.isscalar(T0mean):
        T0mean = np.full(len(x), T0mean)
    else:
        T0mean = np.asarray(T0mean)

    mbs = np.zeros(len(x))
    half_l = int(np.round((1.6 * T0mean[0]) / 2))

    step_exp = 3
    step = 2 ** step_exp

    # Calculate mean-based signal
    for m in range(half_l, len(x) - half_l, step):
        if len(T0mean) == 1:
            half_l = int(np.round((1.7 * T0mean[0]) / 2))
        else:
            half_l = int(np.round((1.7 * T0mean[m]) / 2))

        black_win = blackman(2 * half_l + 1)

        start = int(np.round(m - half_l))
        stop = int(np.round(m + half_l))

        if stop >= len(x):
            break

        if start > 0:
            vec = x[start:stop + 1]
            vec = vec * black_win
            mbs[m] = np.mean(vec)

    # Interpolate non-zero values
    t = np.where(mbs != 0)[0]
    if len(t) > 0:
        f = interp1d(t, mbs[t], kind='linear', bounds_error=False, fill_value=0)
        mbs = f(np.arange(len(x)))

    # High-pass filter
    mbs = zero_phase_hp_filt(mbs, fs, 70, 10, plots=False)

    # Normalize
    max_val = np.max(np.abs(mbs))
    if max_val > 0:
        mbs = mbs / max_val

    # Smooth (7-point moving average to match MATLAB smooth())
    # MATLAB's smooth() uses different edge handling
    from scipy.ndimage import uniform_filter1d
    mbs = uniform_filter1d(mbs, size=7, mode='reflect')

    return mbs

"""
Signal processing utility functions.

This module contains Python implementations of signal processing functions from
the Voice Analysis Toolkit, originally implemented in MATLAB by John Kane and
Thomas Drugman.
"""

import numpy as np
from scipy import signal as sp_signal
from scipy.signal import butter, buttord, filtfilt, lfilter, lfilter_zi, windows


def integrat(x, fs):
    """
    Apply a simple integrator to the signal.

    Parameters
    ----------
    x : array_like
        Input signal
    fs : float
        Sampling frequency (Hz)

    Returns
    -------
    y : ndarray
        Integrated signal
    """
    x = np.asarray(x).flatten()
    y = np.zeros(len(x))
    ts = 1.0 / fs

    y[0] = ts * x[0]

    # Difference equation
    for n in range(1, len(x)):
        y[n] = (ts * x[n]) + y[n - 1]

    return y


def zero_phase_hp_filt(x, fs, f_p, f_s, plots=False):
    """
    Apply forward and backward Butterworth high-pass filter for zero phase.

    Parameters
    ----------
    x : array_like
        Input signal
    fs : float
        Sampling frequency (Hz)
    f_p : float
        Pass-band frequency (Hz)
    f_s : float
        Stop-band frequency (Hz)
    plots : bool, optional
        If True, plot the filter's frequency response (default: False)

    Returns
    -------
    y : ndarray
        Filtered signal
    """
    x = np.asarray(x).flatten()

    rp = 0.5  # Passband ripple
    rs = 40   # Stopband attenuation

    wp = f_p / (fs / 2)
    ws = f_s / (fs / 2)

    # Design Butterworth filter
    n, wn = buttord(wp, ws, rp, rs)
    b, a = butter(n, wn, btype='high')

    # Apply zero-phase filtering
    y = filtfilt(b, a, x)

    if plots:
        import matplotlib.pyplot as plt
        w, h = sp_signal.freqz(b, a)
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.plot(w * fs / (2 * np.pi), 20 * np.log10(abs(h)))
        plt.ylabel('Magnitude (dB)')
        plt.xlabel('Frequency (Hz)')
        plt.title('Frequency Response')
        plt.grid()
        plt.subplot(2, 1, 2)
        plt.plot(w * fs / (2 * np.pi), np.angle(h))
        plt.ylabel('Phase (radians)')
        plt.xlabel('Frequency (Hz)')
        plt.grid()
        plt.show()

    return y


def zero_phase_lp_filt(x, fs, f_p, f_s, plots=False):
    """
    Apply forward and backward Butterworth low-pass filter for zero phase.

    Parameters
    ----------
    x : array_like
        Input signal
    fs : float
        Sampling frequency (Hz)
    f_p : float
        Pass-band frequency (Hz)
    f_s : float
        Stop-band frequency (Hz)
    plots : bool, optional
        If True, plot the filter's frequency response (default: False)

    Returns
    -------
    y : ndarray
        Filtered signal
    """
    x = np.asarray(x).flatten()

    rp = 0.5  # Passband ripple
    rs = 6    # Stopband attenuation

    wp = f_p / (fs / 2)
    ws = f_s / (fs / 2)

    # Design Butterworth filter
    n, wn = buttord(wp, ws, rp, rs)
    b, a = butter(n, wn, btype='low')

    # Apply zero-phase filtering
    y = filtfilt(b, a, x)

    if plots:
        import matplotlib.pyplot as plt
        w, h = sp_signal.freqz(b, a)
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.plot(w * fs / (2 * np.pi), 20 * np.log10(abs(h)))
        plt.ylabel('Magnitude (dB)')
        plt.xlabel('Frequency (Hz)')
        plt.title('Frequency Response')
        plt.grid()
        plt.subplot(2, 1, 2)
        plt.plot(w * fs / (2 * np.pi), np.angle(h))
        plt.ylabel('Phase (radians)')
        plt.xlabel('Frequency (Hz)')
        plt.grid()
        plt.show()

    return y


def get_lpc_residual(wave, L, shift, order, gci=None, type_vec=None, t0=None):
    """
    Extract LPC residual signal.

    Parameters
    ----------
    wave : array_like
        Speech signal
    L : int
        Window length in samples (typically 25ms)
    shift : int
        Window shift in samples (typically 5ms)
    order : int
        LPC order
    gci : array_like, optional
        GCI positions in samples
    type_vec : array_like, optional
        Voicing decisions (0=unvoiced, 1=voiced)
    t0 : array_like, optional
        Period values in samples

    Returns
    -------
    res : ndarray
        LPC residual signal
    lpc_coeff : ndarray
        LPC coefficients for each frame
    """
    from .lpc_utils import lpcauto

    wave = np.asarray(wave).flatten()
    do_ps = gci is not None

    if not do_ps:
        start = 0
        stop = start + L
        res = np.zeros(len(wave))
        lpc_coeff = np.zeros((order + 1, int(np.round(len(wave) / shift))))
        n = 0

        while stop < len(wave):
            segment = wave[start:stop]

            # Use Hanning window (sym=True to match MATLAB)
            segment_win = segment * windows.hann(len(segment), sym=True)

            # Compute LPC coefficients using lpcauto (matches MATLAB better)
            ar, e, _ = lpcauto(segment_win, order)
            a = ar[0]  # Take first frame
            e_val = e[0]

            if n < lpc_coeff.shape[1]:
                lpc_coeff[:, n] = a

            # Inverse filter
            inv = lfilter(a, 1, segment)

            # Normalize energy
            segment_energy = np.sum(segment**2)
            inv_energy = np.sum(inv**2)
            if inv_energy > 0 and segment_energy > 0:
                inv = inv * np.sqrt(segment_energy / inv_energy)

            res[start:stop] += inv

            start += shift
            stop += shift
            n += 1

        # Normalize final residual
        max_res = np.max(np.abs(res))
        if max_res > 0:
            res = res / max_res

    else:
        # Pitch-synchronous version would require CompleteFramingWithType
        # For now, implement basic version
        raise NotImplementedError(
            "Pitch-synchronous LPC residual extraction not yet implemented. "
            "Use the non-pitch-synchronous version by not providing gci parameter."
        )

    return res, lpc_coeff


def calc_residual(x, x_lpc, ord_lpc, gci):
    """
    Perform LPC analysis and inverse filtering (used in IAIF).

    Parameters
    ----------
    x : array_like
        Signal to be inverse filtered
    x_lpc : array_like
        Signal to perform LPC analysis on
    ord_lpc : int
        LPC prediction order
    gci : array_like
        Glottal closure instants (in samples)

    Returns
    -------
    vector_res : ndarray
        Residual after inverse filtering
    ar_lpc : ndarray
        LPC coefficients for each GCI frame
    e_lpc : ndarray
        Prediction errors for each GCI frame
    """
    from .lpc_utils import lpcauto
    from scipy.signal import lfiltic

    x = np.asarray(x).flatten()
    x_lpc = np.asarray(x_lpc).flatten()
    gci = np.asarray(gci, dtype=int).flatten()

    # Initialize
    vector_res = np.zeros(len(x))
    ze_lpc = np.zeros(ord_lpc)

    ar_lpc = np.zeros((ord_lpc + 1, len(gci)))
    e_lpc = np.zeros(len(gci))

    frame_res = None
    residual = None

    for n in range(len(gci)):
        # Get framing information
        if n > 0:
            T0_cur = gci[n] - gci[n - 1]
        else:
            T0_cur = gci[n + 1] - gci[n] if n + 1 < len(gci) else gci[n]

        start = gci[n] - T0_cur
        stop = gci[n] + T0_cur

        if start > 0 and stop <= len(x):
            # Do LPC analysis
            frame_lpc = x_lpc[start:stop]

            if len(frame_lpc) > ord_lpc * 1.5:
                # Use Hamming window (sym=True to match MATLAB)
                frame_wind = frame_lpc * windows.hamming(len(frame_lpc), sym=True)

                # Compute LPC coefficients using lpcauto
                ar, e, _ = lpcauto(frame_wind, ord_lpc)
                ar = np.real(ar[0])  # Take first frame
                e_val = np.real(e[0])

                ar_lpc[:, n] = ar
                e_lpc[n] = e_val

                # Do inverse filtering
                if n > 0 and frame_res is not None and residual is not None:
                    # Compute initial conditions using filtic (matches MATLAB better)
                    try:
                        last_input = np.flipud(frame_res)
                        last_output = np.flipud(residual)
                        ze_lpc = lfiltic(ar, np.sqrt(e_val), last_output[:ord_lpc], last_input[:ord_lpc])
                    except:
                        ze_lpc = np.zeros(ord_lpc)

                frame_res = x[start:stop]

                # Calculate LPC residual
                try:
                    if n > 0 and len(ze_lpc) == max(len(ar)-1, 0):
                        residual = lfilter(ar, np.sqrt(e_val), frame_res, zi=ze_lpc)[0]
                    else:
                        residual = lfilter(ar, np.sqrt(e_val), frame_res)
                except:
                    residual = frame_res.copy()
                    print('Problem with IAIF filtering')

                # Apply window
                residual_win = residual * windows.hamming(len(residual), sym=True)
                vector_res[start:stop] += residual_win

    return vector_res, ar_lpc, e_lpc

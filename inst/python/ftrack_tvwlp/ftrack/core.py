"""
Main formant tracking function using TVWLP/TVLP.
"""

import numpy as np
from scipy import signal
from scipy.ndimage import median_filter

from .tvlp import tvlp_l1, tvlp_l2, tvwlp_l1, tvwlp_l2
from .formants import tvlptoformants_akitofi
from .gloat.pitch import srh_pitch_tracking
from .gloat.gci import sedreams_gci_detection


def ftrack_tvwlp(s, fs, lptype='tvwlp_l2', nwin=None, nshift=None, p=8, q=3,
                 npeaks=3, preemp=0.97, fint=80, plot_flag=False):
    """
    Compute continuous time-varying formant tracks using TVWLP/TVLP.

    Parameters
    ----------
    s : ndarray
        Speech signal (will be converted to 8kHz)
    fs : int
        Sampling rate of input signal (Hz)
    lptype : str, optional
        LP method: 'tvwlp_l2', 'tvwlp_l1', 'tvlp_l2', 'tvlp_l1' (default: 'tvwlp_l2')
    nwin : int, optional
        Window/block size for TVLP analysis in samples at input fs (default: 1600 @ 8kHz = 200ms)
    nshift : int, optional
        Window shift for TVLP in samples at input fs (default: 1600 @ 8kHz = 200ms)
    p : int, optional
        LP/TVLP order (default: 8)
    q : int, optional
        Polynomial order (default: 3)
    npeaks : int, optional
        Number of formants to track (default: 3)
    preemp : float, optional
        Pre-emphasis factor (default: 0.97)
    fint : int, optional
        Formant interval in samples; fs/fint gives output formant rate (default: 80 = 10ms)
    plot_flag : bool, optional
        Plot spectrogram with formant tracks (default: False)

    Returns
    -------
    Fi : ndarray
        Formant tracks [npeaks x Ns_downsampled] in Hz
    Ak : ndarray
        Time-varying LPCs [p+1 x Ns] (before downsampling)
    """
    s = np.asarray(s, dtype=np.float64)
    fs_ref = 8000
    n1ms = int(np.floor(fs_ref / 1000))

    # Set defaults
    if nwin is None:
        nwin = 200 * n1ms  # 200ms
    if nshift is None:
        nshift = 200 * n1ms  # 200ms

    # Resample to 8kHz if needed
    if fs != fs_ref:
        n_samples_new = int(len(s) * fs_ref / fs)
        s = signal.resample(s, n_samples_new)
        nwin = int(np.floor(nwin * fs_ref / fs))
        nshift = int(np.floor(nshift * fs_ref / fs))
        fs = fs_ref

    # Pre-emphasis
    if preemp is not None and preemp != 0:
        shat = signal.lfilter([1, -preemp], 1, s)
    else:
        shat = s.copy()

    Ns = len(shat)

    lptype_lower = lptype.lower()

    # QCP weight function computation (for TVWLP methods)
    if lptype_lower in ['tvwlp_l2', 'tvwlp_l1']:
        # Parameters for QCP weight
        wmin = 0.00001  # Minimum value of weighting function
        DQ = 0.7  # Duration Quotient
        PQ = 0.05  # Position Quotient
        Nramp = 3  # Length of linear ramp

        f0_min = 60
        f0_max = 600

        # Pitch tracking
        f0, vuv, srh_val = srh_pitch_tracking(s, fs, f0_min, f0_max)

        # Compute mean F0 from voiced regions
        f0_tmp = f0 * vuv
        pos = f0_tmp != 0
        if np.sum(pos) > 0:
            f0_mean = np.mean(f0_tmp[pos])
        else:
            f0_mean = 150  # Default

        # GCI detection
        gc, mbs = sedreams_gci_detection(s, fs, f0_mean)

        # Compute QCP weight function
        w = _qcp_wt(s, p, DQ, PQ, wmin, Nramp, gc, fs)
    else:
        # Dummy weight (not used for TVLP methods)
        w = np.ones(len(s), dtype=np.float64)

    # Initialize output arrays
    # Pre-allocate for the middle portion
    Fi_list = []
    Ak_list = []

    # Sliding window analysis
    j = 0
    first_frame = True

    while j <= Ns - nwin:
        if j <= Ns - nwin - nshift:
            x = shat[j:j + nwin]
            wj = w[j:j + nwin]
        else:
            # Last window: take all remaining samples
            x = shat[j:]
            wj = w[j:]

        # Select method
        if lptype_lower == 'tvlp_l2':
            aki = tvlp_l2(x, p, q)
        elif lptype_lower == 'tvwlp_l2':
            aki = tvwlp_l2(x, p, q, wj)
        elif lptype_lower == 'tvlp_l1':
            aki = tvlp_l1(x, p, q)
        elif lptype_lower == 'tvwlp_l1':
            aki = tvwlp_l1(x, p, q, wj)
        else:
            raise ValueError(f"Unknown lptype: {lptype}")

        Nx = len(x)
        fi, ak = tvlptoformants_akitofi(aki, Nx, npeaks, fs)

        # Handle frame boundaries
        half_overlap = int(np.floor((nwin - nshift) / 2))

        if first_frame:
            # First frame: take from 0 to half_overlap
            Fi_list.append(fi[:half_overlap, :])
            Ak_list.append(ak[:, :half_overlap])
            first_frame = False

        # Always append the middle portion
        if j <= Ns - nwin - nshift:
            # Regular frame
            Fi_list.append(fi[half_overlap:half_overlap + nshift, :])
            Ak_list.append(ak[:, half_overlap:half_overlap + nshift])
        else:
            # Last frame: take remaining samples
            Fi_list.append(fi[nshift + half_overlap:, :])
            Ak_list.append(ak[:, nshift + half_overlap:])

        j += nshift

    # Concatenate all frames
    if len(Fi_list) > 0:
        Fi = np.vstack(Fi_list)
        Ak = np.hstack(Ak_list)
    else:
        Fi = np.zeros((Ns, npeaks))
        Ak = np.zeros((p + 1, Ns))

    # Pad to match signal length
    if Fi.shape[0] < Ns:
        pad_rows = Ns - Fi.shape[0]
        Fi = np.vstack([Fi, np.zeros((pad_rows, npeaks))])
        Ak = np.hstack([Ak, np.zeros((p + 1, pad_rows))])

    # Downsample formant tracks
    Fi = Fi[fint - 1::fint, :]  # MATLAB: fint:fint:end

    # Transpose to [npeaks x time]
    Fi = Fi.T

    # Apply median filter (5-tap) along time axis
    # scipy's median_filter with size=(1, 5) applies along axis 1
    Fi = median_filter(Fi, size=(1, 5), mode='reflect')

    # Plotting
    if plot_flag:
        import matplotlib.pyplot as plt
        n1ms_plot = int(np.floor(fs / 1000))
        f, t, Sxx = signal.spectrogram(s, fs, window=signal.windows.hamming(20 * n1ms_plot),
                                       nperseg=20 * n1ms_plot, noverlap=15 * n1ms_plot,
                                       nfft=1024)
        plt.figure(figsize=(12, 6))
        plt.pcolormesh(t, f / 1000, 20 * np.log10(np.abs(Sxx) + 1e-10), shading='gouraud')
        plt.ylabel('Frequency (kHz)')
        plt.xlabel('Time (s)')
        plt.colorbar(label='dB')

        # Overlay formant tracks
        t_formants = np.arange(fint, Ns + 1, fint) / fs
        for i in range(npeaks):
            plt.plot(t_formants, Fi[i, :] / 1000, '.k', markersize=2)
        plt.title(f'Spectrogram with Formant Tracks ({lptype})')
        plt.show()

    return Fi, Ak


def _qcp_wt(x, p, DQ, PQ, d, Nramp, gci_ins, fs):
    """
    Create a QCP (Quasi-Closed-Phase) weight function.

    Parameters
    ----------
    x : ndarray
        Signal
    p : int
        LPC order (for padding)
    DQ : float
        Duration Quotient (0 to 1)
    PQ : float
        Position Quotient (0 to 1)
    d : float
        Minimum value of weight function
    Nramp : int
        Length of linear ramp (samples)
    gci_ins : ndarray
        Glottal Closure Instants (samples, 0-indexed)
    fs : int
        Sampling frequency

    Returns
    -------
    w : ndarray
        Weight function (length N+p)
    """
    N = len(x)

    # Create ramps
    if Nramp > 0:
        up_ramp = np.linspace(d, 1, 2 + Nramp)[1:-1]
        down_ramp = up_ramp[::-1]
    else:
        up_ramp = np.array([])
        down_ramp = np.array([])

    # Adjust DQ and PQ if needed
    if DQ + PQ > 1:
        DQ = 1 - PQ

    # Initialize weight to minimum value
    w = d * np.ones(N + p, dtype=np.float64)

    if len(gci_ins) < 2:
        # No GCIs, return minimum weight
        return w[:N]

    # Process each pitch period
    for i in range(len(gci_ins) - 1):
        T = gci_ins[i + 1] - gci_ins[i]
        T1 = int(np.round(DQ * T))
        T2 = int(np.round(PQ * T))

        while T1 + T2 > T:
            T1 = T1 - 1
            if T1 < 0:
                break

        if T1 <= 0:
            continue

        # Set closed phase to 1
        start_idx = gci_ins[i] + T2
        end_idx = gci_ins[i] + T2 + T1

        if start_idx < len(w) and end_idx <= len(w):
            w[start_idx:end_idx] = 1.0

            # Apply ramps
            if Nramp > 0 and len(up_ramp) > 0:
                ramp_end = min(start_idx + Nramp, end_idx)
                w[start_idx:ramp_end] = up_ramp[:ramp_end - start_idx]

                ramp_start = max(end_idx - Nramp, start_idx)
                if ramp_start < end_idx:
                    w[ramp_start:end_idx] = down_ramp[:end_idx - ramp_start]

    # Handle last period (if there's space after last GCI)
    i = len(gci_ins) - 1
    Nend = N - (T2 + gci_ins[i])

    if T2 + gci_ins[i] < N:
        T = gci_ins[i] - gci_ins[i - 1] if i > 0 else Nend + T2

        if T1 + T2 < Nend:
            start_idx = gci_ins[i] + T2
            end_idx = gci_ins[i] + T2 + T1

            if start_idx < len(w) and end_idx <= len(w):
                w[start_idx:end_idx] = 1.0

                if Nramp > 0 and len(up_ramp) > 0:
                    ramp_end = min(start_idx + Nramp, end_idx)
                    w[start_idx:ramp_end] = up_ramp[:ramp_end - start_idx]

                    ramp_start = max(end_idx - Nramp, start_idx)
                    if ramp_start < end_idx:
                        w[ramp_start:end_idx] = down_ramp[:end_idx - ramp_start]
        else:
            T1_new = Nend - T2
            if T1_new > 0:
                start_idx = gci_ins[i] + T2
                end_idx = gci_ins[i] + T2 + T1_new

                if start_idx < len(w) and end_idx <= len(w):
                    w[start_idx:end_idx] = 1.0

                    if Nramp > 0 and len(up_ramp) > 0:
                        ramp_end = min(start_idx + Nramp, end_idx)
                        w[start_idx:ramp_end] = up_ramp[:ramp_end - start_idx]

    # Return only the first N samples
    return w[:N]

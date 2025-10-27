"""
Fully optimized formant tracking with vectorization + Numba JIT.

This combines:
- Phase 1: Vectorized TVLP matrix construction (11-16x)
- Phase 2A: Numba JIT for GCI detection (50-300x on loops)

Expected performance: ~2.0s (1.8x real-time)
"""

import numpy as np
from scipy import signal
from scipy.signal import medfilt

from .tvlp_optimized import tvlp_l1_opt, tvlp_l2_opt, tvwlp_l1_opt, tvwlp_l2_opt
from .formants import tvlptoformants_akitofi
from .gloat.pitch import srh_pitch_tracking
from .gloat.gci_fast import sedreams_gci_detection_fast


def ftrack_tvwlp_fast(s, fs, lptype='tvwlp_l2', nwin=None, nshift=None, p=8, q=3,
                      npeaks=3, preemp=0.97, fint=80, plot_flag=False, use_numba=True):
    """
    Fully optimized formant tracking (Vectorization + Numba JIT).

    This is the fastest single-core Python implementation, combining:
    - Vectorized TVLP matrix construction (Phase 1)
    - Numba JIT compilation of hot loops (Phase 2A)

    Expected speedup: 1.9-2.0x vs original
    Expected runtime: ~2.0s for 3.6s audio (1.8x real-time)

    Parameters
    ----------
    s : ndarray
        Speech signal
    fs : int
        Sampling rate (Hz)
    lptype : str
        LP method: 'tvwlp_l2', 'tvwlp_l1', 'tvlp_l2', 'tvlp_l1'
    nwin : int, optional
        Window size (default: 1600 @ 8kHz)
    nshift : int, optional
        Window shift (default: 1600 @ 8kHz)
    p : int
        LP order (default: 8)
    q : int
        Polynomial order (default: 3)
    npeaks : int
        Number of formants (default: 3)
    preemp : float
        Pre-emphasis factor (default: 0.97)
    fint : int
        Output interval in samples (default: 80 = 10ms)
    plot_flag : bool
        Plot results (default: False)
    use_numba : bool
        Use Numba JIT if available (default: True)

    Returns
    -------
    Fi : ndarray
        Formant tracks [npeaks x time] in Hz
    Ak : ndarray
        Time-varying LPC coefficients [p+1 x time]
    """
    s = np.asarray(s, dtype=np.float64)
    fs_ref = 8000
    n1ms = int(np.floor(fs_ref / 1000))

    # Set defaults
    if nwin is None:
        nwin = 200 * n1ms
    if nshift is None:
        nshift = 200 * n1ms

    # Resample to 8kHz
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
        wmin = 0.00001
        DQ = 0.7
        PQ = 0.05
        Nramp = 3

        f0_min = 60
        f0_max = 600

        # Pitch tracking
        f0, vuv, srh_val = srh_pitch_tracking(s, fs, f0_min, f0_max)

        # Compute mean F0
        f0_tmp = f0 * vuv
        pos = f0_tmp != 0
        if np.sum(pos) > 0:
            f0_mean = np.mean(f0_tmp[pos])
        else:
            f0_mean = 150

        # GCI detection with Numba optimization
        gc, mbs = sedreams_gci_detection_fast(s, fs, f0_mean, use_numba=use_numba)

        # Compute QCP weight function
        w = _qcp_wt(s, p, DQ, PQ, wmin, Nramp, gc, fs)
    else:
        w = np.ones(len(s), dtype=np.float64)

    # Initialize output
    Fi_list = []
    Ak_list = []
    j = 0
    first_frame = True

    # Sliding window analysis
    while j <= Ns - nwin:
        if j <= Ns - nwin - nshift:
            x = shat[j:j + nwin]
            wj = w[j:j + nwin]
        else:
            x = shat[j:]
            wj = w[j:]

        # Use vectorized TVLP functions (Phase 1 optimization)
        if lptype_lower == 'tvlp_l2':
            aki = tvlp_l2_opt(x, p, q)
        elif lptype_lower == 'tvwlp_l2':
            aki = tvwlp_l2_opt(x, p, q, wj)
        elif lptype_lower == 'tvlp_l1':
            aki = tvlp_l1_opt(x, p, q)
        elif lptype_lower == 'tvwlp_l1':
            aki = tvwlp_l1_opt(x, p, q, wj)
        else:
            raise ValueError(f"Unknown lptype: {lptype}")

        Nx = len(x)
        fi, ak = tvlptoformants_akitofi(aki, Nx, npeaks, fs)

        # Handle frame boundaries
        half_overlap = int(np.floor((nwin - nshift) / 2))

        if first_frame:
            Fi_list.append(fi[:half_overlap, :])
            Ak_list.append(ak[:, :half_overlap])
            first_frame = False

        if j <= Ns - nwin - nshift:
            Fi_list.append(fi[half_overlap:half_overlap + nshift, :])
            Ak_list.append(ak[:, half_overlap:half_overlap + nshift])
        else:
            Fi_list.append(fi[nshift + half_overlap:, :])
            Ak_list.append(ak[:, nshift + half_overlap:])

        j += nshift

    # Concatenate frames
    if len(Fi_list) > 0:
        Fi = np.vstack(Fi_list)
        Ak = np.hstack(Ak_list)
    else:
        Fi = np.zeros((Ns, npeaks))
        Ak = np.zeros((p + 1, Ns))

    # Pad to signal length
    if Fi.shape[0] < Ns:
        pad_rows = Ns - Fi.shape[0]
        Fi = np.vstack([Fi, np.zeros((pad_rows, npeaks))])
        Ak = np.hstack([Ak, np.zeros((p + 1, pad_rows))])

    # Downsample formants
    Fi = Fi[fint - 1::fint, :]
    Fi = Fi.T

    # Median filtering
    Fi_filtered = np.zeros_like(Fi)
    for i in range(npeaks):
        Fi_filtered[i, :] = medfilt(Fi[i, :], kernel_size=5)
    Fi = Fi_filtered

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

        t_formants = np.arange(fint, Ns + 1, fint) / fs
        for i in range(npeaks):
            plt.plot(t_formants, Fi[i, :] / 1000, '.k', markersize=2)
        plt.title(f'Formant Tracks ({lptype}) - FAST (Vectorized + Numba)')
        plt.show()

    return Fi, Ak


def _qcp_wt(x, p, DQ, PQ, d, Nramp, gci_ins, fs):
    """QCP weight function (same as original)."""
    N = len(x)

    if Nramp > 0:
        up_ramp = np.linspace(d, 1, 2 + Nramp)[1:-1]
        down_ramp = up_ramp[::-1]
    else:
        up_ramp = np.array([])
        down_ramp = np.array([])

    if DQ + PQ > 1:
        DQ = 1 - PQ

    w = d * np.ones(N + p, dtype=np.float64)

    if len(gci_ins) < 2:
        return w[:N]

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

    # Handle last period
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

    return w[:N]

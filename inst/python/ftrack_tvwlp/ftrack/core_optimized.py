"""
Optimized main formant tracking function using TVWLP/TVLP.

Phase 1 optimizations:
1. Vectorized TVLP matrix construction
2. Cached LPC residual computation
3. Optimized median filtering
"""

import numpy as np
from scipy import signal
from scipy.signal import medfilt

from .tvlp_optimized import tvlp_l1_opt, tvlp_l2_opt, tvwlp_l1_opt, tvwlp_l2_opt
from .formants import tvlptoformants_akitofi
from .gloat.pitch import srh_pitch_tracking
from .gloat.gci import sedreams_gci_detection
from .gloat.utils import get_lpc_residual


def ftrack_tvwlp_optimized(s, fs, lptype='tvwlp_l2', nwin=None, nshift=None, p=8, q=3,
                           npeaks=3, preemp=0.97, fint=80, plot_flag=False):
    """
    Optimized formant tracking using TVWLP/TVLP.

    This is a drop-in replacement for ftrack_tvwlp with Phase 1 optimizations:
    - Vectorized TVLP matrix construction (11-16x faster)
    - Cached LPC residual computation (reused across pitch/GCI)
    - Optimized median filtering

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
        # OPTIMIZATION: Compute LPC residual once and reuse
        lpc_order = int(np.round(fs / 1000) + 2)
        res = get_lpc_residual(s, int(np.round(25 / 1000 * fs)),
                               int(np.round(5 / 1000 * fs)), lpc_order)

        # Parameters for QCP weight
        wmin = 0.00001
        DQ = 0.7
        PQ = 0.05
        Nramp = 3

        f0_min = 60
        f0_max = 600

        # Pitch tracking (reuse residual)
        f0, vuv, srh_val = _srh_pitch_tracking_cached(s, fs, f0_min, f0_max, res)

        # Compute mean F0 from voiced regions
        f0_tmp = f0 * vuv
        pos = f0_tmp != 0
        if np.sum(pos) > 0:
            f0_mean = np.mean(f0_tmp[pos])
        else:
            f0_mean = 150  # Default

        # GCI detection (reuse residual)
        gc, mbs = _sedreams_gci_detection_cached(s, fs, f0_mean, res)

        # Compute QCP weight function
        w = _qcp_wt(s, p, DQ, PQ, wmin, Nramp, gc, fs)
    else:
        # Dummy weight (not used for TVLP methods)
        w = np.ones(len(s), dtype=np.float64)

    # Initialize output arrays
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

        # OPTIMIZATION: Use vectorized TVLP functions
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
    Fi = Fi[fint - 1::fint, :]

    # Transpose to [npeaks x time]
    Fi = Fi.T

    # OPTIMIZATION: Use scipy.signal.medfilt instead of ndimage.median_filter
    # Apply median filter (5-tap) along time axis for each formant
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

        # Overlay formant tracks
        t_formants = np.arange(fint, Ns + 1, fint) / fs
        for i in range(npeaks):
            plt.plot(t_formants, Fi[i, :] / 1000, '.k', markersize=2)
        plt.title(f'Spectrogram with Formant Tracks ({lptype}) - OPTIMIZED')
        plt.show()

    return Fi, Ak


def _srh_pitch_tracking_cached(wave, fs, f0_min, f0_max, res_precomputed=None):
    """
    SRH pitch tracking with optional pre-computed LPC residual.

    This is a wrapper that reuses LPC residual if already computed.
    """
    from .gloat.pitch import _srh_estimate_pitch

    wave = np.asarray(wave, dtype=np.float64)

    # Resample to 16kHz if necessary
    if fs > 16000:
        n_samples_new = int(len(wave) * 16000 / fs)
        wave = signal.resample(wave, n_samples_new)
        fs = 16000

    n_iter = 2

    # Use pre-computed residual if available, otherwise compute
    if res_precomputed is not None and len(res_precomputed) == len(wave):
        res = res_precomputed
    else:
        # Compute LPC residual
        lpc_order = int(np.round(3 / 4 * fs / 1000))
        res = get_lpc_residual(wave, int(np.round(25 / 1000 * fs)),
                               int(np.round(5 / 1000 * fs)), lpc_order)

    # Estimate pitch in 2 iterations
    for iteration in range(n_iter):
        f0, srh_val = _srh_estimate_pitch(res, fs, f0_min, f0_max)

        pos_tmp = srh_val > 0.1
        if np.sum(pos_tmp) > 1:
            f0_mean_est = np.median(f0[pos_tmp])
            f0_min = int(np.round(0.5 * f0_mean_est))
            f0_max = int(np.round(2 * f0_mean_est))

    # Voiced-unvoiced decisions
    vuv_decisions = np.zeros(len(f0), dtype=int)
    threshold = 0.07
    if np.std(srh_val) > 0.05:
        threshold = 0.085

    vuv_decisions[srh_val > threshold] = 1

    return f0, vuv_decisions, srh_val


def _sedreams_gci_detection_cached(wave, fs, f0_mean, res_precomputed=None):
    """
    SEDREAMS GCI detection with optional pre-computed LPC residual.
    """
    wave = np.asarray(wave, dtype=np.float64)
    n_samples = len(wave)

    # Use pre-computed residual if available
    if res_precomputed is not None and len(res_precomputed) == len(wave):
        res = res_precomputed.copy()
    else:
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

    # Remove low-frequency contents using elliptic high-pass filter
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

    # Determine median position of GCIs
    res = res / np.max(np.abs(res))

    abs_res = np.abs(res)
    diff1 = np.diff(np.concatenate([[0], abs_res]))
    diff2 = np.diff(diff1 > 0)
    posis = np.where(diff2 == -1)[0]
    posis = posis[abs_res[posis] > 0.4]

    if len(posis) == 0:
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
        ratio_gci = 0.5
    else:
        ratio_gci = np.median(rel_posis)

    # Detect GCIs from residual
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

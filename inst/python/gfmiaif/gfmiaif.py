"""
GFM-IAIF: Glottal Flow Model-based Iterative Adaptive Inverse Filtering

Optimized implementation for source-filter separation of speech signals.

References
----------
[1] O. Perrotin and I. V. McLoughlin (2019)
    "A spectral glottal flow model for source-filter separation of speech",
    IEEE ICASSP 2019, pp. 7160-7164.

Copyright (c) 2019 Univ. Grenoble Alpes, CNRS, Grenoble INP, GIPSA-lab
Python implementation (c) 2025

License: GNU Lesser General Public License v3.0 or later
"""

import numpy as np
from scipy import signal
from .lpc import lpc_fast


def gfmiaif_fast(s_gvl, nv=48, ng=3, d=0.99, win=None):
    """
    Fast GFM-IAIF implementation.

    Estimates linear prediction coefficients for vocal tract, glottis,
    and lip radiation filters from a speech signal frame.

    Parameters
    ----------
    s_gvl : array_like
        Speech signal frame, shape (N,)
    nv : int, optional
        Order of LP analysis for vocal tract (default: 48)
    ng : int, optional
        Order of LP analysis for glottal source (default: 3)
        Note: ng=3 is highly recommended
    d : float, optional
        Leaky integration coefficient for lip radiation (default: 0.99)
    win : array_like, optional
        Window used before LPC analysis, shape (N,)
        If None, a Hanning window is used (default: None)

    Returns
    -------
    av : ndarray
        LP coefficients of vocal tract, shape (nv+1,)
    ag : ndarray
        LP coefficients of glottis, shape (ng+1,)
    al : ndarray
        LP coefficients of lip radiation, shape (2,)
    """
    # Convert input to numpy array
    s_gvl = np.asarray(s_gvl, dtype=np.float64)
    if s_gvl.ndim > 1:
        s_gvl = s_gvl.flatten()

    N = len(s_gvl)

    # Set default window
    if win is None:
        win = np.hanning(N)
    else:
        win = np.asarray(win, dtype=np.float64).flatten()

    # ----- Addition of pre-frame --------------------------------------------
    Lpf = nv + 1  # Pre-frame length

    # Create ramp
    ramp = np.linspace(-s_gvl[0], s_gvl[0], Lpf)
    x_gvl = np.concatenate([ramp, s_gvl])

    # Pre-compute slice for removing pre-frame
    idx_pf = slice(Lpf, len(x_gvl))

    # ----- Cancel lip radiation contribution --------------------------------
    al = np.array([1.0, -d])

    # Integration using lfilter
    s_gv = signal.lfilter([1.0], al, s_gvl)
    x_gv = signal.lfilter([1.0], al, x_gvl)

    # ----- Gross glottis estimation -----------------------------------------
    # First 1st order LPC estimation
    ag1, _ = lpc_fast(s_gv * win, 1)

    # Iterate ng-1 times
    for i in range(ng - 1):
        x_v1x = signal.lfilter(ag1, [1.0], x_gv)
        s_v1x = x_v1x[idx_pf]
        ag1x, _ = lpc_fast(s_v1x * win, 1)
        ag1 = np.convolve(ag1, ag1x)

    # ----- Gross vocal tract estimation -------------------------------------
    x_v1 = signal.lfilter(ag1, [1.0], x_gv)
    s_v1 = x_v1[idx_pf]
    av1, _ = lpc_fast(s_v1 * win, nv)

    # ----- Fine glottis estimation ------------------------------------------
    x_g1 = signal.lfilter(av1, [1.0], x_gv)
    s_g1 = x_g1[idx_pf]
    ag, _ = lpc_fast(s_g1 * win, ng)

    # ----- Fine vocal tract estimation --------------------------------------
    x_v = signal.lfilter(ag, [1.0], x_gv)
    s_v = x_v[idx_pf]
    av, _ = lpc_fast(s_v * win, nv)

    return av, ag, al


def gfmiaif_frame_based(audio, sample_rate, window_shift=0.010, window_size=0.032,
                        nv=48, ng=3, d=0.99, window_type='hann'):
    """
    Process audio signal frame-by-frame with GFM-IAIF.

    This function is optimized for integration with R via reticulate,
    processing an entire audio signal and returning time-aligned tracks.

    Parameters
    ----------
    audio : array_like
        Audio signal (mono, float32 or float64)
    sample_rate : int
        Sampling rate in Hz
    window_shift : float, optional
        Frame shift in seconds (default: 0.010 = 10ms)
    window_size : float, optional
        Frame size in seconds (default: 0.032 = 32ms)
    nv : int, optional
        Vocal tract LPC order (default: 48)
    ng : int, optional
        Glottis LPC order (default: 3, highly recommended)
    d : float, optional
        Leaky integration coefficient (default: 0.99)
    window_type : str, optional
        Window function type: 'hann', 'hamming', 'blackman' (default: 'hann')

    Returns
    -------
    dict
        Dictionary with keys:
        - 'av': Vocal tract coefficients, shape (n_frames, nv+1)
        - 'ag': Glottis coefficients, shape (n_frames, ng+1)
        - 'al': Lip radiation coefficients, shape (n_frames, 2)
        - 'timestamps': Frame center times in seconds, shape (n_frames,)
        - 'sample_rate_hz': Frame rate in Hz (1/window_shift)
        - 'n_frames': Number of frames processed
    """
    # Convert to numpy array
    audio = np.asarray(audio, dtype=np.float64).flatten()

    # Calculate frame parameters
    frame_length = int(window_size * sample_rate)
    frame_shift = int(window_shift * sample_rate)

    # Create window function
    if window_type.lower() == 'hann':
        window = np.hanning(frame_length)
    elif window_type.lower() == 'hamming':
        window = np.hamming(frame_length)
    elif window_type.lower() == 'blackman':
        window = np.blackman(frame_length)
    else:
        raise ValueError(f"Unknown window type: {window_type}")

    # Calculate number of frames
    n_frames = (len(audio) - frame_length) // frame_shift + 1

    if n_frames < 1:
        raise ValueError(f"Audio too short: {len(audio)} samples, need at least {frame_length}")

    # Pre-allocate output arrays
    av_frames = np.zeros((n_frames, nv + 1), dtype=np.float64)
    ag_frames = np.zeros((n_frames, ng + 1), dtype=np.float64)
    al_frames = np.zeros((n_frames, 2), dtype=np.float64)
    timestamps = np.zeros(n_frames, dtype=np.float64)

    # Process frames
    for i in range(n_frames):
        # Extract frame
        start_idx = i * frame_shift
        end_idx = start_idx + frame_length
        frame = audio[start_idx:end_idx]

        # Process with GFM-IAIF
        av, ag, al = gfmiaif_fast(frame, nv=nv, ng=ng, d=d, win=window)

        # Store results
        av_frames[i, :] = av
        ag_frames[i, :] = ag
        al_frames[i, :] = al

        # Calculate timestamp (center of frame)
        timestamps[i] = (start_idx + frame_length / 2.0) / sample_rate

    # Return as dictionary (easy to convert to R list)
    return {
        'av': av_frames,
        'ag': ag_frames,
        'al': al_frames,
        'timestamps': timestamps,
        'sample_rate_hz': 1.0 / window_shift,
        'n_frames': n_frames,
        'frame_shift_sec': window_shift,
        'frame_size_sec': window_size
    }

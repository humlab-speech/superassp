"""
Main SE_VQ algorithm for GCI detection.

This module implements the SE_VQ (Spectral Envelope - Voice Quality) algorithm,
an improved version of SEDREAMS for detecting glottal closure instants across
various voice qualities.
"""

import numpy as np
from ..general.pitch import srh_pitch_tracking
from ..general.signal_utils import get_lpc_residual
from .get_mbs import get_mbs
from .get_mbs_gci_intervals import get_mbs_gci_intervals
from .rcvd_reson_gci import rcvd_reson_gci
from .search_res_interval_peaks import search_res_interval_peaks
from .reson_dyprog_mat import reson_dyprog_mat
from .gci_creak_postproc import gci_creak_postproc


def se_vq(x, fs, f0=None, vuv=None, creak=None):
    """
    Extract GCIs using the SE_VQ algorithm.

    The SE_VQ algorithm is an adapted version of SEDREAMS optimized for
    non-modal voice qualities. It uses dynamic programming to select optimal
    GCI candidates from the LP-residual signal.

    Parameters
    ----------
    x : array_like
        Speech signal (samples)
    fs : float
        Sampling frequency (Hz)
    f0 : array_like, optional
        Fundamental frequency contour (Hz). If not provided, SRH will be used.
    vuv : array_like, optional
        Voiced/unvoiced decisions (binary). If not provided, SRH will be used.
    creak : array_like, optional
        Creaky voice binary decision for each sample

    Returns
    -------
    gci : ndarray
        Detected glottal closure instants (sample indices)
    rep : ndarray
        Resonator output signal
    res : ndarray
        LP-residual signal
    mbs : ndarray
        Mean-based signal

    References
    ----------
    Kane, J., Gobl, C., (2013) 'Evaluation of glottal closure instant detection
    in a range of voice qualities', Speech Communication 55(2), pp. 295-314.

    Notes
    -----
    Original MATLAB implementation by John Kane @ Phonetics and Speech Lab,
    Trinity College Dublin. SEDREAMS algorithm by Thomas Drugman, University
    of Mons. Python conversion: 2025.
    """
    x = np.asarray(x).flatten()

    # Settings
    F0min = 20
    F0max = 500

    # Get F0 and voicing if not provided
    if f0 is None or vuv is None:
        print("Computing F0 and voicing using SRH algorithm...")
        f0, vuv, _, _ = srh_pitch_tracking(x, fs, F0min, F0max)

    f0 = np.asarray(f0)
    vuv = np.asarray(vuv)

    # Calculate mean F0
    valid_f0 = f0[(f0 > F0min) & (f0 < F0max) & (vuv == 1)]
    if len(valid_f0) > 0:
        F0mean = np.median(valid_f0)
    else:
        F0mean = 100  # Default

    # Update F0max
    voiced_f0 = f0[vuv == 1]
    if len(voiced_f0) > 0:
        # Apply median filter
        from scipy.signal import medfilt
        F0max = np.max(medfilt(voiced_f0, kernel_size=min(13, len(voiced_f0) if len(voiced_f0) % 2 == 1 else len(voiced_f0)-1)))
    else:
        F0max = 500

    # Handle low F0 (likely creak)
    if F0mean < 70:
        print('Utterance likely to contain creak')
        F0mean = 80

    T0mean = fs / F0mean  # Rough period length

    # Algorithm parameters
    win_len = 25  # window length in ms
    win_shift = 5  # window shift in ms
    lpc_ord = int((fs / 1000) + 2)  # LPC order
    Ncand = 5  # Number of candidate GCI peaks

    trans_wgt = 1  # Transition cost weight
    rel_amp_wgt = 0.3  # Local cost weight

    rep_num = 2
    remove_thresh = 0.4  # Threshold for removing false GCIs
    search_reg = int(1.3 / 1000 * fs)

    # Calculate LP-residual
    print("Computing LP-residual...")
    L = int(win_len / 1000 * fs)
    shift = int(win_shift / 1000 * fs)
    res, _ = get_lpc_residual(x, L, shift, lpc_ord)

    # Get resonator output
    print("Computing resonator output...")
    rep = rcvd_reson_gci(res, fs, F0mean)

    # Extract mean-based signal
    print("Computing mean-based signal...")
    mbs = get_mbs(x, fs, T0mean)

    # Define search intervals
    print("Finding GCI search intervals...")
    interval = get_mbs_gci_intervals(mbs, fs, T0mean, F0max)

    # Find residual peaks
    print(f"Searching for {Ncand} candidate peaks per interval...")
    gci_n, gci_rel_amp = search_res_interval_peaks(res, interval, Ncand)

    # Dynamic programming
    print("Running dynamic programming to select optimal GCI path...")
    gci = reson_dyprog_mat(gci_rel_amp.T, gci_n.T, F0mean, x, fs,
                           trans_wgt, rel_amp_wgt)

    # Post-processing for creaky voice
    if creak is not None and len(creak) == len(x):
        print('Post-processing in detected creaky voice regions...')
        gci = gci_creak_postproc(gci, creak, search_reg, rep,
                                 remove_thresh, rep_num)

    return gci, rep, res, mbs


def se_vq_var_f0(x, fs, f0=None, vuv=None, creak=None):
    """
    SE_VQ algorithm for signals with highly varying F0.

    This variant is optimized for speech with rapid F0 variations.

    Parameters
    ----------
    x : array_like
        Speech signal (samples)
    fs : float
        Sampling frequency (Hz)
    f0 : array_like, optional
        Fundamental frequency contour (Hz)
    vuv : array_like, optional
        Voiced/unvoiced decisions (binary)
    creak : array_like, optional
        Creaky voice binary decision

    Returns
    -------
    gci : ndarray
        Detected glottal closure instants (sample indices)
    rep : ndarray
        Resonator output signal
    res : ndarray
        LP-residual signal
    mbs : ndarray
        Mean-based signal

    Notes
    -----
    This function currently uses the same implementation as se_vq.
    Future versions may include specific adaptations for varying F0.
    """
    # For now, use the same implementation
    # Future: Add specific handling for varying F0
    return se_vq(x, fs, f0, vuv, creak)

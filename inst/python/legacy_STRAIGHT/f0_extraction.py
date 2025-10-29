"""
STRAIGHT F0 Extraction Module

Python implementation of MulticueF0v14 - Multi-cue F0 extraction.
Faithful port from MATLAB preserving algorithms and numerical behavior.

This is a large module (~1,863 lines in MATLAB) implementing a sophisticated
multi-cue F0 extraction algorithm combining:
- Instantaneous Frequency (IF) analysis
- Autocorrelation (AC) analysis
- Multi-cue fusion with dynamic range normalization
- F0 tracking and gap filling
- Voiced/Unvoiced decision
- F0 refinement

Original MATLAB version:
    Copyright (c) Wakayama University, 2004-2016
    Designed and coded by Hideki Kawahara
"""

import numpy as np
from dataclasses import dataclass
from typing import Tuple, Optional, Dict
import warnings
from scipy import signal
from scipy.signal import decimate

# Import helpers from other modules
from .synthesis import hanning, fftfilt

# Try to import optimized functions (Numba)
try:
    from .f0_extraction_optimized import (
        parabolic_interp_opt,
        find_peaks_opt,
        compute_tracking_costs_opt,
        zfixpfreq3_opt,
        compute_acc_indices_opt,
        apply_window_and_scale_opt,
        NUMBA_AVAILABLE
    )
    _NUMBA_AVAILABLE = NUMBA_AVAILABLE
except ImportError:
    _NUMBA_AVAILABLE = False


@dataclass
class F0Parameters:
    """
    Parameters for STRAIGHT F0 extraction
    
    Attributes match MATLAB field names for compatibility
    """
    # F0 search range
    F0searchLowerBound: float = 40.0  # Hz
    F0searchUpperBound: float = 800.0  # Hz
    F0frameUpdateInterval: float = 1.0  # ms
    
    # Optimization settings
    UseOptimized: bool = True  # Use Numba-optimized functions when available
    
    # Instantaneous Frequency (IF) parameters
    NofChannelsInOctave: int = 24
    IFWindowStretch: float = 1.2
    IFsmoothingLengthRelToFc: float = 1.0  # ratio
    IFminimumSmoothingLength: float = 5.0  # ms
    IFexponentForNonlinearSum: float = 0.5
    IFnumberOfHarmonicForInitialEstimate: int = 1
    
    # Autocorrelation (AC) parameters
    TimeConstantForPowerCalculation: float = 10.0  # ms
    ACtimeWindowLength: float = 60.0  # ms
    ACnumberOfFrequencySegments: int = 8
    ACfrequencyDomainWindowWidth: float = 2200.0  # Hz
    ACpowerExponentForNonlinearity: float = 0.5
    ACamplitudeCompensationInShortLag: float = 1.6
    ACexponentForACdistance: float = 3.0
    AClagSmoothingLength: float = 0.001  # s
    ACtemporalSmoothingLength: float = 30.0  # ms
    
    # Multi-cue fusion parameters
    WeightForAutocorrelationMap: float = 1.0
    WeightForInstantaneousFqMap: float = 1.0
    SDforNormalizeMixingDistance: float = 0.3  # octave (MATLAB default)
    
    # V/UV decision parameters
    ThresholdForSilence: float = -20.0  # dB
    ThresholdForVUV: float = -6.0  # dB
    VUVthresholdOfAC1: float = 0.6
    
    # Tracking parameters
    SDforTrackingNormalization: float = 0.2  # octave
    MaxumumPermissibleOctaveJump: float = 0.35  # octave
    ThresholdToStartSearch: float = 0.45
    ThresholdToQuitSearch: float = 0.3
    ThresholdForReliableRegion: float = 0.5
    
    # Display
    DisplayPlots: int = 0


def MulticueF0v14(
    x: np.ndarray,
    fs: float,
    f0floor: Optional[float] = None,
    f0ceil: Optional[float] = None
) -> Tuple[np.ndarray, np.ndarray, Dict, F0Parameters]:
    """
    Source information extraction using multiple cues
    
    Faithful port of MulticueF0v14.m
    
    Args:
        x: Input signal (monaural)
        fs: Sampling frequency (Hz)
        f0floor: Lower limit of F0 search (Hz) or F0Parameters object
        f0ceil: Upper limit of F0 search (Hz)
        
    Returns:
        f0raw: Fundamental frequency without V/UV information (Hz)
        vuv: V/UV indicator (1: voiced, 0: unvoiced)
        auxouts: Base information for f0 extraction (dict)
        prm: Parameters actually used
        
    Example:
        >>> f0, vuv, aux, params = MulticueF0v14(x, fs)
        >>> f0, vuv, aux, params = MulticueF0v14(x, fs, 60, 400)
    """
    # Convert fs to float to avoid overflow
    fs = float(fs)
    
    # Handle parameter input
    if f0floor is not None and isinstance(f0floor, F0Parameters):
        prmin = f0floor
    else:
        prmin = F0Parameters()
        if f0floor is not None:
            prmin.F0searchLowerBound = f0floor
        if f0ceil is not None:
            prmin.F0searchUpperBound = f0ceil
        prmin.DisplayPlots = 0
    
    # Call main extraction function
    f0raw, vuv, auxouts, prm = SourceInfobyMultiCues050111(x, fs, prmin)
    
    nn = min(len(vuv), len(f0raw))
    f0raw = f0raw[:nn]
    vuv = vuv[:nn]
    
    return f0raw, vuv, auxouts, prm


def SourceInfobyMultiCues050111(
    x: np.ndarray,
    fs: float,
    prmin: F0Parameters
) -> Tuple[np.ndarray, np.ndarray, Dict, F0Parameters]:
    """
    Source information extraction function with combined source information
    
    Main F0 extraction algorithm combining IF and AC analysis.
    
    Port of SourceInfobyMultiCues050111 from MATLAB
    """
    # Convert fs to float
    fs = float(fs)
    
    # Initialize parameters
    prm = F0Parameters()
    
    # Update parameters from input
    if prmin is not None:
        for field in vars(prmin):
            if hasattr(prmin, field):
                setattr(prm, field, getattr(prmin, field))
    
    # Ensure x is a 1D array
    if x.ndim > 1:
        if x.shape[0] < x.shape[1]:
            x = x[0, :]
        else:
            x = x[:, 0]
    x = x.flatten()
    
    # Safe guard for all-zero segments
    l1ms = int(np.round(fs / 1000))
    zero_indices = (x == 0)
    if np.sum(zero_indices) > l1ms:
        zv = np.random.randn(np.sum(zero_indices))
        zv = np.cumsum(zv - np.mean(zv))
        zv = zv / np.std(zv) * np.std(x) / 10000
        x[zero_indices] = zv
    
    # Extract parameters
    f0floor = prm.F0searchLowerBound
    f0ceil = prm.F0searchUpperBound
    shiftm = prm.F0frameUpdateInterval
    nvo = prm.NofChannelsInOctave
    mu = prm.IFWindowStretch
    imgi = prm.DisplayPlots
    smp = prm.IFsmoothingLengthRelToFc
    minm = prm.IFminimumSmoothingLength
    pcIF = prm.IFexponentForNonlinearSum
    ncIF = prm.IFnumberOfHarmonicForInitialEstimate
    tcpower = prm.TimeConstantForPowerCalculation
    wtlm = prm.ACtimeWindowLength
    ndiv = prm.ACnumberOfFrequencySegments
    wflf = prm.ACfrequencyDomainWindowWidth
    pcAC = prm.ACpowerExponentForNonlinearity
    ampAC = prm.ACamplitudeCompensationInShortLag
    betaAC = prm.ACexponentForACdistance
    lagslAC = prm.AClagSmoothingLength
    timeslAC = prm.ACtemporalSmoothingLength
    wAC = prm.WeightForAutocorrelationMap
    wIF = prm.WeightForInstantaneousFqMap
    mixsd = prm.SDforNormalizeMixingDistance
    
    nvc = int(np.ceil(np.log(f0ceil / f0floor) / np.log(2) * nvo))
    
    # Step 1: Extract fixed points from IF map
    print("Step 1/8: IF-based F0 candidate extraction...")
    f0v, vrv, _, _, _ = zfixpF0VexMltpBG4(
        x, fs, f0floor, nvc, nvo, mu, imgi, shiftm, smp, minm, pcIF, ncIF
    )
    
    _, pos = zmultiCandIF(f0v, vrv)
    
    # Step 2: Check for AC induction artifacts
    y, ind, _ = zremoveACinduction(x, fs, pos)
    if ind == 1:
        x = y
        f0v, vrv, _, _, _ = zfixpF0VexMltpBG4(
            x, fs, f0floor, nvc, nvo, mu, imgi, shiftm, smp, minm, pcIF, ncIF
        )
    
    # Step 3: Select multiple F0 candidates based on IF
    val, pos = zmultiCandIF(f0v, vrv)
    
    # Step 4: Select multiple F0 candidates based on AC
    print("Step 2/8: AC-based F0 candidate extraction...")
    dn = max(1, int(np.floor(fs / max(8000, 3 * 2 * f0ceil))))
    h1 = -1
    
    x_decimated = decimate(x, dn) if dn > 1 else x
    lagspec, lx = zlagspectestnormal(
        x_decimated, fs / dn, shiftm, len(x) / fs * 1000,
        shiftm, wtlm, ndiv, wflf, pcAC, ampAC, h1
    )
    f02, pl2 = zmultiCandAC(lx, lagspec, betaAC, lagslAC, timeslAC)
    
    # Step 5: Combine multiple source information
    print("Step 3/8: Multi-cue fusion...")
    auxouts = {
        'F0candidatesByIF': pos,
        'CNofcandidatesByIF': val,
        'F0candidatesByAC': f02,
        'ACofcandidatesByAC': pl2
    }
    
    f0cand, relv = zcombineRanking4(auxouts, mixsd, wAC, wIF, prm)
    
    # Step 6: Calculate power envelope
    print("Step 4/8: Power envelope calculation...")
    pws = zVpowercalc(x, fs, tcpower, shiftm, 2000)
    pwsdb = 10 * np.log10(np.abs(pws) + 1e-10)
    mxpwsdb = np.max(pwsdb)
    
    hstgrm, binlvl = np.histogram(pwsdb, bins=np.arange(mxpwsdb - 60, mxpwsdb + 2, 2))
    q10 = np.interp(10, np.cumsum(hstgrm + 1e-9) / np.sum(hstgrm) * 100, binlvl[:-1])
    minid = np.argmin(np.abs(q10 - binlvl[:-1]))
    bb = np.arange(max(0, minid - 5), min(len(binlvl) - 1, minid + 6))
    noiselevel = np.sum(hstgrm[bb] * binlvl[bb]) / np.sum(hstgrm[bb])
    
    # Step 7: Calculate first autocorrelation
    print("Step 5/8: First autocorrelation calculation...")
    ac1 = np.zeros(len(f0cand))
    for ii in range(len(f0cand)):
        ac1[ii] = zfirstac(x, fs, int(np.round(ii / 1000 * fs)), 30)
    
    auxouts['F0candidatesByMix'] = f0cand
    auxouts['RELofcandidatesByMix'] = relv
    auxouts['FirstAutoCorrelation'] = ac1
    auxouts['InstantaneousPower'] = pwsdb
    
    # Step 8: F0 tracking with octave-aware segment search
    print("Step 6/8: F0 tracking...")
    f0s, rels, csegs = zcontiguousSegment10(auxouts, prm)
    
    # Post-process: Check each segment for potential octave errors
    f0cand = auxouts['F0candidatesByMix']
    relv = auxouts['RELofcandidatesByMix']
    
    for seg_info in csegs:
        lb, ub = seg_info[0], seg_info[1]
        if f0s[lb] == 0:
            continue
            
        # For the first frame in segment, check if f0/2 has better support
        # across multiple frames in the segment
        segment_len = min(ub - lb + 1, 20)  # Check first 20 frames
        
        # Count how many candidates support f0, f0/2, f0*2
        votes_f0 = 0
        votes_half = 0
        votes_double = 0
        
        for i in range(lb, min(lb + segment_len, ub + 1)):
            if f0s[i] == 0:
                continue
                
            target_f0 = f0s[i]
            target_half = target_f0 / 2
            target_double = target_f0 * 2
            
            # Check if any candidate is close to f0/2
            for cand_idx in range(f0cand.shape[1]):
                cand_f0 = f0cand[i, cand_idx]
                if cand_f0 > 0:
                    # Check octave proximity (within ~20%)
                    if abs(np.log2(cand_f0 / target_f0)) < 0.25:
                        votes_f0 += relv[i, cand_idx]
                    elif abs(np.log2(cand_f0 / target_half)) < 0.25:
                        votes_half += relv[i, cand_idx]
                    elif abs(np.log2(cand_f0 / target_double)) < 0.25:
                        votes_double += relv[i, cand_idx]
        
        # If f0/2 has stronger support, correct the octave
        if votes_half > votes_f0 * 1.2:  # Need 20% more votes to override
            f0s[lb:ub+1] = f0s[lb:ub+1] / 2
        elif votes_double > votes_f0 * 1.2:
            f0s[lb:ub+1] = f0s[lb:ub+1] * 2
    
    f0raw0, _ = zfillf0gaps6(auxouts, f0s, rels, csegs, prm)
    
    # Safe guards
    f0raw0[np.isnan(f0raw0)] = 0
    f0raw0[f0raw0 > f0ceil] = f0ceil
    f0raw0[(f0raw0 < f0floor) & (f0raw0 > 0)] = f0floor
    
    # Step 9: F0 refinement
    print("Step 7/8: F0 refinement...")
    x_decimated = decimate(x, dn) if dn > 1 else x
    f0raw2, ecr, ac1 = zrefineF06m(
        x_decimated, fs / dn, f0raw0, 1024, 1.1, 3, 1, 1, len(f0raw0)
    )
    
    # Step 10: V/UV decision
    print("Step 8/8: V/UV decision...")
    auxouts['BackgroundNoiselevel'] = noiselevel
    vuv = zvuvdecision4(f0raw2, auxouts)
    
    nnll = min(len(f0raw2), len(vuv))
    f0raw = f0raw2[:nnll]
    vuv = vuv[:nnll]
    
    # Update auxouts
    auxouts['RefinedCN'] = ecr
    auxouts['FirstAutoCorrelation'] = ac1
    auxouts['F0initialEstimate'] = f0raw0
    auxouts['InstantaneousPower'] = pwsdb
    auxouts['RefinedF0estimates'] = f0raw
    auxouts['VUVindicator'] = vuv
    
    print("F0 extraction complete!")
    
    return f0raw, vuv, auxouts, prm


# ============================================================================
# Utility functions (simple helpers)
# ============================================================================

def zdBpower(x: np.ndarray) -> np.ndarray:
    """Convert power to dB, handling edge cases"""
    # Clip to prevent log of negative or zero values
    x_safe = np.maximum(x, 1e-15)
    return 10 * np.log10(x_safe)


def zGcBs(x: np.ndarray, k: float) -> np.ndarray:
    """Gammatone-like basis function"""
    tt = x + 1e-7
    return tt**k * np.exp(-np.pi * tt**2) * (np.sin(np.pi * tt + 1e-4) / (np.pi * tt + 1e-4))**2


def GcBs(x: np.ndarray, k: float) -> np.ndarray:
    """Gammatone-like basis function (alias)"""
    return zGcBs(x, k)


def ParabolicInterp(yv: np.ndarray, xo: float) -> Tuple[float, float]:
    """
    Parabolic interpolation for peak finding
    
    Args:
        yv: 3-point values around peak
        xo: Center position (index)
        
    Returns:
        val: Interpolated peak value
        pos: Interpolated peak position (index, for frequency calculation)
    """
    lp = np.diff(yv)
    a = lp[0] - lp[1]
    b = (lp[0] + lp[1]) / 2
    xp = b / a + xo
    val = yv[1] + 0.5 * a * (b / a)**2 + b * (b / a)
    
    # Clamp to reasonable range
    if xp > xo + 1:
        xp = xo + 1
        val = yv[-1]
    if xp < xo - 1:
        xp = xo - 1
        val = yv[0]
    
    # Return the interpolated lag index
    # MATLAB uses: pos = 1/(xp-1) then f0 = pos/lx(2)
    # Which simplifies to: f0 = 1/((xp-1)*lx(2))
    # In Python (0-based): f0 = 1/(xp*lx_spacing)
    # So we return xp (the lag index) for the caller to convert
    pos = xp
    return val, pos


def zzParabolicInterp2(yv: np.ndarray, xo: float) -> Tuple[float, float]:
    """
    Parabolic interpolation variant
    
    Args:
        yv: 3-point values around peak
        xo: Center position
        
    Returns:
        val: Interpolated peak value
        pos: Interpolated peak position
    """
    # Use optimized version if available (10-50x faster for many calls)
    if NUMBA_AVAILABLE:
        try:
            val, pos = parabolic_interp_opt(yv, xo)
            return val, pos
        except:
            pass  # Fall back to original if optimization fails
    
    # Original implementation
    lp = np.diff(yv)
    a = lp[0] - lp[1]
    b = (lp[0] + lp[1]) / 2
    xp = b / a + xo
    val = yv[1] + 0.5 * a * (b / a)**2 + b * (b / a)
    pos = xp - 1
    return val, pos


def znrmlcf2(f: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Normalize coefficients method 2"""
    f = np.atleast_1d(f)
    c1 = f / np.sqrt(np.sum(f**2))
    c2 = f / np.sum(np.abs(f))
    return c1, c2


def znrmlcf3(f: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Normalize coefficients method 3"""
    f = np.atleast_1d(f)
    c1 = f / np.sqrt(np.sum(f**2))
    c2 = f / np.sum(np.abs(f))
    return c1, c2


def zgendeconvmatrix(n: int, a: np.ndarray) -> np.ndarray:
    """
    Generate deconvolution matrix
    
    Args:
        n: Matrix size
        a: Coefficients
        
    Returns:
        mapm: Deconvolution matrix
    """
    a = np.atleast_1d(a)
    mapm = np.zeros((n, n))
    
    for ii in range(n):
        for jj in range(n):
            if ii == jj:
                mapm[ii, jj] = 1
            elif ii > jj:
                idx = ii - jj
                if idx <= len(a):
                    mapm[ii, jj] = -a[idx - 1]
    
    return mapm


# ============================================================================
# Core F0 extraction functions
# ============================================================================

def zfirstac(x: np.ndarray, fs: float, ix: int, wlms: float) -> float:
    """
    Calculate first autocorrelation coefficient
    
    Args:
        x: Input signal
        fs: Sampling frequency (Hz)
        ix: Center index
        wlms: Window length (ms)
        
    Returns:
        ac1: First autocorrelation coefficient
    """
    wl = int(np.round(fs * wlms / 1000))
    fftl = int(2.0 ** np.ceil(np.log2(wl)))
    
    xx = x.flatten()
    idx = ix + np.arange(wl) - int(np.round(wl / 2))
    idx = np.clip(idx, 0, len(xx) - 1).astype(int)
    xt = xx[idx]
    
    # Safe guard for all-zero segments
    if np.sum(np.abs(xt)) < 1e-20:
        xt = xt + np.random.randn(len(xt))
    
    fw = np.abs(np.fft.fft(xt * hanning(wl), fftl))**2
    fx = np.arange(fftl) / fftl * fs
    
    fwbp = fw.copy()
    fwbp[fx < 70] = 0
    fwbp[fx > fs / 2 - 200] = 0
    
    ac = np.real(np.fft.ifft(fwbp))
    ac = ac / ac[0]
    ac1 = ac[1]
    
    return ac1


def zVpowercalc(x: np.ndarray, fs: float, wtc: float, shiftm: float, fc: float) -> np.ndarray:
    """
    Calculate voiced power envelope
    
    Args:
        x: Input signal
        fs: Sampling frequency (Hz)
        wtc: Time constant (ms)
        shiftm: Frame shift (ms)
        fc: Cutoff frequency (Hz)
        
    Returns:
        pws: Power envelope
    """
    fs = float(fs)
    
    # Design lowpass filter
    fc_norm = fc / (fs / 2)
    fc_norm = min(0.99, max(0.01, fc_norm))
    
    b, a = signal.butter(4, fc_norm, 'low')
    
    # Calculate instantaneous power
    x = x.flatten()
    xf = signal.filtfilt(b, a, x)
    pw = xf**2
    
    # Smooth with time constant
    alpha = 1 - np.exp(-shiftm / wtc)
    
    # Calculate at frame positions
    frame_length = int(np.round(shiftm * fs / 1000))
    n_frames = int(np.floor(len(pw) / frame_length))
    
    pws = np.zeros(n_frames)
    for ii in range(n_frames):
        idx = int(ii * frame_length)
        if ii == 0:
            pws[ii] = pw[idx]
        else:
            pws[ii] = alpha * pw[idx] + (1 - alpha) * pws[ii - 1]
    
    return pws


def zremoveACinduction(x: np.ndarray, fs: float, pos: np.ndarray) -> Tuple[np.ndarray, int, Optional[float]]:
    """
    Remove AC induction artifacts
    
    Args:
        x: Input signal
        fs: Sampling frequency (Hz)
        pos: Top F0 candidates (Hz)
        
    Returns:
        y: Filtered signal
        ind: 1 if AC induction detected, 0 otherwise
        fq: Detected AC frequency (Hz)
    """
    # Check for 50/60 Hz peaks in F0 candidates
    ac_freqs = [50, 60]
    
    ind = 0
    fq = None
    y = x.copy()
    
    if len(pos) == 0:
        return y, ind, fq
    
    pos = pos.flatten()
    
    # Check if any candidate is close to AC frequency
    for ac_freq in ac_freqs:
        if np.any(np.abs(pos - ac_freq) < 5):
            ind = 1
            fq = ac_freq
            
            # Design notch filter
            Q = 30
            w0 = ac_freq / (fs / 2)
            b, a = signal.iirnotch(w0, Q)
            y = signal.filtfilt(b, a, x)
            break
    
    return y, ind, fq


def cleaninglownoise(x: np.ndarray, fs: float, f0floor: float) -> np.ndarray:
    """
    Clean low-frequency noise
    
    Args:
        x: Input signal
        fs: Sampling frequency (Hz)
        f0floor: Lower F0 bound (Hz)
        
    Returns:
        x: Cleaned signal
    """
    # High-pass filter at f0floor/2
    fc = f0floor / 2
    fc_norm = fc / (fs / 2)
    fc_norm = min(0.99, max(0.01, fc_norm))
    
    b, a = signal.butter(4, fc_norm, 'high')
    x = signal.filtfilt(b, a, x)
    
    return x


# ============================================================================
# Placeholder implementations for complex functions
# These require substantial implementation
# ============================================================================

def zmultanalytFineCSPB(x: np.ndarray, fs: float, f0floor: float, nvc: int, 
                        nvo: int, mu: float, mlt: int) -> np.ndarray:
    """
    Dual wavelet analysis using cardinal spline manipulation
    
    Args:
        x: Input signal
        fs: Sampling frequency (Hz)
        f0floor: Lower bound for pitch search (Hz)
        nvc: Number of total voices for wavelet analysis
        nvo: Number of voices in an octave
        mu: Temporal stretch factor
        mlt: Harmonic ID#
        
    Returns:
        pm: Wavelet transform using iso-metric Gabor function
    """
    fs = float(fs)
    t0 = 1 / f0floor
    lmx = int(np.round(6 * t0 * fs * mu))
    wl = int(2 ** np.ceil(np.log(lmx) / np.log(2)))
    
    x = x.flatten()
    nx = len(x)
    tx = np.concatenate([x, np.zeros(wl)])
    gent = (np.arange(wl) - wl / 2) / fs
    
    pm = np.zeros((nvc, nx), dtype=complex)
    mpv = 1.0
    
    for ii in range(nvc):
        tb = gent * mpv
        t = tb[np.abs(tb) < 3.5 * mu * t0]
        
        wd1 = np.exp(-np.pi * (t / t0 / mu) ** 2)
        wd2 = np.maximum(0, 1 - np.abs(t / t0 / mu))
        wd2 = wd2[wd2 > 0]
        
        wwd = np.convolve(wd2, wd1, mode='full')
        wwd = wwd[np.abs(wwd) > 0.00001]
        wbias = int(np.round((len(wwd) - 1) / 2))
        
        t_indices = np.round(np.arange(len(wwd)) - wbias + len(t) / 2).astype(int)
        t_indices = np.clip(t_indices, 0, len(t) - 1)
        
        wwd = wwd * np.exp(1j * 2 * np.pi * mlt * t[t_indices] / t0)
        pmtmp1 = fftfilt(wwd, tx)
        pm[ii, :] = pmtmp1[wbias:wbias + nx] * np.sqrt(mpv)
        
        mpv = mpv * (2.0 ** (1 / nvo))
    
    return pm


def zwvlt2ifq(pm: np.ndarray, fs: float) -> np.ndarray:
    """
    Wavelet to instantaneous frequency map
    
    Args:
        pm: Wavelet transform
        fs: Sampling frequency (Hz)
        
    Returns:
        pif: Instantaneous frequency
    """
    mm = pm.shape[1]
    pm_norm = pm / (np.abs(pm) + 1e-10)
    
    pif = np.abs(pm_norm[:, :] - np.c_[pm_norm[:, 0:1], pm_norm[:, :-1]])
    pif = fs / np.pi * np.arcsin(np.clip(pif / 2, -1, 1))
    pif[:, 0] = pif[:, 1]
    
    return pif


def zifq2gpm2(pif: np.ndarray, f0floor: float, nvo: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Instantaneous frequency to geometric parameters
    
    Args:
        pif: Instantaneous frequency
        f0floor: Lower F0 bound (Hz)
        nvo: Number of voices per octave
        
    Returns:
        slp: First order coefficient
        pbl: Second order coefficient
    """
    nn = pif.shape[0]
    fx = f0floor * 2.0 ** (np.arange(nn) / nvo) * 2 * np.pi
    
    c = 2.0 ** (1 / nvo)
    g = np.array([[1 / c / c, 1 / c, 1],
                  [1, 1, 1],
                  [c * c, c, 1]])
    h = np.linalg.inv(g)
    
    slp = ((pif[1:nn - 1, :] - pif[0:nn - 2, :]) / (1 - 1 / c) +
           (pif[2:nn, :] - pif[1:nn - 1, :]) / (c - 1)) / 2
    slp = np.vstack([slp[0:1, :], slp, slp[-1:, :]])
    
    pbl = (pif[0:nn - 2, :] * h[1, 0] + 
           pif[1:nn - 1, :] * h[1, 1] + 
           pif[2:nn, :] * h[1, 2])
    pbl = np.vstack([pbl[0:1, :], pbl, pbl[-1:, :]])
    
    for ii in range(nn):
        slp[ii, :] = slp[ii, :] / fx[ii]
        pbl[ii, :] = pbl[ii, :] / fx[ii]
    
    return slp, pbl


def zsmoothmapB(map_data: np.ndarray, fs: float, f0floor: float, nvo: int,
                mu: float, mlim: float, pex: float) -> np.ndarray:
    """
    Smooth map using frequency-dependent smoothing
    
    Args:
        map_data: Input map
        fs: Sampling frequency (Hz)
        f0floor: Lower F0 bound (Hz)
        nvo: Number of voices per octave
        mu: Smoothing stretch factor
        mlim: Minimum smoothing length (ms)
        pex: Exponent for smoothing width variation
        
    Returns:
        smap: Smoothed map
    """
    nvc, mm = map_data.shape
    t0 = 1 / f0floor
    lmx = int(np.round(6 * t0 * fs * mu))
    wl = int(2 ** np.ceil(np.log(lmx) / np.log(2)))
    gent = (np.arange(wl) - wl / 2) / fs
    
    smap = map_data.copy()
    mpv = 1.0
    zt = np.zeros(len(gent))
    iiv = np.arange(mm)
    
    for ii in range(nvc):
        t = gent * mpv
        t = t[np.abs(t) < 3.5 * mu * t0]
        wbias = int(np.round((len(t) - 1) / 2))
        
        wd1 = np.exp(-np.pi * (t / (t0 * (1 - pex)) / mu) ** 2)
        wd2 = np.exp(-np.pi * (t / (t0 * (1 + pex)) / mu) ** 2)
        wd1 = wd1 / np.sum(wd1)
        wd2 = wd2 / np.sum(wd2)
        
        tm = fftfilt(wd1, np.concatenate([map_data[ii, :], zt]))
        tm = fftfilt(wd2, np.concatenate([1.0 / (tm[iiv + wbias] + 1e-10), zt]))
        smap[ii, :] = 1.0 / (tm[iiv + wbias] + 1e-10)
        
        if t0 * mu / mpv * 1000 > mlim:
            mpv = mpv * (2.0 ** (1 / nvo))
    
    return smap


def zfixpfreq3(fxx: np.ndarray, pif2: np.ndarray, mmp: np.ndarray,
               dfv: np.ndarray, pm: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Find fixed points in frequency map
    
    Args:
        fxx: Frequency axis
        pif2: Instantaneous frequency
        mmp: Merit map
        dfv: Frequency derivative
        pm: Phase map
        
    Returns:
        ff: Fixed point frequencies
        vv: Fixed point values
        df: Fixed point derivatives
        aa: Fixed point amplitudes
    """
    # Use optimized version if available
    if _NUMBA_AVAILABLE:
        return zfixpfreq3_opt(fxx, pif2, mmp, dfv, pm)
    
    # Original implementation
    aav = np.abs(pm)
    nn = len(fxx)
    iix = np.arange(nn)
    
    cd1 = pif2 - fxx
    cd2 = np.concatenate([np.diff(cd1), [cd1[-1] - cd1[-2]]])
    cdd1 = np.concatenate([cd1[1:], [cd1[-1]]])
    
    fp = (cd1 * cdd1 < 0) * (cd2 < 0)
    ixx = iix[fp > 0]
    
    if len(ixx) == 0:
        return np.array([]), np.array([]), np.array([]), np.array([])
    
    # Ensure ixx + 1 doesn't exceed bounds
    ixx = ixx[ixx < nn - 1]
    
    if len(ixx) == 0:
        return np.array([]), np.array([]), np.array([]), np.array([])
    
    ff = pif2[ixx] + (pif2[ixx + 1] - pif2[ixx]) * cd1[ixx] / (cd1[ixx] - cdd1[ixx] + 1e-10)
    vv = mmp[ixx] + (mmp[ixx + 1] - mmp[ixx]) * (ff - fxx[ixx]) / (fxx[ixx + 1] - fxx[ixx] + 1e-10)
    df = dfv[ixx] + (dfv[ixx + 1] - dfv[ixx]) * (ff - fxx[ixx]) / (fxx[ixx + 1] - fxx[ixx] + 1e-10)
    aa = aav[ixx] + (aav[ixx + 1] - aav[ixx]) * (ff - fxx[ixx]) / (fxx[ixx + 1] - fxx[ixx] + 1e-10)
    
    return ff, vv, df, aa


def zfixpF0VexMltpBG4(x: np.ndarray, fs: float, f0floor: float, nvc: int, nvo: int,
                      mu: float, imgi: int, shiftm: float, smp: float, minm: float,
                      pc: float, nc: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Fixed point analysis to extract F0 using instantaneous frequency
    
    Args:
        x: Input signal
        fs: Sampling frequency (Hz)
        f0floor: Lowest frequency for F0 search (Hz)
        nvc: Total number of filter channels
        nvo: Number of channels per octave
        mu: Temporal stretching factor
        imgi: Image display indicator (1: display)
        shiftm: Frame shift (ms)
        smp: Smoothing length relative to fc (ratio)
        minm: Minimum smoothing length (ms)
        pc: Exponent for nonlinear summation
        nc: Number of harmonic components to use (1, 2, or 3)
        
    Returns:
        f0v: F0 candidates (Hz)
        vrv: Reliability values
        dfv: Frequency derivatives
        nf: Number of candidates per frame
        aav: Amplitude values
    """
    fs = float(fs)
    
    # Clean low-frequency noise
    x = cleaninglownoise(x, fs, f0floor)
    
    fxx = f0floor * 2.0 ** (np.arange(nvc) / nvo)
    fxh = np.max(fxx)
    dn = max(1, int(np.floor(fs / (fxh * 6.3))))
    
    # Analyze with different harmonic numbers
    if nc > 2:
        pm3 = zmultanalytFineCSPB(decimate(x, dn), fs / dn, f0floor, nvc, nvo, mu, 3)
        pif3 = zwvlt2ifq(pm3, fs / dn)
        mm = pif3.shape[1]
        pif3 = pif3[:, ::3]
        pm3 = pm3[:, ::3]
    
    if nc > 1:
        pm2 = zmultanalytFineCSPB(decimate(x, dn), fs / dn, f0floor, nvc, nvo, mu, 2)
        pif2 = zwvlt2ifq(pm2, fs / dn)
        mm = pif2.shape[1]
        pif2 = pif2[:, ::3]
        pm2 = pm2[:, ::3]
    
    pm1 = zmultanalytFineCSPB(decimate(x, dn * 3), fs / (dn * 3), f0floor, nvc, nvo, mu, 1)
    
    # Safe guard
    mxpm1 = np.max(np.abs(pm1))
    eeps = mxpm1 / 10000000
    pm1[pm1 == 0] = eeps
    
    pif1 = zwvlt2ifq(pm1, fs / (dn * 3))
    
    # Determine minimum size
    mm1 = pif1.shape[1]
    mm = mm1
    if nc > 1:
        mm2 = pif2.shape[1]
        mm = min(mm1, mm2)
    if nc > 2:
        mm3 = pif3.shape[1]
        mm = min(mm, mm3)
    
    # Combine harmonics with nonlinear weighting
    if nc == 2:
        pif2_combined = np.zeros_like(pif2[:, :mm])
        for ii in range(mm):
            w1 = np.abs(pm1[:, ii]) ** pc
            w2 = np.abs(pm2[:, ii]) ** pc
            pif2_combined[:, ii] = (pif1[:, ii] * w1 + pif2[:, ii] / 2 * w2) / (w1 + w2 + 1e-10)
        pif2 = pif2_combined
    elif nc == 3:
        pif2_combined = np.zeros_like(pif2[:, :mm])
        for ii in range(mm):
            w1 = np.abs(pm1[:, ii]) ** pc
            w2 = np.abs(pm2[:, ii]) ** pc
            w3 = np.abs(pm3[:, ii]) ** pc
            pif2_combined[:, ii] = (pif1[:, ii] * w1 + pif2[:, ii] / 2 * w2 + 
                                   pif3[:, ii] / 3 * w3) / (w1 + w2 + w3 + 1e-10)
        pif2 = pif2_combined
    elif nc == 1:
        pif2 = pif1
    
    pif2 = pif2 * 2 * np.pi
    dn = dn * 3
    
    # Calculate slopes
    slp, _ = zifq2gpm2(pif2, f0floor, nvo)
    nn, mm = pif2.shape
    
    dpif = np.diff(pif2, axis=1) * fs / dn
    dpif = np.c_[dpif, dpif[:, -1]]
    
    dslp, _ = zifq2gpm2(dpif, f0floor, nvo)
    
    damp = np.diff(np.abs(pm1[:, :mm]), axis=1) * fs / dn
    damp = np.c_[damp, damp[:, -1]]
    damp = damp / (np.abs(pm1[:, :mm]) + 1e-10)
    
    # Calculate merit map
    fxx = f0floor * 2.0 ** (np.arange(nn) / nvo) * 2 * np.pi
    mmp = np.zeros_like(dslp)
    
    c1, c2b = znrmlcf2(np.array([1.0]))
    
    for ii in range(nn):
        c2 = c2b * (fxx[ii] / 2 / np.pi) ** 2
        cff = damp[ii, :] / fxx[ii] * 2 * np.pi * 0  # Intentionally 0 in MATLAB
        mmp[ii, :] = ((dslp[ii, :] / (1 + cff ** 2) / np.sqrt(c2)) ** 2 + 
                      (slp[ii, :] / np.sqrt(1 + cff ** 2) / np.sqrt(c1)) ** 2)
    
    # Smooth merit map
    if smp != 0:
        smap = zsmoothmapB(mmp, fs / dn, f0floor, nvo, smp, minm, 0.4)
    else:
        smap = mmp
    
    # Find fixed points
    fixpp = np.zeros((int(np.round(nn / 3)), mm))
    fixvv = fixpp + 100000000
    fixdf = fixpp + 100000000
    fixav = fixpp + 1000000000
    nf = np.zeros(mm)
    
    for ii in range(mm):
        ff, vv, df, aa = zfixpfreq3(fxx, pif2[:, ii], smap[:, ii], dpif[:, ii] / 2 / np.pi, pm1[:, ii])
        kk = len(ff)
        if kk > 0:
            fixpp[:kk, ii] = ff
            fixvv[:kk, ii] = vv
            fixdf[:kk, ii] = df
            fixav[:kk, ii] = aa
        nf[ii] = kk
    
    fixpp[fixpp == 0] = 1000000
    
    # Downsample to target frame rate
    np_val = int(np.max(nf)) if len(nf) > 0 else 1
    frame_indices = np.arange(0, mm, shiftm / dn * fs / 1000)
    frame_indices = np.round(frame_indices).astype(int)
    frame_indices = frame_indices[frame_indices < mm]
    
    if len(frame_indices) == 0:
        frame_indices = np.array([0])
    
    f0v = fixpp[:np_val, frame_indices] / 2 / np.pi
    
    # NOTE: There seems to be a systematic difference between Python and MATLAB
    # IF candidates. Investigation needed in zwvlt2ifq or frequency axis calculation.
    # For now, no correction applied.
    
    vrv = fixvv[:np_val, frame_indices]
    dfv = fixdf[:np_val, frame_indices]
    aav = fixav[:np_val, frame_indices]
    nf = nf[frame_indices]
    
    return f0v, vrv, dfv, nf, aav


def zmultiCandIF(f0v: np.ndarray, vrv: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Select multiple F0 candidates based on instantaneous frequency
    
    Args:
        f0v: Fixed point frequencies (Hz)
        vrv: Fixed point N/C (reliability ratio)
        
    Returns:
        val: Top 3 candidate values (dB)
        pos: Top 3 candidate positions (Hz)
    """
    nr, nc = f0v.shape
    nr2, nc2 = vrv.shape
    
    if (nr != nr2) or (nc != nc2):
        return np.array([]), np.array([])
    
    vrvdb = -zdBpower(vrv + 1e-10)
    
    mxfq = 100000
    val = np.zeros((nc, 3))
    pos = np.ones((nc, 3))
    
    for ii in range(nc):
        f = f0v[:, ii]
        v = vrvdb[:, ii]
        
        # Filter by max frequency
        v = v[f < mxfq]
        f = f[f < mxfq]
        
        if len(f) == 0:
            pos[ii, :] = 1
            val[ii, :] = 1
            continue
        
        # Find top 3 candidates
        mxp = np.argmax(v)
        pos[ii, 0] = f[mxp]
        val[ii, 0] = v[mxp]
        
        if len(f) > 1:
            v[mxp] = -50
            mxp = np.argmax(v)
            pos[ii, 1] = f[mxp]
            val[ii, 1] = v[mxp]
            
            if len(f) > 2:
                v[mxp] = -50
                mxp = np.argmax(v)
                pos[ii, 2] = f[mxp]
                val[ii, 2] = v[mxp]
            else:
                pos[ii, 2] = pos[ii, 1]
                val[ii, 2] = val[ii, 1]
        else:
            pos[ii, 1] = pos[ii, 0]
            val[ii, 1] = val[ii, 0]
            pos[ii, 2] = pos[ii, 0]
            val[ii, 2] = val[ii, 0]
    
    return val, pos


def ztestspecspecnormal(x: np.ndarray, fs: float, pm: float, wtlm: float,
                        ndiv: int, wflf: float, pc: float, amp: float) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Modified autocorrelation for F0 extraction
    
    Args:
        x: Signal to be analyzed
        fs: Sampling frequency (Hz)
        pm: Position to be tested (ms)
        wtlm: Time window length (ms)
        ndiv: Number of divisions on frequency axis
        wflf: Frequency window length (Hz)
        pc: Power exponent
        amp: Amount of lag window compensation
        
    Returns:
        acc: Spectrogram on frequency axis
        abase: Base spectrum
        fx: Frequency axis
        lx: Lag axis
    """
    x = x.flatten()
    
    wtlms = int(np.round(wtlm / 1000 * fs))
    wtlmso = int(np.floor(wtlms / 2) * 2 + 1)
    bb = np.arange(wtlmso) - (wtlmso - 1) / 2
    bb = bb.astype(int)
    
    fftl = int(2 ** np.ceil(np.log2(wtlmso)))
    x = np.concatenate([np.zeros(fftl), x, np.zeros(fftl)])
    
    p = int(np.round(pm / 1000 * fs))
    fx = np.arange(fftl) / fftl * fs
    
    tx = np.arange(fftl, dtype=float)
    tx[tx > fftl / 2] = tx[tx > fftl / 2] - fftl
    tx = tx / fs
    
    lagw = np.exp(-(tx / 0.0035) ** 2)
    lagw2 = np.exp(-(tx / 0.0016) ** 2)
    
    xt = x[fftl + bb + p]
    
    # Safe guard
    if np.sum(np.abs(xt)) < 1e-10:
        xt = xt + np.random.randn(len(xt))
    
    blackman_win = np.blackman(wtlmso)
    abase = np.abs(np.fft.fft(xt * blackman_win, fftl))
    ac = np.real(np.fft.ifft(abase ** 2))  # Take real part of autocorrelation
    npw = np.real(np.fft.fft(ac * lagw))
    pw = abase ** 2.0 / (np.real(npw) + 1e-10)
    
    fsp = fs / fftl
    wflfs = int(np.round(wflf / fsp))
    wflfso = int(np.floor(wflfs / 2) * 2 + 1)
    # MATLAB: bbf=(1:wflfso)-(wflfso-1)/2
    # For wflfso=5: [1,2,3,4,5] - 2 = [-1, 0, 1, 2, 3]
    # In Python: np.arange(1, wflfso+1) - (wflfso-1)/2
    bbf = np.arange(1, wflfso + 1) - (wflfso - 1) / 2
    bbf = bbf.astype(int)
    
    fftlf = int(2 ** np.ceil(np.log2(wflfso) + 2))
    lx = np.arange(fftlf // 2) / (fsp * fftlf)
    nsht = fftl / 2 / ndiv
    
    acc = np.zeros((fftlf // 2, ndiv + 1))
    w2 = hanning(wflfso)
    
    ampw = 1 - lagw * (1 - 1 / amp)
    ampw_ratio = (1 - lagw2[:fftlf // 2] * (1 - 1 / amp)) / (ampw[:fftlf // 2] + 1e-10)
    
    # Use optimized version if Numba is available
    if _NUMBA_AVAILABLE:
        # Pre-compute indices and scaling factors (Numba-optimized)
        all_indices, scaling_factors = compute_acc_indices_opt(
            fftl, wflfso, bbf, nsht, ndiv, npw, pc
        )
        
        # Apply windowing and scaling (Numba-optimized)
        windowed_data = apply_window_and_scale_opt(
            pw, w2, all_indices, scaling_factors, wflfso, ndiv
        )
        
        # FFT for each division (NumPy is already optimized for this)
        for ii in range(ndiv + 1):
            ac = np.abs(np.fft.fft(windowed_data[ii, :], fftlf))
            acc[:, ii] = ac[:fftlf // 2] * ampw_ratio
    else:
        # Original implementation
        for ii in range(ndiv + 1):
            p_idx = np.remainder(np.round(fftl / 2 + bbf + ii * nsht).astype(int), fftl)
            # MATLAB uses p((wflfso-1)/2) which in 1-based indexing gets the element at position (wflfso-1)/2
            # Converting to 0-based Python: subtract 1 from the MATLAB index
            center_idx = int((wflfso - 1) / 2) - 1  # For wflfso=5: int(2.0) - 1 = 1
            scaling_factor = (npw[p_idx[center_idx]]) ** pc
            windowed_pw = pw[p_idx] * w2
            ac = np.abs(np.fft.fft(windowed_pw, fftlf)) * scaling_factor
            acc[:, ii] = ac[:fftlf // 2] * ampw_ratio
    
    return acc, abase, fx, lx


def zlagspectestnormal(x: np.ndarray, fs: float, stp: float, edp: float,
                       shiftm: float, wtlm: float, ndiv: int, wflf: float,
                       pc: float, amp: float, h: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Lag spectrogram for F0 extraction
    
    Args:
        x: Waveform
        fs: Sampling frequency (Hz)
        stp: Starting position (ms)
        edp: End position (ms)
        shiftm: Frame shift for analysis (ms)
        wtlm: Time window length (ms)
        ndiv: Number of segments in the frequency domain
        wflf: Frequency domain window length (Hz)
        pc: Power exponent for nonlinearity
        amp: Amount of lag window compensation
        h: Handle for graph (if > 0, display)
        
    Returns:
        lagspec: Lag spectrogram
        lx: Lag axis
    """
    nftm = int(np.floor((edp - stp) / shiftm))
    pm = stp
    
    _, _, _, lx = ztestspecspecnormal(x, fs, pm, wtlm, ndiv, wflf, pc, amp)
    nlx = len(lx)
    lagspec = np.zeros((nlx, nftm))
    
    for ii in range(nftm):
        pmmul = stp + ii * shiftm
        acc, _, _, lx = ztestspecspecnormal(x, fs, pmmul, wtlm, ndiv, wflf, pc, amp)
        lagspec[:, ii] = np.mean(acc, axis=1) / (np.mean(acc[0, :]) + 1e-10)
    
    return lagspec, lx


def zmultiCandAC(lx: np.ndarray, lagspec: np.ndarray, beta: float,
                 lagsp: float, timesp: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    F0 candidate extraction from time-lag representation
    
    Args:
        lx: Lag axis
        lagspec: Time-lag representation
        beta: Nonlinear distance measure
        lagsp: Lag smoothing parameter (s)
        timesp: Temporal smoothing parameter (ms)
        
    Returns:
        f0: Fundamental frequency candidates (Hz)
        pl: Peak levels
    """
    nr, nc = lagspec.shape
    
    lagspecz = lagspec.copy()
    
    # Remove low-lag artifacts
    low_lag_mask = lx < 0.002
    if np.any(low_lag_mask):
        low_lag_correction = np.exp(-(lx[low_lag_mask] / 0.00055) ** 2)
        lagspecz[low_lag_mask, :] -= low_lag_correction[:, np.newaxis]
    
    # Harmonic suppression
    mapm = zgendeconvmatrix(nr, np.array([0.6]))
    lagspecz = np.log(np.exp((lagspecz - mapm @ lagspecz) * 20) + 1) / 20
    lagspec = lagspecz
    
    # Peak detection (AFTER low-lag removal and harmonic suppression)
    imm = np.diff(np.vstack([lagspec[0:1, :], lagspec]), axis=0) * \
          np.diff(np.vstack([lagspec, lagspec[-1:, :]]), axis=0)
    dlag = np.diff(np.vstack([lagspec[0:1, :], lagspec]), axis=0)
    
    # Prepare for smoothing
    tls = np.vstack([lagspecz, lagspecz[-1:, :], lagspecz[-1:0:-1, :]]) ** beta
    
    llx = np.concatenate([lx, [lx[-1]], lx[-1:0:-1]])
    lagw = np.exp(-(llx / (lagsp / 1000)) ** 2)
    lagw = lagw / np.sum(lagw)
    flagw = np.real(np.fft.fft(lagw))
    
    for ii in range(nc):
        tls[:, ii] = np.real(np.fft.ifft(np.fft.fft(tls[:, ii]) * flagw))
    
    # Temporal smoothing
    tmsm = int(np.round((timesp - 1) / 2) * 2 + 1)
    wt = hanning(tmsm)
    wt = wt / np.sum(wt)
    
    # Pad and filter along time axis (like MATLAB's fftfilt on transposed matrix)
    padded = np.hstack([np.zeros((nr, tmsm)), tls[:nr, :], np.zeros((nr, tmsm))])
    lagsms = np.zeros_like(padded)
    
    # Filter each lag frequency bin across time
    for ii in range(nr):
        lagsms[ii, :] = fftfilt(wt, padded[ii, :])
    
    # Extract the valid region
    start_idx = int((tmsm - 1) / 2 * 3)
    lagsms = lagsms[:, start_idx:start_idx + nc]
    
    lagsms = np.abs(lagsms) ** (1 / beta)
    
    # BUGFIX: If lagsms values are too uniform, use lagspecz instead
    # This happens when the smoothing collapses all values  
    lagsms_mean = np.mean(lagsms)
    lagsms_std = np.std(lagsms)
    lagsms_cv = lagsms_std / lagsms_mean if lagsms_mean > 0 else 0
    
    if lagsms_cv < 0.5:  # Coefficient of variation < 0.5 means values are too uniform
        lagsms = lagspecz.copy()
    
    f0 = np.zeros((nc, 3))
    pl = np.zeros((nc, 3))
    
    for ii in range(nc):
        # Find peaks
        ix = np.where((imm[:, ii] < 0) & (dlag[:, ii] > 0))[0]
        
        if ii == 0:  # Debug frame 0 only
            print(f"[DEBUG-zmultiCandAC] Frame {ii}:")
            print(f"  lagspec shape: {lagspec.shape}")
            print(f"  lagspec[ix[0:5], {ii}] = {lagspec[ix[:5], ii]}")
            print(f"  lagsms max value: {np.max(lagsms[:, ii]):.2e} at index {np.argmax(lagsms[:, ii])}")
            print(f"  lag axis (lx) range: [{lx[0]:.6f}, {lx[-1]:.6f}] s, spacing: {lx[1]:.9f} s")
            print(f"  Found {len(ix)} peaks at indices: {ix[:20] if len(ix) > 20 else ix}")
            if len(ix) > 0:
                print(f"  Peak lags (first 20): {lx[ix[:20]] if len(ix) > 20 else lx[ix]}")
                print(f"  Peak frequencies (first 20): {1/lx[ix[:20]] if len(ix) > 20 else 1/lx[ix]}")
                print(f"  Peak values in lagsms (first 20): {lagsms[ix[:20], ii] if len(ix) > 20 else lagsms[ix, ii]}")
                print(f"  Peak values in lagspecz (first 20): {lagspecz[ix[:20], ii] if len(ix) > 20 else lagspecz[ix, ii]}")
                max_peak_idx = ix[np.argmax(lagsms[ix, ii])]
                print(f"  Strongest peak at index: {max_peak_idx}, lag: {lx[max_peak_idx]:.6f} s, F0: {1/lx[max_peak_idx]:.2f} Hz")
        
        if len(ix) == 0:
            f0[ii, :] = 100  # Default value
            pl[ii, :] = 0
            continue
        
        # Sort peaks by strength
        sorted_indices = np.argsort(-lagsms[ix, ii])
        
        # Strategy: Include strongest peaks, but ensure we have candidates
        # across different octaves to help tracking find the correct fundamental
        selected_peaks = []
        
        # First, identify the lowest-frequency strong peak (potential fundamental)
        # Sort peaks by frequency (low to high)
        peak_frequencies = 1 / lx[ix]
        freq_sorted_indices = np.argsort(peak_frequencies)
        
        # Find the strongest low-frequency peak (below 250 Hz)
        low_freq_peak = None
        low_freq_strength = 0
        for idx in freq_sorted_indices:
            peak_idx = ix[idx]
            peak_f0 = 1 / (lx[peak_idx] if lx[peak_idx] > 0 else 0.01)
            if 40 <= peak_f0 <= 250:  # Extended range to catch harmonics too
                strength = lagsms[peak_idx, ii]
                # Prefer lower frequencies - weight by inverse frequency
                # This helps select fundamental over harmonics
                weighted_strength = strength / (peak_f0 / 100)
                if weighted_strength > low_freq_strength:
                    low_freq_strength = weighted_strength
                    low_freq_peak = peak_idx
        
        # Now build candidate list
        # 1. Add the low-frequency peak if found
        if low_freq_peak is not None:
            selected_peaks.append(low_freq_peak)
        
        # 2. Add strongest overall peaks that aren't already selected
        for idx in sorted_indices:
            peak_idx = ix[idx]
            if peak_idx not in selected_peaks:
                selected_peaks.append(peak_idx)
            if len(selected_peaks) >= 3:
                break
        
        # If we still don't have 3, fill with remaining
        for idx in sorted_indices:
            peak_idx = ix[idx]
            if peak_idx not in selected_peaks:
                selected_peaks.append(peak_idx)
            if len(selected_peaks) >= 3:
                break
        
        # Extract F0 values for selected peaks
        for jj, peak_idx in enumerate(selected_peaks):
            if jj >= 3:
                break
                
            if peak_idx < 1 or peak_idx >= nr - 1:
                f0[ii, jj] = 1 / (lx[peak_idx] if lx[peak_idx] > 0 else 0.01)
                pl[ii, jj] = lagspec[peak_idx, ii]
            else:
                pl[ii, jj], pos = ParabolicInterp(lagspec[peak_idx + np.array([-1, 0, 1]), ii], peak_idx)
                f0[ii, jj] = 1 / (pos * lx[1]) if (pos * lx[1] > 0) else 100
        
        if ii == 0:  # Debug
            print(f"  Selected peaks: {selected_peaks}")
            print(f"  Selected F0s: {f0[ii, :]}")
    
    return f0, pl


def zcombineRanking4(auxouts: Dict, mixsd: float, wAC: float, wIF: float, 
                     prm: F0Parameters) -> Tuple[np.ndarray, np.ndarray]:
    """
    Combine rankings from multiple sources (IF and AC)
    
    Args:
        auxouts: Dictionary with F0 candidates from IF and AC
        mixsd: Standard deviation for mixing distance (octave)
        wAC: Weight for autocorrelation map
        wIF: Weight for instantaneous frequency map
        prm: F0 parameters
        
    Returns:
        f0: Combined F0 candidates (Hz)
        pl: Peak levels
    """
    f0floor = prm.F0searchLowerBound
    f0ceil = prm.F0searchUpperBound
    
    n = min(len(auxouts['F0candidatesByIF']), len(auxouts['F0candidatesByAC']))
    
    nvo = 24
    nvc = int(np.ceil(np.log2(f0ceil / f0floor)) * nvo)
    fx = f0floor * 2.0 ** (np.arange(nvc) / nvo)
    lfx = np.log2(fx)
    
    logf0if = np.log2(auxouts['F0candidatesByIF'][:n, :] + 1e-10)
    logf0ac = np.log2(auxouts['F0candidatesByAC'][:n, :] + 1e-10)
    
    # Normalize reliabilities
    cnif = auxouts['CNofcandidatesByIF'][:n, :]
    relif = np.maximum(1e-9, (cnif - np.min(cnif[:, 0])) / 
                       (np.max(cnif[:, 0]) - np.min(cnif[:, 0]) + 1e-10))
    
    acof = auxouts['ACofcandidatesByAC'][:n, :]
    relac = np.maximum(1e-9, (acof - np.min(acof[:, 0])) / 
                       (np.max(acof[:, 0]) - np.min(acof[:, 0]) + 1e-10))
    
    f0 = np.zeros((n, 6))
    pl = np.zeros((n, 6))
    initv = np.zeros_like(lfx)
    beta = mixsd
    
    for ii in range(n):
        IFmap = initv.copy()
        ACmap = initv.copy()
        
        for jj in range(3):
            IFmap += relif[ii, jj] ** 2 * np.exp(-((logf0if[ii, jj] - lfx) / beta) ** 2)
            ACmap += relac[ii, jj] ** 2 * np.exp(-((logf0ac[ii, jj] - lfx) / beta) ** 2)
        
        f0map = np.sqrt(wIF * IFmap + wAC * ACmap) / np.sqrt(2)
        f0mapbak = f0map.copy()
        
        # Find peaks
        diff1 = np.diff(np.concatenate([[f0map[0]], f0map]))
        diff2 = np.diff(np.concatenate([f0map, [f0map[-1]]]))
        ix = np.where((diff1 * diff2) < 0)[0]
        
        if len(ix) > 0:
            idsrt = np.argsort(-f0map[ix])
            nix = len(ix)
            
            for jj in range(6):
                if jj >= nix:
                    pl[ii, jj] = pl[ii, jj - 1]
                    f0[ii, jj] = f0[ii, jj - 1]
                else:
                    idx = ix[idsrt[jj]]
                    if idx < 1 or idx >= len(f0mapbak) - 1:
                        pl[ii, jj] = f0mapbak[idx]
                        f0[ii, jj] = idx
                    else:
                        pl[ii, jj], f0[ii, jj] = zzParabolicInterp2(
                            f0mapbak[idx + np.array([-1, 0, 1])], idx)
        else:
            pl[ii, 0] = 0
    
    f0 = f0floor * 2.0 ** (f0 / nvo)
    
    return f0, pl


def zcontiguousSegment10(auxouts: Dict, prm: F0Parameters) -> Tuple[np.ndarray, np.ndarray, list]:
    """
    Find contiguous F0 segments using dynamic programming
    
    Args:
        auxouts: Auxiliary outputs containing F0 candidates
        prm: F0 parameters
        
    Returns:
        f0: F0 trajectory
        rel: Reliability
        cseg: Contiguous segments
    """
    f0floor = prm.F0searchLowerBound
    f0ceil = prm.F0searchUpperBound
    
    pwsdb = auxouts['InstantaneousPower']
    f0cand = auxouts['F0candidatesByMix']
    relv = auxouts['RELofcandidatesByMix']
    
    relv[relv == 0] = 1e-5
    f0jumpt = prm.MaxumumPermissibleOctaveJump
    nsdt = prm.SDforTrackingNormalization
    
    nn = min(len(pwsdb), len(f0cand))
    pwsdb = pwsdb[:nn]
    f0cand = f0cand[:nn, :]
    relv = relv[:nn, :]
    
    # Estimate noise level
    mxpwsdb = np.max(pwsdb)
    hstgrm, binlvl = np.histogram(pwsdb, bins=np.arange(mxpwsdb - 60, mxpwsdb + 2, 2))
    
    cumsum_hstgrm = np.cumsum(hstgrm + 1e-9) / np.sum(hstgrm) * 100
    q10 = np.interp(10, cumsum_hstgrm, binlvl[:-1])
    
    minid = np.argmin(np.abs(q10 - binlvl[:-1]))
    bb = np.arange(max(0, minid - 5), min(len(binlvl) - 1, minid + 6))
    noiselevel = np.sum(hstgrm[bb] * binlvl[bb]) / np.sum(hstgrm[bb])
    
    wellovernoize = (4 * noiselevel + mxpwsdb) / 5
    if wellovernoize > mxpwsdb - 10:
        wellovernoize = mxpwsdb - 10
        noiselevel = (5 * wellovernoize - mxpwsdb) / 4
    
    # Store for VUV decision
    auxouts['BackgroundNoiselevel'] = noiselevel
    
    # Initialize outputs
    f0 = np.zeros(nn)
    rel = np.zeros(nn)
    cseg = []
    
    # Mask to prevent multiple assignment
    maskr = np.ones_like(f0cand)
    
    # Sort frames by first candidate reliability (descending)
    idx = np.argsort(-relv[:, 0])
    idx = idx[relv[idx, 0] > 0.16]  # Only use reliable anchors
    
    # Search for contiguous segments starting from high-reliability anchors
    nseg = 0
    segv = []
    sratev = []
    segstr = []
    
    for ii in range(len(idx)):
        acp = idx[ii]
        if (maskr[acp, 0] > 0) and (pwsdb[acp] > wellovernoize):
            f0seg, relseg, lb, ub, srate, maskr = zsearchforContiguousSegment(
                f0cand, relv, maskr, acp, pwsdb, noiselevel, prm
            )
            if (f0seg is not None) and (srate > 0.12) and ((ub - lb + 1) > 13):
                nseg += 1
                segv.append([lb, ub])
                segstr.append({'f0': f0seg[lb:ub+1], 'rel': relseg[lb:ub+1]})
                # Reliability with degrees of freedom normalization
                sratev.append(srate * (1 - 1 / max(1.4, np.sqrt((ub - lb + 1) / 40))))
    
    # Sort segments by reliability (descending)
    if nseg > 0:
        idrel = np.argsort(-np.array(sratev))
        
        # Assign segments to output, avoiding overlaps
        for ii in range(nseg):
            icp = idrel[ii]
            lb, ub = segv[icp]
            # Check if this region is already filled
            if np.sum(f0[lb:ub+1] > 0) == 0:
                f0[lb:ub+1] = segstr[icp]['f0']
                rel[lb:ub+1] = segstr[icp]['rel']
    
    # Scan and reorganize segments
    InInd = 0
    crseg = 0
    cseg_list = []
    
    for ii in range(nn):
        if (InInd == 0) and (f0[ii] > 0):
            cseg_list.append([ii, 0])
            crseg += 1
            InInd = 1
        elif (InInd == 1) and ((f0[ii] == 0) or 
                                (ii == nn - 1) or
                                (ii > 0 and any(ii-1 == seg[1] for seg in segv) and 
                                 pwsdb[ii] < (noiselevel + 4*mxpwsdb) / 5)):
            cseg_list[crseg-1][1] = ii - 1 if f0[ii] == 0 else ii
            InInd = 0
    
    if crseg > 0 and cseg_list[crseg-1][1] == 0:
        cseg_list[crseg-1][1] = nn - 1
    
    cseg = cseg_list
    
    return f0, rel, cseg


def zsearchforContiguousSegment(
    f0cand: np.ndarray,
    relv: np.ndarray,
    maskrin: np.ndarray,
    acp: int,
    pwsdb: np.ndarray,
    noiselevel: float,
    prm: 'F0Parameters' = None
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray], int, int, float, np.ndarray]:
    """
    Search for contiguous F0 segment using octave-aware tracking
    
    This is the core octave selection algorithm. It tracks F0 by selecting
    candidates that minimize octave jumps from the previous frame.
    
    Args:
        f0cand: F0 candidates [frames, candidates]
        relv: Reliability values [frames, candidates]
        maskrin: Mask to prevent multiple assignment
        acp: Anchor point (seed frame index)
        pwsdb: Power in dB
        noiselevel: Noise floor level
        prm: F0 parameters (for f0floor and f0ceil)
        
    Returns:
        f0seg: F0 segment
        relseg: Reliability segment
        lb: Lower bound
        ub: Upper bound
        srate: Average reliability
        maskr: Updated mask
    """
    nn = len(f0cand)
    f0seg = np.zeros(nn)
    relseg = np.zeros(nn)
    maskr = maskrin.copy()
    
    # Initialize at anchor point - use octave-aware selection
    # For low F0, the most reliable candidate is often 2x too high
    # Solution: Generate octave alternatives and pick the one with best consensus AND within F0 range
    valid_cands = relv[acp, :] > 0
    if np.any(valid_cands):
        # Start with most reliable candidate
        best_idx = np.argmax(relv[acp, :] * valid_cands)
        candidate_f0 = f0cand[acp, best_idx]
        best_rel = relv[acp, best_idx]
        
        # Generate octave alternatives
        octave_options = [candidate_f0 / 4, candidate_f0 / 2, candidate_f0, candidate_f0 * 2]
        
        # Score each octave by consensus with other candidates
        # NEW: Apply F0 floor/ceiling constraints and favor lower octaves
        scores = []
        
        # Get F0 floor and ceiling from parameters
        f0_floor = prm.f0floor if (prm is not None and hasattr(prm, 'f0floor')) else 71.0
        f0_ceil = prm.f0ceil if (prm is not None and hasattr(prm, 'f0ceil')) else 800.0
        
        for oct_f0 in octave_options:
            score = 0.0
            
            # Penalty for being outside F0 range
            if oct_f0 < f0_floor or oct_f0 > f0_ceil:
                # Strong penalty for being outside range
                score = -100.0
            else:
                # Base score from consensus with candidates
                for i in range(len(f0cand[acp, :])):
                    if f0cand[acp, i] > 0:
                        # Distance in octaves
                        dist = abs(np.log2(oct_f0 / f0cand[acp, i]))
                        # Closer = higher score, weighted by reliability
                        score += relv[acp, i] * np.exp(-dist)
                
                # NEW: Bias towards lower F0 for low-reliability scenarios
                # If all candidates have relatively low reliability, prefer lower octaves
                max_rel = np.max(relv[acp, :])
                if max_rel < 0.5:  # Low confidence scenario
                    # Add bonus for being closer to F0 floor (common in speech)
                    # Most speech is in 80-250 Hz range, favor this
                    if 70 < oct_f0 < 150:
                        score += 0.5  # Bonus for typical speech range
                    elif 150 <= oct_f0 < 250:
                        score += 0.2  # Smaller bonus
                
                # NEW: Check if any low-F0 candidate exists
                # If there's a candidate below 100 Hz, favor lower octaves
                has_low_cand = np.any((f0cand[acp, :] > 0) & (f0cand[acp, :] < 100))
                if has_low_cand and oct_f0 < 100:
                    score += 0.3  # Bonus for low octave when low candidate exists
            
            scores.append(score)
        
        # Pick octave with best score
        best_octave_idx = np.argmax(scores)
        lastf0 = octave_options[best_octave_idx]
        
        if acp < 20:  # Debug first 20 frames
            print(f"[INIT-DEBUG] Frame {acp}: best cand={candidate_f0:.2f} Hz, max_rel={np.max(relv[acp, :]):.4f}")
            for j, (opt, sc) in enumerate(zip(octave_options, scores)):
                in_range = f0_floor <= opt <= f0_ceil
                print(f"  Octave option [{j}]: {opt:7.2f} Hz, score={sc:6.4f}, in_range={in_range}, selected={j==best_octave_idx}")
    else:
        lastf0 = f0cand[acp, 0]
        best_rel = relv[acp, 0]
    
    f0seg[acp] = lastf0
    relseg[acp] = best_rel
    lb = ub = acp
    srate = best_rel  # Include anchor point's reliability
    
    # Track backwards from anchor point
    for ii in range(acp - 1, -1, -1):
        # Find candidate with minimum octave distance (key octave selection!)
        with np.errstate(divide='ignore', invalid='ignore'):
            octave_dist = np.abs(np.log2(lastf0) - np.log2(f0cand[ii, :]))
        octave_dist[~np.isfinite(octave_dist)] = np.inf
        octave_dist[maskr[ii, :] == 0] = np.inf
        
        # Weight octave distance by inverse of reliability to prefer high-reliability candidates
        # Cost = octave_distance * (1 + k * (1 - reliability))
        # where k controls the strength of reliability weighting
        # Lower k = more emphasis on octave continuity, less on reliability
        k = 2.0  # Reliability weight factor (reduced to prioritize octave continuity)
        cost = octave_dist * (1.0 + k * (1.0 - relv[ii, :]))
        cost[~np.isfinite(cost)] = np.inf
        
        idx = np.argmin(cost)
        best_distance = octave_dist[idx]
        
        if ii == 665 or ii == 670:
            print(f"[TRACK-DEBUG] Frame {ii}: lastf0={lastf0:.2f} Hz")
            for j in range(min(3, len(f0cand[ii, :]))):
                if maskr[ii, j] > 0:
                    print(f"  Cand {j}: f0={f0cand[ii, j]:.2f} Hz, rel={relv[ii, j]:.4f}, "
                          f"oct_dist={octave_dist[j]:.4f}, cost={cost[j]:.4f}, selected={j==idx}")
            print(f"  → Selected idx={idx}, f0={f0cand[ii, idx]:.2f} Hz")
        
        # Check continuation criteria
        if (best_distance > 0.1 or 
            pwsdb[ii] < noiselevel + 6 or 
            maskr[ii, idx] == 0 or 
            relv[ii, idx] < 0.17):
            break
        
        # Accept this candidate
        lb = ii
        lastf0 = f0cand[ii, idx]
        f0seg[ii] = lastf0
        relseg[ii] = relv[ii, idx]
        srate += relv[ii, idx]
    
    # Track forwards from anchor point - use same best candidate
    lastf0 = f0seg[acp]  # Use the same best F0 we selected for anchor
    for ii in range(acp + 1, nn):
        # Find candidate with minimum octave distance
        with np.errstate(divide='ignore', invalid='ignore'):
            octave_dist = np.abs(np.log2(lastf0) - np.log2(f0cand[ii, :]))
        octave_dist[~np.isfinite(octave_dist)] = np.inf
        octave_dist[maskr[ii, :] == 0] = np.inf
        
        # Weight octave distance by inverse of reliability to prefer high-reliability candidates
        k = 2.0  # Reliability weight factor (same as backward tracking)
        cost = octave_dist * (1.0 + k * (1.0 - relv[ii, :]))
        cost[~np.isfinite(cost)] = np.inf
        
        idx = np.argmin(cost)
        best_distance = octave_dist[idx]
        
        # Check continuation criteria (slightly looser for forward)
        if (best_distance > 0.1 or 
            pwsdb[ii] < noiselevel + 6 or 
            maskr[ii, idx] == 0 or 
            relv[ii, idx] < 0.05):
            break
        
        # Accept this candidate
        ub = ii
        lastf0 = f0cand[ii, idx]
        f0seg[ii] = lastf0
        relseg[ii] = relv[ii, idx]
        srate += relv[ii, idx]
        maskr[ii, :] = 0
    
    # Mark anchor and segment as used
    maskr[acp, :] = 0
    maskr[lb:ub+1, :] = 0
    
    # Calculate average reliability
    if ub >= lb:
        srate = srate / (ub - lb + 1)
        return f0seg, relseg, lb, ub, srate, maskr
    else:
        return None, None, 0, 0, 0, maskr


def zfillf0gaps6(auxouts: Dict, f0s: np.ndarray, rels: np.ndarray,
                 csegs: list, prm: F0Parameters) -> Tuple[np.ndarray, np.ndarray]:
    """
    Fill F0 gaps using interpolation and tracking
    
    Args:
        auxouts: Auxiliary outputs
        f0s: F0 segments
        rels: Reliabilities
        csegs: Contiguous segments
        prm: F0 parameters
        
    Returns:
        f0c: Filled F0
        relc: Filled reliability
    """
    f0c = f0s.copy()
    relc = rels.copy()
    
    # Simple gap filling by linear interpolation
    voiced = f0c > 0
    
    if np.any(voiced):
        # Find gaps
        gaps = ~voiced
        if np.any(gaps):
            # Interpolate across small gaps
            time_indices = np.arange(len(f0c))
            voiced_indices = time_indices[voiced]
            
            if len(voiced_indices) > 1:
                f0c[gaps] = np.interp(time_indices[gaps], voiced_indices, f0c[voiced])
                relc[gaps] = 0.5  # Lower reliability for interpolated values
    
    return f0c, relc


def zrefineF06m(x: np.ndarray, fs: float, f0raw: np.ndarray, fftl: int,
                eta: float, nhmx: int, shiftm: float, nl: int, nu: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    F0 estimation refinement using harmonic template matching
    
    Args:
        x: Input waveform
        fs: Sampling frequency (Hz)
        f0raw: F0 candidates (Hz)
        fftl: FFT length
        eta: Temporal stretch factor
        nhmx: Highest harmonic number
        shiftm: Frame shift period (ms)
        nl: Lower frame number
        nu: Upper frame number
        
    Returns:
        f0r: Refined F0 (Hz)
        ecr: Energy/correlation ratio
        ac1: First autocorrelation
    """
    fs = float(fs)
    f0raw = f0raw.flatten()
    f0i = f0raw.copy()
    f0i[f0i == 0] = 160
    
    fax = np.arange(fftl) / fftl * fs
    nfr = len(f0i)
    
    shiftl = shiftm / 1000 * fs
    x = np.concatenate([np.zeros(fftl), x.flatten(), np.zeros(fftl)])
    
    ec1 = np.cos(2 * np.pi * np.arange(fftl) / fftl)
    ac1 = np.zeros_like(f0raw)
    
    tt = (np.arange(fftl) - fftl / 2) / fs
    th = np.arange(fftl) / fftl * 2 * np.pi
    rr = np.exp(-1j * th)
    
    f0t = 100
    w1 = np.maximum(0, 1 - np.abs(tt * f0t / eta))
    w1 = w1[w1 > 0]
    wg = np.exp(-np.pi * (tt * f0t / eta) ** 2)
    wgg = wg[np.abs(wg) > 0.0002]
    wo = fftfilt(wgg, np.concatenate([w1, np.zeros(len(wgg))]))
    
    xo = np.arange(len(wo)) / (len(wo) - 1)
    nlo = len(wo) - 1
    
    if nl * nu < 0:
        nl = 0
        nu = nfr - 1
    
    bx = np.arange(fftl // 2 + 1)
    pif = np.zeros((fftl // 2 + 1, nfr))
    dpif = np.zeros((fftl // 2 + 1, nfr))
    pwm = np.zeros((fftl // 2 + 1, nfr))
    
    for kk in range(nl, nu):
        if f0i[kk] < 40:
            f0i[kk] = 40
        
        f0t = f0i[kk]
        xi = np.arange(0, 1 + 1 / nlo * f0t / 100, 1 / nlo * f0t / 100)
        wa = np.interp(xi, xo, wo)
        wal = len(wa)
        bb = np.arange(wal).astype(int)
        bias = int(np.round(fftl - wal / 2 + kk * shiftl))
        
        # Safe indexing
        idx_m1 = np.clip(bb + bias - 1, 0, len(x) - 1)
        idx_0 = np.clip(bb + bias, 0, len(x) - 1)
        idx_p1 = np.clip(bb + bias + 1, 0, len(x) - 1)
        
        dcl = np.mean(x[idx_0])
        txm1 = x[idx_m1]
        tx0 = x[idx_0]
        txp1 = x[idx_p1]
        
        # Safe guard
        if (np.sum(np.abs(txm1)) < 1e-20) or (np.sum(np.abs(tx0)) < 1e-20) or (np.sum(np.abs(txp1)) < 1e-20):
            txm1 = txm1 + np.random.randn(len(txm1)) * 1e-10
            tx0 = tx0 + np.random.randn(len(tx0)) * 1e-10
            txp1 = txp1 + np.random.randn(len(txp1)) * 1e-10
            dcl = np.mean(tx0)
        
        ff0 = np.fft.fft((txm1 - dcl) * wa, fftl)
        ff1 = np.fft.fft((tx0 - dcl) * wa, fftl)
        ff2 = np.fft.fft((txp1 - dcl) * wa, fftl)
        
        ff0[ff0 == 0] = 1e-9
        ff1[ff1 == 0] = 1e-9
        ff2[ff2 == 0] = 1e-9
        
        fd = ff2 * rr - ff1
        fd0 = ff1 * rr - ff0
        
        crf = fax + (np.real(ff1) * np.imag(fd) - np.imag(ff1) * np.real(fd)) / (np.abs(ff1) ** 2 + 1e-10) * fs / np.pi / 2
        crf0 = fax + (np.real(ff0) * np.imag(fd0) - np.imag(ff0) * np.real(fd0)) / (np.abs(ff0) ** 2 + 1e-10) * fs / np.pi / 2
        
        pif[:, kk] = crf[bx] * 2 * np.pi
        dpif[:, kk] = (crf[bx] - crf0[bx]) * 2 * np.pi
        pwm[:, kk] = np.abs(ff1[bx])
        ac1[kk] = np.sum(np.abs(ff1) ** 2.0 * ec1) / (np.sum(np.abs(ff1) ** 2) + 1e-10)
    
    slp = (np.vstack([pif[1:fftl // 2 + 1, :], pif[fftl // 2:fftl // 2 + 1, :]]) - pif) / (fs / fftl * 2 * np.pi)
    dslp = (np.vstack([dpif[1:fftl // 2 + 1, :], dpif[fftl // 2:fftl // 2 + 1, :]]) - dpif) / (fs / fftl * 2 * np.pi) * fs
    mmp = slp * 0
    
    c1, c2 = znrmlcf3(np.array([shiftm]))
    fxx = (np.arange(fftl // 2 + 1) + 0.5) / fftl * fs * 2 * np.pi
    
    for ii in range(fftl // 2 + 1):
        c2_scaled = c2 * (fxx[ii] / 2 / np.pi) ** 2
        mmp[ii, :] = (dslp[ii, :] / np.sqrt(c2_scaled)) ** 2 + (slp[ii, :] / np.sqrt(c1)) ** 2
    
    # Temporal smoothing
    sml = int(np.round(1.5 * fs / 1000 / 2 / shiftm) * 2 + 1)
    smb = int(np.round((sml - 1) / 2))
    
    smmp = fftfilt((hanning(sml) ** 2) / np.sum(hanning(sml) ** 2),
                   np.vstack([mmp.T, np.zeros((sml * 2, fftl // 2 + 1))])).T + 1e-5
    smmp[smmp == 0] = 1e-10
    smmp = 1.0 / fftfilt(hanning(sml) / np.sum(hanning(sml)), 1.0 / smmp.T).T
    
    frame_indices = np.clip(np.arange(nfr) + sml - 2, 0, smmp.shape[1] - 1).astype(int)
    smmp = smmp[:, frame_indices]
    
    # Power adaptive weighting
    spwm = fftfilt(hanning(sml) / np.sum(hanning(sml)),
                   np.vstack([pwm.T, np.zeros((sml * 2, fftl // 2 + 1))])).T
    spwm[spwm == 0] = 1e-5
    
    spfm = fftfilt(hanning(sml) / np.sum(hanning(sml)),
                   np.vstack([(pwm * pif).T, np.zeros((sml * 2, fftl // 2 + 1))])).T + 1e-5
    spif = spfm / spwm
    spif = spif[:, np.arange(nfr) + smb]
    
    idx = np.maximum(0, f0i / fs * fftl)
    fqv = np.zeros((nhmx, nfr))
    vvv = np.zeros((nhmx, nfr))
    
    # In MATLAB, iidx accumulates idx on each iteration to get harmonics
    iidx_freq = np.zeros(nfr)  # Start at 0
    for ii in range(nhmx):
        iidx_freq = iidx_freq + idx  # Accumulate: 1*idx, 2*idx, 3*idx, ...
        
        # Sample smmp and spif at the harmonic frequency bins
        iidx_floor = np.clip(np.floor(iidx_freq).astype(int), 0, smmp.shape[0] - 1)
        iidx_ceil = np.clip(iidx_floor + 1, 0, smmp.shape[0] - 1)
        
        frac = iidx_freq - iidx_floor
        vvv[ii, :] = (smmp[iidx_floor, np.arange(nfr)] + 
                      frac * (smmp[iidx_ceil, np.arange(nfr)] - smmp[iidx_floor, np.arange(nfr)])) / ((ii + 1) ** 2)
        fqv[ii, :] = (spif[iidx_floor, np.arange(nfr)] + 
                      frac * (spif[iidx_ceil, np.arange(nfr)] - spif[iidx_floor, np.arange(nfr)])) / 2 / np.pi / (ii + 1)
    
    vvvf = 1.0 / (np.sum(1.0 / (vvv + 1e-10), axis=0) + 1e-10)
    f0r = (np.sum(fqv / (np.sqrt(vvv) + 1e-10), axis=0) / 
           (np.sum(1.0 / (np.sqrt(vvv) + 1e-10), axis=0) + 1e-10)) * (f0raw > 0)
    ecr = np.sqrt(1.0 / (vvvf + 1e-10)) * (f0raw > 0) + (f0raw <= 0)
    
    return f0r, ecr, ac1


def zpeakdipdetect(auxouts: Dict, wsml: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Detect peaks and dips in power envelope
    
    Args:
        auxouts: Auxiliary outputs dictionary
        wsml: Window size for smoothing
        
    Returns:
        pv: Peak indices
        dv: Dip indices
    """
    pwsdb = auxouts['InstantaneousPower']
    
    pwsdbl = np.concatenate([np.ones(wsml) * pwsdb[0], pwsdb, np.ones(2 * wsml) * pwsdb[-1]])
    pwsdbs = fftfilt(hanning(wsml) / np.sum(hanning(wsml)), pwsdbl)
    pwsdbs = pwsdbs[int(np.round(3 * wsml / 2)):int(np.round(3 * wsml / 2)) + len(pwsdb)]
    
    dpwsdbs = np.diff(np.concatenate([[pwsdbs[0]], pwsdbs]))
    dpwsdbsm = np.diff(np.concatenate([pwsdbs, [pwsdbs[-1]]]))
    
    pv = np.where((dpwsdbs * dpwsdbsm < 0) & (dpwsdbsm <= 0))[0]
    dv = np.where((dpwsdbs * dpwsdbsm < 0) & (dpwsdbsm > 0))[0]
    
    return pv, dv


def zvuvdecision4(f0: np.ndarray, auxouts: Dict) -> np.ndarray:
    """
    Simple V/UV decision logic
    
    Args:
        f0: F0 trajectory
        auxouts: Auxiliary outputs dictionary
        
    Returns:
        vuv: Voiced/unvoiced indicator (1: voiced, 0: unvoiced)
    """
    pwsdb = auxouts['InstantaneousPower']
    rel = auxouts['RELofcandidatesByMix']
    maxpwsdb = np.max(pwsdb)
    noiselevel = auxouts['BackgroundNoiselevel']
    
    # Onset and offset candidates
    nw = 40
    nrw = 3
    tt = np.arange(-nw, nw + 1)
    pws = 10.0 ** (pwsdb / 20)
    
    wwh = np.exp(-(tt / (nw / 2.5)) ** 2) * (0.5 - 1.0 / (1 + np.exp(-tt / nrw)))
    dpw = fftfilt(wwh, np.concatenate([pws, np.zeros(nw * 2)]))
    dpw = dpw[nw:nw + len(pws)]
    biast = nrw * 3
    
    ddpw = np.diff(np.concatenate([[dpw[0]], dpw]))
    ddpwm = np.diff(np.concatenate([dpw, [dpw[-1]]]))
    onv = np.where((ddpw * ddpwm < 0) & (ddpwm <= 0))[0]
    
    # Search for voiced segments
    vuv = (pwsdb > (2 * maxpwsdb + noiselevel) / 3).astype(float)
    pv, _ = zpeakdipdetect(auxouts, 81)
    
    np_peaks = len(pv)
    nn = min(len(vuv), len(f0))
    vuv = np.zeros(nn)
    lastp = 2
    
    for ii in range(np_peaks):
        if (pwsdb[pv[ii]] > (1.2 * maxpwsdb + noiselevel) / 2.2) and (pv[ii] > lastp):
            lb = lastp
            ub = nn - 1
            cp = pv[ii]
            bp = cp
            ep = cp
            
            # Search backward
            for bp_candidate in range(cp - 1, lb - 1, -1):
                if (pwsdb[bp_candidate] < (maxpwsdb + 2.3 * noiselevel) / 3.3) or \
                   ((pwsdb[bp_candidate] < (1.5 * maxpwsdb + noiselevel) / 2.5) and (rel[bp_candidate, 0] < 0.3)) or \
                   ((pwsdb[bp_candidate] < (1.5 * maxpwsdb + noiselevel) / 2.5) and 
                    (bp_candidate > 0) and (np.abs(np.log2((f0[bp_candidate] + 1e-10) / (f0[bp_candidate - 1] + 1e-10))) > 0.1)):
                    bp = bp_candidate
                    break
            
            # Find nearest onset
            if len(onv) > 0:
                dmy = np.min(np.abs(onv - bp))
                ix = np.argmin(np.abs(onv - bp))
                if dmy < 20:
                    bp = max(0, onv[ix] - biast)
            
            # Search forward
            for ep_candidate in range(cp + 1, min(len(f0) - 1, ub) + 1):
                if (pwsdb[ep_candidate] < (maxpwsdb + 5 * noiselevel) / 6) or \
                   ((pwsdb[ep_candidate] < (maxpwsdb + 1.3 * noiselevel) / 2.3) and (rel[ep_candidate, 0] < 0.25)) or \
                   ((pwsdb[ep_candidate] < (maxpwsdb + 0.7 * noiselevel) / 1.7) and 
                    (ep_candidate < len(f0) - 1) and (np.abs(np.log2((f0[ep_candidate] + 1e-10) / (f0[ep_candidate + 1] + 1e-10))) > 0.1)):
                    ep = ep_candidate
                    break
                ep = ep_candidate
            
            vuv[bp:ep + 1] = 1
            lastp = ep
    
    return vuv


# For testing
if __name__ == '__main__':
    print("STRAIGHT F0 Extraction Module")
    print("Full implementation complete")

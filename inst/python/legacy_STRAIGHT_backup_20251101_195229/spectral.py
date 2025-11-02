"""
STRAIGHT Spectral Analysis Module

Python implementation of exstraightspec, straightBodyC03ma, and related functions.
Faithful port from MATLAB preserving algorithms and numerical behavior.

Original MATLAB version:
    Copyright (c) Wakayama University, 2004-2006
    Designed and coded by Hideki Kawahara
"""

import numpy as np
from scipy.signal import butter, lfilter
from scipy.signal.windows import hann
from dataclasses import dataclass
from typing import Tuple, Optional
import warnings

# Import helpers from synthesis
from .synthesis import hanning, fftfilt


@dataclass
class SpectralParameters:
    """
    Parameters for STRAIGHT spectral analysis
    
    Attributes match MATLAB field names for compatibility
    """
    DisplayPlots: int = 0  # 0: no graphics, 1: show graphics
    defaultFrameLength: float = 80.0  # frame length (ms)
    spectralUpdateInterval: float = 1.0  # frame shift (ms)
    spectralTimeWindowStretch: float = 1.0  # time window stretch factor (eta)
    spectralExponentForNonlinearity: float = 0.6  # exponent for nonlinearity (pc)
    spectralTimeDomainCompensation: float = 0.2  # TD compensation (mag)


def exstraightspec(
    x: np.ndarray,
    f0raw: np.ndarray,
    fs: float,
    optionalParamsSP: Optional[SpectralParameters] = None
) -> Tuple[np.ndarray, SpectralParameters]:
    """
    Spectral information extraction for STRAIGHT
    
    Faithful port of exstraightspec.m
    
    Args:
        x: Input signal (only first channel analyzed)
        f0raw: Fundamental frequency (Hz) in 1 ms temporal resolution, 0 for aperiodic
        fs: Sampling frequency (Hz)
        optionalParamsSP: Optional spectrum analysis parameters
        
    Returns:
        n3sgram: Smoothed time-frequency representation (spectrogram)
        analysisParamsSp: Actually used parameters
        
    Example:
        >>> n3sgram, params = exstraightspec(x, f0raw, fs)
    """
    # Initialize parameters
    if optionalParamsSP is None:
        prm = SpectralParameters()
    else:
        prm = optionalParamsSP
    
    # Extract parameters
    imageOn = prm.DisplayPlots
    framem = prm.defaultFrameLength
    shiftm = prm.spectralUpdateInterval
    eta = prm.spectralTimeWindowStretch
    pc = prm.spectralExponentForNonlinearity
    mag = prm.spectralTimeDomainCompensation
    
    framel = int(framem * fs / 1000)
    fftl = 1024  # default FFT length
    
    if fftl < framel:
        fftl = 2 ** int(np.ceil(np.log(framel) / np.log(2)))
    
    # Ensure x is a column vector
    if x.ndim > 1:
        nr, nc = x.shape
        if nr > nc:
            xold = x[:, 0]
        else:
            xold = x[0, :]
    else:
        xold = x
    
    xold = xold.flatten()
    
    # Spectral estimation with normalization
    xamp = np.std(xold)
    scaleconst = 2200  # magic number for compatibility
    xold = xold / xamp * scaleconst
    
    # Obsolete dummy variables (kept for compatibility)
    f0var = 1
    f0varL = 1
    
    # Main spectral analysis
    n2sgrambk, _ = straightBodyC03ma(
        xold, fs, shiftm, fftl, f0raw, f0var, f0varL, eta, pc, imageOn
    )
    
    # Time domain compensation
    if mag > 0:
        n3sgram = specreshape(fs, n2sgrambk, eta, pc, mag, f0raw, imageOn)
    else:
        n3sgram = n2sgrambk
    
    # Denormalize
    n3sgram = n3sgram / scaleconst * xamp
    
    analysisParamsSp = prm
    
    return n3sgram, analysisParamsSp


def straightBodyC03ma(
    x: np.ndarray,
    fs: float,
    shiftm: float,
    fftl: int,
    f0raw: np.ndarray,
    f0var: float,
    f0varL: float,
    eta: float,
    pc: float,
    imgi: int
) -> Tuple[np.ndarray, np.ndarray]:
    """
    STRAIGHT body: Interpolation using adaptive Gaussian weighting
    
    Faithful port of straightBodyC03ma.m
    
    Args:
        x: Input waveform
        fs: Sampling frequency (Hz)
        shiftm: Frame shift (ms)
        fftl: Length of FFT
        f0raw: Pitch information to guide analysis
        f0var: Expected f0 variance including zerocross information
        f0varL: Expected f0 variance
        eta: Time window stretch factor
        pc: Exponent for nonlinearity
        imgi: Display indicator (1: on, 0: off)
        
    Returns:
        n2sgram: Smoothed spectrogram
        nsgram: Isometric spectrogram
    """
    # Convert fs to float to avoid uint16 overflow issues
    fs = float(fs)
    
    # Dummy variable usage (kept for compatibility)
    f0l = f0raw.flatten() + 0 * f0var + 0 * f0varL
    
    framem = 80
    framel = int(np.round(framem * fs / 1000))
    
    if fftl < framel:
        warnings.warn('Warning! fftl is too small.')
        fftl = 2 ** int(np.ceil(np.log(framel) / np.log(2)))
        warnings.warn(f'New length: {fftl} is used.')
    
    x = x.flatten()
    shiftl = shiftm * fs / 1000
    
    # High-pass filtering using 70Hz cut-off butterworth filter
    # Note: MATLAB uses filter() (causal), not filtfilt() (zero-phase)
    b, a = butter(6, 70 / (fs / 2), 'high')
    from scipy.signal import lfilter
    xh = lfilter(b, a, x)
    rmsp = np.std(xh)
    
    b, a = butter(6, 300 / (fs / 2), 'high')
    xh2 = lfilter(b, a, x)
    
    # High-pass filter using 3000Hz cut-off butterworth filter
    b, a = butter(6, 3000 / (fs / 2), 'high')
    xhh = lfilter(b, a, x)
    
    # Add noise padding
    tx = np.concatenate([
        np.random.randn(int(np.round(framel / 2))) * rmsp / 4000,
        xh,
        np.random.randn(framel) * rmsp / 4000
    ])
    
    nframe = min(len(f0l), int(np.round(len(x) / shiftl)))
    
    nsgram = np.zeros((fftl // 2 + 1, nframe))
    n2sgram = np.zeros((fftl // 2 + 1, nframe))
    
    tt = (np.arange(framel) - framel / 2) / fs
    bbase = np.arange(fftl // 2 + 1)
    
    ist = 0  # Python 0-based indexing
    f0x = np.zeros_like(f0l)
    
    # Optimum blending table for interference free spec
    cfv = np.array([0.36, 0.30, 0.26, 0.21, 0.17, 0.14, 0.10])
    muv = np.array([1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6])
    
    bcf = np.interp(eta, muv, cfv)
    
    # Calculate optimum smoothing function coefficients
    ovc = optimumsmoothing(eta, pc)
    
    # Designing pitch synchronized gaussian
    fNominal = 40
    wGaussian = np.exp(-np.pi * (tt * fNominal / eta) ** 2)
    wSynchronousBartlett = np.maximum(0, 1 - np.abs(tt * fNominal))
    wPSGSeed = fftfilt(
        wSynchronousBartlett[wSynchronousBartlett > 0],
        np.concatenate([wGaussian, np.zeros(len(tt))])
    )
    wPSGSeed = wPSGSeed / np.max(wPSGSeed)
    maxLocation = np.argmax(wPSGSeed)
    tNominal = (np.arange(len(wPSGSeed)) - maxLocation) / fs
    
    ttm = np.concatenate([[0.00001], np.arange(1, fftl // 2 + 1), np.arange(-fftl // 2 + 1, 0)]) / fs
    
    # Safeguard
    lft = 1.0 / (1 + np.exp(-(np.abs(np.arange(fftl) - fftl / 2) - fftl / 30) / 2))
    
    # Main analysis loop - F0 adaptive time-frequency analysis
    for ii in range(nframe):
        f0 = f0l[max(0, ii)]
        if f0 == 0:
            f0 = 160
        
        f0x[ii] = f0
        t0 = 1 / f0
        
        # Interpolate pitch-synchronized window
        wxe = np.interp(tt * f0 / fNominal, tNominal, wPSGSeed, left=0, right=0)
        wxe[np.isnan(wxe)] = 0
        wxe = wxe / np.sqrt(np.sum(wxe ** 2))
        wxd = bcf * wxe * np.sin(np.pi * tt / t0)
        
        # Extract frame
        # MATLAB: iix=round(ist:ist+framel-1)
        # This extracts exactly framel samples
        iix_start = int(np.round(ist))
        iix_end = iix_start + framel
        iix = np.arange(iix_start, iix_end)
        
        # Ensure we don't go out of bounds
        if iix_end > len(tx):
            # Pad tx if needed (shouldn't happen with proper sizing)
            pad_needed = iix_end - len(tx)
            tx = np.pad(tx, (0, pad_needed), mode='constant')
        
        tx_frame = tx[iix]
        tx_frame = tx_frame - np.mean(tx_frame)
        
        # Compute power spectrum
        pw = np.sqrt(
            np.abs(np.fft.fft(tx_frame * wxe, fftl)) ** 2 +
            np.abs(np.fft.fft(tx_frame * wxd, fftl)) ** 2
        ) ** pc
        
        nsgram[:, ii] = pw[bbase]
        
        # F0-adaptive smoothing
        f0p2 = int(np.floor((f0 / fs * fftl) / 2 + 1))
        f0p = int(np.ceil((f0 / fs * fftl) + 1))
        f0pr = f0 / fs * fftl + 1
        
        # MATLAB: tmppw=interp1(1:f0p,pw(1:f0p),f0pr-((1:f0p2)-1))
        # In Python 0-based: interpolate pw[0:f0p] at positions f0pr - (0:f0p2-1) - 1
        indices_matlab = f0pr - np.arange(f0p2)  # MATLAB (1:f0p2) becomes 0:f0p2 in Python
        x_vals = np.arange(1, f0p + 1)  # MATLAB 1:f0p
        indices_matlab = np.clip(indices_matlab, 1, f0p)
        tmppw = np.interp(indices_matlab, x_vals, pw[:f0p])
        pw[:f0p2] = tmppw
        # MATLAB: pw(fftl:-1:fftl-f0p2+2)=pw(2:f0p2)
        # Python: pw[fftl-1:fftl-f0p2:-1] = pw[1:f0p2]
        if fftl - f0p2 + 1 < fftl:
            pw[fftl-1:fftl-f0p2:-1] = pw[1:f0p2]
        
        # Local level equalization
        ww2t = (np.sin(ttm / (t0 / 3) * np.pi) / (ttm / (t0 / 3) * np.pi + 1e-10)) ** 2
        spw2 = np.real(np.fft.ifft(ww2t * np.fft.fft(pw) * lft))
        spw2[spw2 == 0] = np.finfo(float).eps
        
        # Optimum weighting
        wwt = ((np.sin(ttm / t0 * np.pi) / (ttm / t0 * np.pi + 1e-10)) ** 2 *
               (ovc[0] + ovc[1] * 2 * np.cos(ttm / t0 * 2 * np.pi) +
                ovc[2] * 2 * np.cos(ttm / (t0 / 2) * 2 * np.pi)))
        spw = np.real(np.fft.ifft(wwt * np.fft.fft(pw / spw2))) / wwt[0]
        
        # Smooth half wave rectification
        n2sgram[:, ii] = (spw2[bbase] *
                          (0.25 * (np.log(2 * np.cosh(spw[bbase] * 4 / 1.4)) * 1.4 +
                                   spw[bbase] * 4) / 2))
        
        ist = ist + shiftl
    
    nsgram = nsgram ** (1 / pc)
    n2sgram = n2sgram ** (2 / pc)
    
    # Time constant control for unvoiced parts
    ttlv = np.sum(n2sgram)
    ncw = int(np.round(2 * fs / 1000))
    
    lbb = int(np.round(300 / fs * fftl))
    
    h3 = np.convolve(hanning(int(np.round(fs / 1000))),
                     np.exp(-1400 / fs * np.arange(ncw * 2 + 1)))
    pwc = fftfilt(h3, np.abs(np.concatenate([xh2, np.zeros(ncw * 10)])) ** 2)
    # MATLAB: pwc(round(1:fs/(1000/shiftm):length(pwc)))
    # Start at index 0 in Python (MATLAB index 1), then step
    step = int(np.round(fs / (1000 / shiftm)))
    indices = np.arange(0, len(pwc), step)  # Start at 0 (MATLAB's 1)
    pwc = pwc[indices]
    
    nn, mm = n2sgram.shape
    pwc = pwc[:mm]
    pwc = pwc / np.sum(pwc) * np.sum(n2sgram[lbb:nn, :])
    
    pwch = fftfilt(h3, np.abs(np.concatenate([xhh, np.zeros(ncw * 10)])) ** 2)
    pwch = pwch[indices[:len(pwch)]] if len(indices) <= len(pwch) else pwch[indices]
    pwch = pwch[:mm]
    pwch = pwch / np.sum(pwch) * ttlv
    
    # Impact detection
    ipwm = 7
    ipl = int(np.round(ipwm / shiftm))
    ww = hanning(ipl * 2 + 1)
    ww = ww / np.sum(ww)
    
    apwt = fftfilt(ww, np.concatenate([pwch, np.zeros(len(ww) * 2)]))
    apwt = apwt[ipl:ipl + len(pwch)]
    
    dpwt = fftfilt(ww, np.concatenate([np.diff(pwch) ** 2, np.zeros(len(ww) * 2)]))
    dpwt = dpwt[ipl:ipl + len(pwch)]
    
    mmaa = np.max(apwt)
    apwt[apwt <= 0] = mmaa
    
    rr = np.sqrt(dpwt) / apwt
    lmbd = 1.0 / (1 + np.exp(-(np.sqrt(rr) - 0.75) * 20))
    
    pwc = pwc * lmbd + (1 - lmbd) * np.sum(n2sgram, axis=0)
    
    # Shaping amplitude envelope
    for ii in range(mm):
        if f0raw[ii] == 0:
            n2sgram[:, ii] = pwc[ii] * n2sgram[:, ii] / np.sum(n2sgram[:, ii])
    
    n2sgram = np.abs(n2sgram + 1e-10)
    n2sgram = np.sqrt(n2sgram)
    
    return n2sgram, nsgram


def specreshape(
    fs: float,
    n2sgram: np.ndarray,
    eta: float,
    pc: float,
    mag: float,
    f0: np.ndarray,
    imgi: int
) -> np.ndarray:
    """
    Spectral compensation using Time Domain technique
    
    Port of specreshape.m
    
    Args:
        fs: Sampling frequency (Hz)
        n2sgram: STRAIGHT smoothed spectrogram
        eta: Temporal stretch factor
        pc: Power exponent for nonlinearity
        mag: Magnification factor of TD compensation
        f0: Fundamental frequency (Hz)
        imgi: Display indicator
        
    Returns:
        n2sgram3: Compensated spectrogram
    """
    nn, mm = n2sgram.shape
    fftl = (nn - 1) * 2
    fbb = np.arange(nn)
    rbb = np.arange(nn - 2, 0, -1)
    rbb2 = np.arange(fftl - 1, nn - 1, -1)
    bb3 = np.arange(1, nn - 1)
    n2sgram3 = np.zeros_like(n2sgram)
    
    ovc = optimumsmoothing(eta, pc)
    hh = np.array([[1, 1, 1, 1],
                   [0, 1/2, 2/3, 3/4],
                   [0, 0, 1/3, 2/4],
                   [0, 0, 0, 1/4]])
    bb = np.linalg.solve(hh, ovc)
    
    tt = np.arange(fftl) / fs
    pb2 = (np.pi / (eta ** 2) +
           (np.pi ** 2) / 3 * (bb[0] + 4 * bb[1] + 9 * bb[2] + 16 * bb[3])) * tt ** 2
    
    for ii in range(mm):
        ffs = np.concatenate([n2sgram[:, ii], n2sgram[rbb, ii]])
        ccs2 = np.real(np.fft.fft(ffs)) * np.minimum(20, 1 + mag * pb2 * f0[ii] ** 2)
        ccs2[rbb2] = ccs2[bb3]
        ngg = np.real(np.fft.ifft(ccs2))
        n2sgram3[:, ii] = ngg[fbb]
    
    n2sgram3 = (np.abs(n2sgram3) + n2sgram3) / 2 + 0.1
    
    return n2sgram3


def optimumsmoothing(eta: float, pc: float) -> np.ndarray:
    """
    Calculate optimum smoothing function
    
    Port of optimumsmoothing.m
    
    Args:
        eta: Temporal stretch factor
        pc: Power exponent for nonlinearity
        
    Returns:
        ovc: Coefficients for 2nd order cardinal B-spline
    """
    fx = np.arange(-8, 8.05, 0.05)
    cb = np.maximum(0, 1 - np.abs(fx))
    gw = np.exp(-np.pi * (fx * eta * 1.4) ** 2) ** pc
    cmw = np.convolve(cb, gw, mode='full')
    
    bb = np.arange(len(cb))
    bbc = bb + (len(cb) - 1) // 2
    bbc = np.clip(bbc, 0, len(cmw) - 1).astype(int)
    cmw = cmw[bbc] / np.max(cmw)
    
    # Find samples near integers
    ss = np.where(np.abs(fx - np.round(fx)) < 0.025)[0]
    cmws = cmw[ss]
    
    nn = len(cmws)
    idv = np.arange(nn)
    
    hh = np.zeros((2 * nn, nn))
    for ii in range(nn):
        indices = (ii) + idv
        indices = indices[indices < 2 * nn]
        hh[indices, ii] = cmws[:len(indices)]
    
    bv = np.zeros(2 * nn)
    bv[nn] = 1  # Unit impulse
    
    h = hh.T @ hh
    ov = np.linalg.solve(h, hh.T @ bv)
    
    idc = (nn - 1) // 2 + 1
    ovc = ov[idc:idc + 4]
    
    return ovc


# For testing
if __name__ == '__main__':
    print("STRAIGHT Spectral Analysis Module")
    print("Use: from straight.spectral import exstraightspec")

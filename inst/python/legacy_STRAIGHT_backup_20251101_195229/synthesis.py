"""
STRAIGHT Synthesis Module

Python implementation of exstraightsynth and straightSynthTB07ca.
Faithful port from MATLAB preserving algorithms and numerical behavior.

Original MATLAB version:
    Copyright (c) ATR Human Info. Proc. Res. Labs. 1996-2016
    Designed and coded by Hideki Kawahara
"""

import numpy as np
from scipy.signal.windows import hann
from dataclasses import dataclass, field
from typing import Tuple, Optional, Union
import warnings


def hanning(N):
    """
    Wrapper for scipy's hann function to match MATLAB behavior
    
    MATLAB's hanning(N) is equivalent to hann(N+2, sym=True)[1:-1]
    This ensures non-zero values at endpoints.
    """
    return hann(N + 2, sym=True)[1:-1]


def fftfilt(b, x):
    """
    FFT-based FIR filtering (like MATLAB's fftfilt)
    
    Args:
        b: Filter coefficients (1D)
        x: Input signal (1D or 2D)
        
    Returns:
        Filtered signal (same shape as x)
    """
    b = np.asarray(b).flatten()
    x = np.asarray(x)
    
    # Handle 2D arrays by filtering each column
    if x.ndim == 2:
        y = np.zeros_like(x)
        for i in range(x.shape[1]):
            y[:, i] = fftfilt(b, x[:, i])
        return y
    
    # 1D case
    x = x.flatten()
    
    # Determine FFT length
    n_min = len(b) + len(x) - 1
    n_fft = 2 ** int(np.ceil(np.log2(n_min)))
    
    # FFT-based convolution
    B = np.fft.fft(b, n_fft)
    X = np.fft.fft(x, n_fft)
    Y = B * X
    y = np.real(np.fft.ifft(Y))
    
    # Return same length as input (like filter)
    return y[:len(x)]


@dataclass
class SynthesisParameters:
    """
    Parameters for STRAIGHT synthesis
    
    Attributes match MATLAB field names for compatibility
    """
    spectralUpdateInterval: float = 1.0  # frame shift (ms)
    groupDelayStandardDeviation: float = 0.5  # std dev of random group delay (ms)
    groupDelaySpatialBandWidth: float = 70.0  # smoothing window length (Hz)
    groupDelayRandomizeCornerFrequency: float = 4000.0  # corner freq for random phase (Hz)
    ratioToFundamentalPeriod: float = 0.2  # fractional group delay (ratio)
    ratioModeIndicator: int = 0  # use fractional group delay if 1
    levelNormalizationIndicator: int = 1  # normalize voiced part level if 1
    headRoomToClip: float = 22.0  # head room from voiced RMS to clipping (dB)
    powerCheckSegmentLength: float = 15.0  # segment length for power check (ms)
    timeAxisMappingTable: Union[int, np.ndarray] = 1  # identity mapping = 1
    fundamentalFrequencyMappingTable: Union[int, np.ndarray] = 1  # identity = 1
    frequencyAxisMappingTable: Union[int, np.ndarray] = 1  # identity = 1
    timeAxisStretchingFactor: float = 1.0  # duration stretch coefficient
    DisplayPlots: int = 0  # 0: no display, 1: display on
    lowestF0: float = 50.0  # lower limit of resynthesized F0 (Hz)
    statusReport: str = 'ok'  # status message


def exstraightsynth(
    f0raw: np.ndarray,
    n3sgram: np.ndarray,
    ap: np.ndarray,
    fs: float,
    optionalParamsS: Optional[SynthesisParameters] = None
) -> Tuple[np.ndarray, SynthesisParameters]:
    """
    Synthesis using STRAIGHT parameters with linear modifications
    
    This is a faithful Python port of exstraightsynth.m
    
    Args:
        f0raw: Fundamental frequency (Hz), shape (n_frames,)
        n3sgram: STRAIGHT spectrogram (absolute value), shape (n_freq, n_frames)
        ap: Aperiodic component (dB re. to total power), shape (n_freq, n_frames)
        fs: Sampling frequency (Hz)
        optionalParamsS: Optional synthesis parameters
        
    Returns:
        sy: Synthesized speech signal
        prmS: Actually used synthesis parameters
        
    Example:
        >>> sy, params = exstraightsynth(f0raw, n3sgram, ap, fs)
    """
    # Initialize parameters
    if optionalParamsS is None:
        prmS = SynthesisParameters()
    else:
        prmS = optionalParamsS
    
    # Extract parameters
    shiftm = prmS.spectralUpdateInterval
    delsp = prmS.groupDelayStandardDeviation
    gdbw = prmS.groupDelaySpatialBandWidth
    cornf = prmS.groupDelayRandomizeCornerFrequency
    delfrac = prmS.ratioToFundamentalPeriod
    delfracind = prmS.ratioModeIndicator
    normalizedOut = prmS.levelNormalizationIndicator
    headRoom = prmS.headRoomToClip
    lsegment = prmS.powerCheckSegmentLength
    imap = prmS.timeAxisMappingTable
    pconv = prmS.fundamentalFrequencyMappingTable
    fconv = prmS.frequencyAxisMappingTable
    sconv = prmS.timeAxisStretchingFactor
    imgi = prmS.DisplayPlots
    lowestF0 = prmS.lowestF0
    
    # Call synthesis engine
    sy, statusReport = straightSynthTB07ca(
        n3sgram, f0raw, shiftm, fs,
        pconv, fconv, sconv, gdbw, delfrac, delsp, cornf, 
        delfracind, ap, imap, imgi, lowestF0
    )
    
    # Normalize output if requested
    if normalizedOut:
        dBsy = powerchk(sy, fs, lsegment)
        cf = (20 * np.log10(32768) - headRoom) - dBsy
        sy = sy * (10.0 ** (cf / 20))
    
    prmS.statusReport = statusReport
    
    return sy, prmS


def straightSynthTB07ca(
    n2sgram: np.ndarray,
    f0raw: np.ndarray,
    shiftm: float,
    fs: float,
    pcnv: Union[float, np.ndarray],
    fconv: Union[float, np.ndarray],
    sconv: float,
    gdbw: float,
    delfrac: float,
    delsp: float,
    cornf: float,
    delfracind: int,
    ap: np.ndarray,
    imap: Union[int, np.ndarray],
    imgi: int,
    lowestF0: float
) -> Tuple[np.ndarray, str]:
    """
    STRAIGHT synthesis with all-pass filter design
    
    Faithful port of straightSynthTB07ca.m
    
    Args:
        n2sgram: Amplitude spectrogram, shape (n_freq, n_frames)
        f0raw: Pitch pattern (Hz), shape (n_frames,)
        shiftm: Frame shift (ms) for spectrogram
        fs: Sampling frequency (Hz)
        pcnv: Pitch stretch factor
        fconv: Frequency stretch factor
        sconv: Speaking duration stretch factor
        gdbw: Finest resolution in group delay (Hz)
        delfrac: Ratio of std dev of group delay in terms of F0
        delsp: Standard deviation of group delay (ms)
        cornf: Lower corner frequency for phase randomization (Hz)
        delfracind: Selector of fixed and proportional group delay
        ap: Aperiodicity measure, shape (n_freq, n_frames)
        imap: Arbitrary mapping from new time (sample) to old time (frame)
        imgi: Display indicator, 1: display on, 0: off
        lowestF0: Lower limit of the resynthesized fundamental frequency (Hz)
        
    Returns:
        sy: Synthesized speech signal
        statusReport: Status message
    """
    statusReport = 'ok'
    
    # Ensure f0raw is 1D
    f0l = f0raw.flatten() if f0raw.ndim > 1 else f0raw.copy()
    
    # Get dimensions
    nii, njj = n2sgram.shape
    njj = min(njj, len(f0raw))
    f0l = f0l[:njj]
    
    # Check minimum F0
    f0_voiced = f0l[f0l > 0]
    if len(f0_voiced) > 0 and np.min(f0_voiced) * pcnv < lowestF0:
        statusReport = f'Minimum synthesized F0 exceeded the lower limit({lowestF0} Hz).'
    
    # Determine FFT length
    fftLengthForLowestF0 = 2 ** int(np.ceil(np.log2(2 * np.round(fs / lowestF0))))
    fftl = nii + nii - 2
    
    if fftl < fftLengthForLowestF0:
        niiNew = fftLengthForLowestF0 // 2 + 1
        statusReport = 'The FFT length was inconsistent and replaced'
        
        # Interpolate spectrogram and aperiodicity
        freq_old = np.arange(nii)
        freq_new = np.arange(niiNew) * (nii - 1) / (niiNew - 1)
        n2sgram = np.array([np.interp(freq_new, freq_old, n2sgram[:, i]) 
                            for i in range(n2sgram.shape[1])]).T
        ap = np.array([np.interp(freq_new, freq_old, ap[:, i]) 
                       for i in range(ap.shape[1])]).T
        
        fftl = fftLengthForLowestF0
        nii = niiNew
    
    # Safeguard for ap mismatch
    if ap.shape[0] != n2sgram.shape[0]:
        apDouble = np.zeros_like(n2sgram)
        freq_old_ap = np.arange(ap.shape[0])
        freq_new_ap = np.arange(n2sgram.shape[0]) * (ap.shape[0] - 1) / (n2sgram.shape[0] - 1)
        
        for ik in range(ap.shape[1]):
            apDouble[:, ik] = np.interp(freq_new_ap, freq_old_ap, ap[:, ik], 
                                        left=ap[0, ik], right=ap[-1, ik])
        ap = apDouble
    
    # Aperiodicity conversion
    aprms = 10.0 ** (ap / 20)
    aprm = np.clip(aprms * 1.6 - 0.015, 0.001, 1.0)
    
    # Frequency axis mapping
    if isinstance(fconv, (int, float)):
        idcv = np.minimum(np.arange(fftl // 2 + 1) / fconv, fftl // 2)
    elif len(fconv) == nii:
        idcv = fconv
    else:
        idcv = np.arange(fftl // 2 + 1)
        statusReport += '\nFrequency axis mapping function is not consistent with lowestF0.'
    
    # Time axis mapping
    if isinstance(imap, int) or len(imap) <= 1:
        sy = np.zeros(int(np.round((njj * shiftm / 1000 * fs) * sconv + 3 * fftl + 1)))
        imap = np.arange(len(sy))
        imap = np.minimum(((imap) / fs * 1000 / shiftm / sconv), len(f0l) - 1)
    else:
        sy = np.zeros(len(imap) + 3 * fftl)
    
    # Safeguard
    imap = np.concatenate([imap, np.ones(int(np.round(fs * 0.2))) * (len(f0l) - 1)])
    ix = np.where(imap >= len(f0l) - 1)[0]
    if len(ix) > 0:
        ix = ix[0]
    else:
        ix = len(imap)
    
    rmap = np.interp(np.arange(len(f0l)), np.arange(ix), imap[:ix])
    
    # Phase function for fractional pitch
    phs = fractpitch2(fftl)
    
    fftl2 = fftl // 2
    nsyn = len(sy)
    idx = 0
    bb = np.arange(fftl)
    rbb2 = np.arange(fftl // 2 - 1, 0, -1)  # fftl/2:-1:2 in MATLAB
    
    # Parameters for noise-based apf design
    t = (np.arange(fftl) - fftl / 2) / fftl * 2
    adjd = 1.0 / (1 + np.exp(-20 * t))
    gw = np.exp(-0.25 * np.pi * (fs * (t / 2) / gdbw) ** 2)
    gw = gw / np.sum(gw)
    fgw = np.real(np.fft.fft(np.fft.fftshift(gw)))
    df = fs / fftl * 2 * np.pi
    fw = np.arange(1, fftl2 + 2) / fftl * fs
    
    trbw = 300
    rho = 1.0 / (1 + np.exp(-(fw - cornf) / trbw))
    
    # Frozen group delay component calculation
    nz = np.random.randn(fftl2 + 1) * rho
    lft = 1 - hanning(fftl) + nz[0] * 0
    lft = 1.0 / (1 + np.exp(-(lft - 0.5) * 60))
    ww = 1.0 / (1 + np.exp(-(hanning(fftl) - 0.3) * 23))
    
    iin = 0
    icntr = 0
    dmx = np.max(n2sgram)
    
    # Voiced part synthesis
    while (idx < nsyn - fftl - 10) and (int(np.ceil(iin)) < len(f0l) - 1):
        icntr += 1
        iix = int(np.round(imap[int(np.round(idx))]))
        ii = np.clip(iix, 0, min(njj, len(f0l)) - 1)
        
        f0 = f0l[ii]
        if f0 == 0:
            f0 = 200
        else:
            f0 = max(lowestF0 / pcnv, f0l[ii])
        f0 = f0 * pcnv
        
        # Look ahead correction of F0
        tnf0 = fs / f0
        tidx = idx + tnf0
        if tidx < len(imap):
            tiix = int(np.round(imap[int(np.round(tidx))]))
            tii = np.clip(tiix, 0, min(njj, len(f0l)) - 1)
            tf0 = f0l[tii]
            
            if (tf0 > 0) and (f0l[ii] > 0):
                mid_ii = int(np.round((ii + tii) / 2))
                mid_ii = np.clip(mid_ii, 0, len(f0l) - 1)
                if f0l[mid_ii] > 0:
                    f0 = max(lowestF0 / pcnv, f0l[mid_ii])
                else:
                    f0 = f0l[ii]
                f0 = f0 * pcnv
        
        # Spectrum extraction and minimum phase conversion
        idcv_int = np.round(idcv).astype(int)
        idcv_int = np.clip(idcv_int, 0, n2sgram.shape[0] - 1)
        ff = np.concatenate([
            n2sgram[idcv_int, ii],
            n2sgram[idcv_int[rbb2], ii]
        ])
        
        ccp = np.real(np.fft.fft(np.log(ff + dmx / 1000000)))
        ccp2 = np.concatenate([
            [ccp[0]],
            2 * ccp[1:fftl // 2],
            np.zeros(fftl - fftl // 2)
        ])
        ffx = np.fft.fft(ccp2 * lft) / fftl
        
        nidx = int(np.round(idx))
        nf0 = fs / f0
        frt = idx - nidx
        frtz = np.exp(1j * phs * frt)
        
        # Random group delay
        nz = np.random.randn(fftl2 + 1) * rho
        nz = np.real(np.fft.ifft(np.fft.fft(np.concatenate([nz, nz[rbb2]])) * fgw))
        nz = nz * np.sqrt(fftl * gdbw / fs)
        
        if delfracind:
            delsp = delfrac * 1000 / f0
        
        nz = nz * delsp * df / 1000
        mz = np.cumsum(np.concatenate([nz[:fftl2 + 1], nz[rbb2]])) - nz[0]
        # MATLAB uses 1-based indexing: mz(fftl) and mz(2)
        # In Python 0-based: mz[fftl-1] and mz[1]
        mmz = -(mz - adjd * (np.remainder(mz[fftl-1] + mz[1], 2 * np.pi) - 2 * np.pi))
        pzr = np.exp(-1j * mmz)
        
        pz = pzr
        wnz = aprm[idcv_int, ii]
        wpr = np.sqrt(np.maximum(0, 1 - wnz * wnz))
        
        # Temporal envelope control of aperiodic component
        nf0n = int(np.round(nf0))
        rx = np.random.randn(nf0n)
        zt0 = nf0 / fs + rx[0] * 0
        ztc = 0.01
        ztp = np.arange(nf0n) / fs
        nev = np.sqrt(2 * zt0 / ztc / (1 - np.exp(-2 * zt0 / ztc))) * np.exp(-ztp / ztc)
        rx = np.random.randn(nf0n)
        wfv = np.fft.fft((rx - np.mean(rx)) * nev, fftl)
        
        # Generate periodic and aperiodic components
        ep = np.zeros_like(ffx, dtype=float)
        gh = hanning(nf0n * 2)
        ep[:nf0n] = gh[nf0n-1::-1]
        ep[-nf0n+1:] = ep[1:nf0n]
        ep = -ep / np.sum(ep)
        ep[0] = ep[0] + 1
        epf = np.fft.fft(ep)
        
        wpr_full = np.concatenate([wpr, wpr[rbb2]])
        wnz_full = np.concatenate([wnz, wnz[rbb2]])
        
        tx = np.fft.fftshift(np.real(np.fft.ifft(
            epf * np.exp(ffx) * pz * frtz * wpr_full
        ))) * ww
        tx2 = np.fft.fftshift(np.real(np.fft.ifft(
            np.exp(ffx) * frtz * wnz_full * wfv
        ))) * ww
        
        if nidx + fftl <= len(sy):
            sy[bb + nidx] += (tx * np.sqrt(nf0) + tx2) * (f0raw[ii] > 0)
        
        idx += nf0
        iin = np.clip(imap[int(np.round(idx))], 0, min(njj, len(f0raw)) - 1)
        
        # V/UV transition handling
        if (f0raw[ii] == 0) and (f0raw[int(iin)] > 0):
            idxo = idx
            voiced_indices = np.where(f0raw[ii:int(iin)+1] > 0)[0]
            if len(voiced_indices) > 0:
                ipos = voiced_indices[0] + ii
                idx = max(idxo - nf0 + 1, rmap[ipos])
            else:
                idx = idxo
    
    # Unvoiced part synthesis
    ii = 0
    idx = 0
    f0 = 1000
    
    while (idx < nsyn - fftl) and (ii < len(f0l) - 1):
        nidx = int(np.round(idx))
        if f0raw[ii] == 0:
            ff = np.concatenate([
                n2sgram[idcv_int, ii],
                n2sgram[idcv_int[rbb2], ii]
            ])
            ccp = np.real(np.fft.fft(np.log(ff + dmx / 100000)))
            ccp2 = np.concatenate([
                [ccp[0]],
                2 * ccp[1:fftl // 2],
                np.zeros(fftl - fftl // 2)
            ])
            ffx = np.fft.fft(ccp2 * lft) / fftl
            nf0 = fs / f0
            tx = np.fft.fftshift(np.real(np.fft.ifft(np.exp(ffx))))
            rx = np.random.randn(int(np.round(nf0)))
            tnx = fftfilt(rx - np.mean(rx), tx)
            
            if nidx + fftl <= len(sy):
                sy[bb + nidx] += tnx[bb] * ww
        
        idx += nf0
        ii = int(np.round(imap[int(np.round(idx))]))
        ii = np.clip(ii, 0, len(f0l) - 1)
    
    # Extract synthesized signal
    sy2 = sy[fftl // 2:fftl // 2 + ix]
    sy = sy2
    
    return sy, statusReport


def fractpitch2(fftl: int) -> np.ndarray:
    """
    Phase rotator for fractional pitch
    
    Port of fractpitch2.m
    
    Args:
        fftl: FFT length
        
    Returns:
        phs: Phase rotator array
    """
    amp = 15
    t = (np.arange(fftl) - fftl / 2) / fftl * 2
    phs = t + (1 - np.exp(amp * t)) / (1 + np.exp(amp * t)) - \
          (1 + (1 - np.exp(amp)) / (1 + np.exp(amp))) * t
    phs[0] = 0
    phs = phs * np.pi
    return phs


def powerchk(x: np.ndarray, fs: float, segms: float) -> float:
    """
    Calculate average power of voiced portion
    
    Port of zpowerchk function from exstraightsynth.m
    
    Args:
        x: Signal
        fs: Sampling frequency (Hz)
        segms: Segment length (ms)
        
    Returns:
        pow: Average power in dB
    """
    x1 = x.flatten().copy()
    
    # Replace NaN with small values
    iv = np.arange(len(x1))
    nan_mask = np.isnan(x1)
    x1[nan_mask] = 1e-10
    
    x2 = x1 * x1
    n = int(np.round(segms / 1000 * fs))
    nw = int(np.ceil(len(x) / n))
    
    # Pad if necessary
    if len(x) % n > 0:
        padding = n * nw - len(x)
        x2 = np.concatenate([x2, (1e-6 * np.random.randn(padding)) ** 2])
    
    x2[x2 == 0] = 1e-6
    
    # Reshape and compute segment power
    pw = np.sum(x2.reshape(nw, n), axis=1) / n
    
    # Average power of high-power segments
    pow = 10 * np.log10(np.mean(pw[pw > (np.mean(pw) / 30)]))
    
    return pow


# For testing
if __name__ == '__main__':
    print("STRAIGHT Synthesis Module")
    print("Use: from straight.synthesis import exstraightsynth")

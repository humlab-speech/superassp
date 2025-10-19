"""
Voicebox Compatibility Layer

This module provides Python implementations of commonly used functions from
Mike Brookes' Voicebox MATLAB toolkit that are dependencies in COVAREP.

The goal is to provide API-compatible replacements using NumPy/SciPy where possible,
while maintaining numerical accuracy.

Original Voicebox: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
"""

import numpy as np
from scipy import signal as sp_signal
from scipy.interpolate import interp1d
from scipy.fft import rfft, irfft, fft, ifft
import warnings

__all__ = [
    'rfft_voicebox',
    'irfft_voicebox',
    'frq2mel',
    'mel2frq',
    'frq2erb',
    'erb2frq',
    'frq2bark',
    'bark2frq',
    'enframe',
    'v_windows',
    'zerocros',
    'activlev',
]


def rfft_voicebox(x, n=None, axis=-1):
    """
    Real FFT compatible with MATLAB's rfft behavior from voicebox
    
    Parameters
    ----------
    x : array_like
        Input signal
    n : int, optional
        FFT length
    axis : int, optional
        Axis along which to compute FFT
        
    Returns
    -------
    X : ndarray
        Complex FFT coefficients
    """
    if n is None:
        n = x.shape[axis]
    return rfft(x, n=n, axis=axis)


def irfft_voicebox(X, n=None, axis=-1):
    """
    Inverse real FFT compatible with voicebox
    
    Parameters
    ----------
    X : array_like
        Complex FFT coefficients
    n : int, optional
        Output length
    axis : int, optional
        Axis along which to compute IFFT
        
    Returns
    -------
    x : ndarray
        Real signal
    """
    if n is None:
        n = (X.shape[axis] - 1) * 2
    return irfft(X, n=n, axis=axis)


def frq2mel(frq):
    """
    Convert frequency in Hz to mel scale
    
    Uses the formula: mel = 2595 * log10(1 + f/700)
    
    Parameters
    ----------
    frq : array_like
        Frequency in Hz
        
    Returns
    -------
    mel : ndarray
        Mel scale frequency
    """
    frq = np.asarray(frq)
    return 2595.0 * np.log10(1.0 + frq / 700.0)


def mel2frq(mel):
    """
    Convert mel scale to frequency in Hz
    
    Inverse of frq2mel
    
    Parameters
    ----------
    mel : array_like
        Mel scale frequency
        
    Returns
    -------
    frq : ndarray
        Frequency in Hz
    """
    mel = np.asarray(mel)
    return 700.0 * (np.power(10.0, mel / 2595.0) - 1.0)


def frq2erb(frq):
    """
    Convert frequency in Hz to ERB (Equivalent Rectangular Bandwidth) scale
    
    Uses Glasberg & Moore formula
    
    Parameters
    ----------
    frq : array_like
        Frequency in Hz
        
    Returns
    -------
    erb : ndarray
        ERB scale frequency
    """
    frq = np.asarray(frq)
    return 11.17268 * np.log(1.0 + 46.06538 * frq / (frq + 14678.49))


def erb2frq(erb):
    """
    Convert ERB scale to frequency in Hz
    
    Inverse of frq2erb
    
    Parameters
    ----------
    erb : array_like
        ERB scale frequency
        
    Returns
    -------
    frq : ndarray
        Frequency in Hz
    """
    erb = np.asarray(erb)
    return 676170.4 * (1.0 / (47.06538 - np.exp(erb * 0.08950404)) - 1.0)


def frq2bark(frq):
    """
    Convert frequency in Hz to Bark scale
    
    Uses Traunmüller formula
    
    Parameters
    ----------
    frq : array_like
        Frequency in Hz
        
    Returns
    -------
    bark : ndarray
        Bark scale frequency
    """
    frq = np.asarray(frq)
    return 26.81 * frq / (1960.0 + frq) - 0.53


def bark2frq(bark):
    """
    Convert Bark scale to frequency in Hz
    
    Inverse of frq2bark
    
    Parameters
    ----------
    bark : array_like
        Bark scale frequency
        
    Returns
    -------
    frq : ndarray
        Frequency in Hz
    """
    bark = np.asarray(bark)
    return 1960.0 * (bark + 0.53) / (26.81 - bark - 0.53)


def enframe(x, win, hop=None, m=None):
    """
    Divide signal into overlapping frames
    
    Parameters
    ----------
    x : array_like
        Input signal
    win : int or array_like
        Window length (int) or window function (array)
    hop : int, optional
        Hop size between frames (default: win//2)
    m : int, optional
        Pad signal to this length (default: no padding)
        
    Returns
    -------
    frames : ndarray
        2D array of frames (n_frames, frame_length)
    """
    x = np.asarray(x).ravel()
    
    # Handle window
    if np.isscalar(win):
        win_len = int(win)
        window = np.ones(win_len)
    else:
        window = np.asarray(win)
        win_len = len(window)
    
    # Handle hop
    if hop is None:
        hop = win_len // 2
    hop = int(hop)
    
    # Pad if needed
    if m is not None:
        if m > len(x):
            x = np.pad(x, (0, m - len(x)), mode='constant')
    
    # Calculate number of frames
    n_frames = 1 + (len(x) - win_len) // hop
    
    if n_frames < 1:
        return np.array([]).reshape(0, win_len)
    
    # Create frame indices
    indices = np.arange(win_len)[None, :] + np.arange(n_frames)[:, None] * hop
    
    # Extract frames and apply window
    frames = x[indices] * window[None, :]
    
    return frames


def v_windows(n, mode='hann', p1=None, p2=None):
    """
    Generate various window functions compatible with voicebox
    
    Parameters
    ----------
    n : int
        Window length
    mode : str
        Window type: 'hann', 'hamming', 'blackman', 'kaiser', etc.
    p1, p2 : float, optional
        Parameters for certain window types
        
    Returns
    -------
    w : ndarray
        Window function
    """
    mode = mode.lower()
    
    if mode in ['hann', 'hanning']:
        return np.hanning(n)
    elif mode == 'hamming':
        return np.hamming(n)
    elif mode == 'blackman':
        return np.blackman(n)
    elif mode == 'bartlett':
        return np.bartlett(n)
    elif mode == 'kaiser':
        beta = p1 if p1 is not None else 8.6
        return np.kaiser(n, beta)
    elif mode == 'rectangular' or mode == 'rect':
        return np.ones(n)
    else:
        warnings.warn(f"Unknown window type '{mode}', using Hann window")
        return np.hanning(n)


def zerocros(x, mode='b'):
    """
    Find zero crossings in a signal
    
    Parameters
    ----------
    x : array_like
        Input signal
    mode : str
        'b' - both directions
        'p' - positive-going only
        'n' - negative-going only
        
    Returns
    -------
    zc : ndarray
        Indices of zero crossings
    """
    x = np.asarray(x).ravel()
    
    # Find sign changes
    signs = np.sign(x)
    signs[signs == 0] = 1  # Treat zeros as positive
    
    # Detect crossings
    crossings = np.diff(signs)
    
    if mode == 'p':
        # Positive-going only
        indices = np.where(crossings > 0)[0]
    elif mode == 'n':
        # Negative-going only
        indices = np.where(crossings < 0)[0]
    else:
        # Both directions
        indices = np.where(crossings != 0)[0]
    
    return indices


def activlev(x, fs, mode='a'):
    """
    Measure active speech level according to ITU-T P.56
    
    Parameters
    ----------
    x : array_like
        Input signal
    fs : float
        Sampling frequency
    mode : str
        Activity detection mode
        'a' - standard
        
    Returns
    -------
    lev : float
        Active level in dB
    af : float
        Activity factor (0-1)
    """
    x = np.asarray(x).ravel()
    
    # Simple energy-based activity detection
    # (Simplified version - full P.56 is more complex)
    frame_len = int(0.016 * fs)  # 16 ms frames
    hop = frame_len // 2
    
    frames = enframe(x, frame_len, hop)
    
    if len(frames) == 0:
        return -np.inf, 0.0
    
    # Compute frame energies
    energies = np.sum(frames**2, axis=1)
    
    # Threshold for activity (simplified)
    threshold = np.percentile(energies, 30)
    active_frames = energies > threshold
    
    # Activity factor
    af = np.sum(active_frames) / len(active_frames)
    
    # Active level
    if af > 0:
        active_energy = np.mean(energies[active_frames])
        lev = 10 * np.log10(active_energy + 1e-10)
    else:
        lev = -np.inf
    
    return lev, af


# Additional utility functions will be added as needed
# This is the foundation for the voicebox compatibility layer

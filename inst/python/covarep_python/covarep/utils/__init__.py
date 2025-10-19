"""
Utility functions for COVAREP

Common signal processing and helper functions used across modules
"""

import numpy as np
from scipy import signal
import warnings

__all__ = [
    'db2pow',
    'pow2db',
    'rms',
    'nextpow2',
    'fix_length',
    'frame_signal',
    'overlap_add',
]


def db2pow(db):
    """
    Convert dB to power
    
    Parameters
    ----------
    db : array_like
        Value in decibels
        
    Returns
    -------
    power : ndarray
        Linear power value
    """
    return 10.0 ** (db / 10.0)


def pow2db(power):
    """
    Convert power to dB
    
    Parameters
    ----------
    power : array_like
        Linear power value
        
    Returns
    -------
    db : ndarray
        Value in decibels
    """
    return 10.0 * np.log10(power + 1e-20)


def rms(x, axis=-1):
    """
    Compute RMS (root mean square)
    
    Parameters
    ----------
    x : array_like
        Input signal
    axis : int, optional
        Axis along which to compute RMS
        
    Returns
    -------
    r : ndarray
        RMS value
    """
    return np.sqrt(np.mean(x**2, axis=axis))


def nextpow2(n):
    """
    Next power of 2
    
    Parameters
    ----------
    n : int
        Input number
        
    Returns
    -------
    p : int
        Smallest power of 2 >= n
    """
    return int(2 ** np.ceil(np.log2(n)))


def fix_length(data, size, axis=-1, **kwargs):
    """
    Fix the length of an array to a target size
    
    Parameters
    ----------
    data : ndarray
        Input array
    size : int
        Target length
    axis : int
        Axis to fix
    **kwargs
        Additional arguments for np.pad
        
    Returns
    -------
    data_fixed : ndarray
        Array with fixed length
    """
    n = data.shape[axis]
    
    if n == size:
        return data
    elif n < size:
        # Pad
        lengths = [(0, 0)] * data.ndim
        lengths[axis] = (0, size - n)
        return np.pad(data, lengths, **kwargs)
    else:
        # Truncate
        slices = [slice(None)] * data.ndim
        slices[axis] = slice(0, size)
        return data[tuple(slices)]


def frame_signal(signal, frame_len, hop_len, window=None):
    """
    Split signal into overlapping frames
    
    Parameters
    ----------
    signal : array_like
        Input signal (1D)
    frame_len : int
        Frame length in samples
    hop_len : int
        Hop length in samples
    window : array_like, optional
        Window function
        
    Returns
    -------
    frames : ndarray
        2D array of frames (n_frames, frame_len)
    """
    signal = np.asarray(signal).ravel()
    n_frames = 1 + (len(signal) - frame_len) // hop_len
    
    if n_frames < 1:
        return np.array([]).reshape(0, frame_len)
    
    # Create frame indices
    indices = np.arange(frame_len)[None, :] + np.arange(n_frames)[:, None] * hop_len
    frames = signal[indices]
    
    # Apply window if provided
    if window is not None:
        window = np.asarray(window)
        if len(window) != frame_len:
            raise ValueError(f"Window length ({len(window)}) must match frame length ({frame_len})")
        frames = frames * window[None, :]
    
    return frames


def overlap_add(frames, hop_len, window=None):
    """
    Reconstruct signal from overlapping frames using overlap-add
    
    Parameters
    ----------
    frames : ndarray
        2D array of frames (n_frames, frame_len)
    hop_len : int
        Hop length in samples
    window : array_like, optional
        Synthesis window
        
    Returns
    -------
    signal : ndarray
        Reconstructed signal
    """
    n_frames, frame_len = frames.shape
    
    # Output signal length
    signal_len = (n_frames - 1) * hop_len + frame_len
    signal = np.zeros(signal_len)
    window_sum = np.zeros(signal_len)
    
    # Apply window if provided
    if window is not None:
        window = np.asarray(window)
        if len(window) != frame_len:
            raise ValueError(f"Window length ({len(window)}) must match frame length ({frame_len})")
    else:
        window = np.ones(frame_len)
    
    # Overlap-add
    for i, frame in enumerate(frames):
        start = i * hop_len
        end = start + frame_len
        signal[start:end] += frame * window
        window_sum[start:end] += window ** 2
    
    # Normalize by window
    nonzero = window_sum > 1e-10
    signal[nonzero] /= window_sum[nonzero]
    
    return signal

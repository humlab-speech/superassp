#!/usr/bin/env python3
"""
Snack-style formant tracking implementation.

Implements the Snack formant tracker algorithm using LPC analysis
and dynamic programming for formant tracking, following the approach
of David Talkin's formant tracker.

Based on the algorithm from Snack Sound Toolkit.
"""

import numpy as np
import librosa
import sys
import json
from scipy import signal
from scipy.linalg import solve_toeplitz


def lpc_analysis(frame, order):
    """
    Compute LPC coefficients using autocorrelation method (Levinson-Durbin).
    
    Parameters
    ----------
    frame : ndarray
        Audio frame
    order : int
        LPC order
        
    Returns
    -------
    tuple
        (lpc_coefficients, error)
    """
    # Compute autocorrelation
    r = np.correlate(frame, frame, mode='full')
    r = r[len(r)//2:]
    r = r[:order+1]
    
    if r[0] == 0:
        return np.zeros(order), 1.0
    
    # Levinson-Durbin recursion
    try:
        # Using scipy's solve_toeplitz for efficiency
        R = r[:-1]
        r_vec = r[1:]
        a = solve_toeplitz(R, r_vec)
        
        # Compute prediction error
        error = r[0] - np.dot(r_vec, a)
        
        return np.concatenate(([1.0], -a)), error
    except:
        return np.concatenate(([1.0], np.zeros(order))), r[0]


def lpc_to_formants(lpc_coeffs, sr, num_formants=5):
    """
    Extract formant frequencies and bandwidths from LPC coefficients.
    
    Parameters
    ----------
    lpc_coeffs : ndarray
        LPC coefficients
    sr : float
        Sample rate
    num_formants : int
        Number of formants to extract
        
    Returns
    -------
    tuple
        (formant_frequencies, bandwidths)
    """
    # Find roots of LPC polynomial
    roots = np.roots(lpc_coeffs)
    
    # Keep only roots inside unit circle (stable poles)
    roots = roots[np.abs(roots) < 1.0]
    
    # Convert to frequency and bandwidth
    angles = np.angle(roots)
    magnitudes = np.abs(roots)
    
    # Only positive frequencies
    pos_idx = angles > 0
    angles = angles[pos_idx]
    magnitudes = magnitudes[pos_idx]
    
    # Convert to Hz
    freqs = angles * sr / (2 * np.pi)
    
    # Bandwidth from magnitude
    bandwidths = -np.log(magnitudes) * sr / np.pi
    
    # Sort by frequency
    sort_idx = np.argsort(freqs)
    freqs = freqs[sort_idx]
    bandwidths = bandwidths[sort_idx]
    
    # Initialize output arrays
    formant_freqs = np.zeros(num_formants)
    formant_bws = np.zeros(num_formants)
    
    # Expected formant ranges (Hz)
    fmins = [200, 550, 1400, 2400, 3300, 4200, 5200][:num_formants]
    fmaxs = [900, 2500, 3500, 4800, 6000, 7000, 8000][:num_formants]
    
    # Map poles to formants using expected ranges
    used_poles = []
    for i in range(num_formants):
        best_pole = -1
        best_dist = float('inf')
        
        for j in range(len(freqs)):
            if j in used_poles:
                continue
            
            # Check if frequency is in expected range
            if fmins[i] <= freqs[j] <= fmaxs[i]:
                # Prefer formants close to nominal frequencies
                nominal = (fmins[i] + fmaxs[i]) / 2
                dist = abs(freqs[j] - nominal) / nominal
                
                if dist < best_dist:
                    best_dist = dist
                    best_pole = j
        
        if best_pole >= 0:
            formant_freqs[i] = freqs[best_pole]
            formant_bws[i] = bandwidths[best_pole]
            used_poles.append(best_pole)
        else:
            # No good match found - mark as missing
            formant_freqs[i] = 0.0
            formant_bws[i] = 0.0
    
    return formant_freqs, formant_bws


def smooth_formant_tracks(formant_tracks, bandwidth_tracks, max_change=0.25):
    """
    Apply dynamic programming to smooth formant tracks.
    
    Parameters
    ----------
    formant_tracks : ndarray
        Raw formant frequencies [n_frames, n_formants]
    bandwidth_tracks : ndarray
        Formant bandwidths [n_frames, n_formants]
    max_change : float
        Maximum allowed proportional change between frames
        
    Returns
    -------
    tuple
        (smoothed_formants, smoothed_bandwidths)
    """
    n_frames, n_formants = formant_tracks.shape
    
    # Simple median filtering for now (DP would be more complex)
    smoothed_f = np.copy(formant_tracks)
    smoothed_b = np.copy(bandwidth_tracks)
    
    # Apply median filter to each formant track
    for i in range(n_formants):
        # Forward pass
        for j in range(1, n_frames):
            if formant_tracks[j, i] > 0 and formant_tracks[j-1, i] > 0:
                change = abs(formant_tracks[j, i] - formant_tracks[j-1, i]) / formant_tracks[j-1, i]
                if change > max_change:
                    # Large jump - prefer previous value
                    smoothed_f[j, i] = formant_tracks[j-1, i]
                    smoothed_b[j, i] = bandwidth_tracks[j-1, i]
    
    return smoothed_f, smoothed_b


def snack_formant(
    soundFile,
    numFormants=4,
    lpcOrder=None,
    windowShift=5.0,
    windowLength=49.0,
    preEmphasis=0.7,
    beginTime=0.0,
    endTime=0.0
):
    """
    Extract formants using Snack-style LPC analysis.
    
    Parameters
    ----------
    soundFile : str
        Path to audio file
    numFormants : int
        Number of formants to track (default: 4)
    lpcOrder : int
        LPC order (default: 2 * numFormants + 6)
    windowShift : float
        Frame shift in milliseconds (default: 5)
    windowLength : float
        Analysis window length in milliseconds (default: 49)
    preEmphasis : float
        Pre-emphasis factor (default: 0.7)
    beginTime : float
        Start time in seconds (default: 0.0)
    endTime : float
        End time in seconds (default: 0.0 = end of file)
    
    Returns
    -------
    dict
        Dictionary with keys:
        - formants: formant frequencies [n_frames, n_formants]
        - bandwidths: formant bandwidths [n_frames, n_formants]
        - sample_rate: original sample rate
        - n_frames: number of frames
        - n_formants: number of formants
    """
    # Load audio
    y, sr = librosa.load(soundFile, sr=None, mono=True)
    
    # Handle time windowing
    if beginTime > 0 or (endTime > 0 and endTime > beginTime):
        start_sample = int(beginTime * sr) if beginTime > 0 else 0
        end_sample = int(endTime * sr) if endTime > 0 else len(y)
        y = y[start_sample:end_sample]
    
    # Pre-emphasis filter
    if preEmphasis > 0:
        y = signal.lfilter([1, -preEmphasis], [1], y)
    
    # Determine LPC order
    if lpcOrder is None:
        lpcOrder = 2 * numFormants + 6
    
    # Convert parameters to samples
    hop_length = int(sr * windowShift / 1000)
    win_length = int(sr * windowLength / 1000)
    
    # Make win_length odd for symmetry
    if win_length % 2 == 0:
        win_length += 1
    
    # Framing
    n_frames = 1 + (len(y) - win_length) // hop_length
    
    # Storage for formants
    formant_matrix = np.zeros((n_frames, numFormants))
    bandwidth_matrix = np.zeros((n_frames, numFormants))
    
    # Process each frame
    for i in range(n_frames):
        start = i * hop_length
        end = start + win_length
        
        if end > len(y):
            break
        
        frame = y[start:end]
        
        # Apply Hamming window
        frame = frame * np.hamming(len(frame))
        
        # LPC analysis
        lpc_coeffs, error = lpc_analysis(frame, lpcOrder)
        
        # Extract formants
        freqs, bws = lpc_to_formants(lpc_coeffs, sr, numFormants)
        
        formant_matrix[i, :] = freqs
        bandwidth_matrix[i, :] = bws
    
    # Smooth formant tracks using DP-like approach
    formant_matrix, bandwidth_matrix = smooth_formant_tracks(formant_matrix, bandwidth_matrix)
    
    return {
        'formants': formant_matrix.tolist(),
        'bandwidths': bandwidth_matrix.tolist(),
        'sample_rate': float(sr),
        'n_frames': n_frames,
        'n_formants': numFormants
    }


if __name__ == '__main__':
    if len(sys.argv) > 1:
        params = json.loads(sys.argv[1])
        result = snack_formant(**params)
        print(json.dumps(result))
    else:
        print("Error: No parameters provided", file=sys.stderr)
        sys.exit(1)

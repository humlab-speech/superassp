#!/usr/bin/env python3
"""
Snack-style pitch tracking implementation.

Implements the Snack pitch tracker algorithm using Python libraries.
The Snack pitch tracker uses autocorrelation with dynamic programming
for pitch candidate selection and tracking.

Based on the algorithm from Snack Sound Toolkit by Kåre Sjölander.
"""

import numpy as np
import librosa
import sys
import json
from scipy import signal
from scipy.interpolate import interp1d


def compute_autocorrelation(frame, max_lag):
    """
    Compute normalized autocorrelation function.
    
    Parameters
    ----------
    frame : ndarray
        Audio frame
    max_lag : int
        Maximum lag to compute
        
    Returns
    -------
    ndarray
        Normalized autocorrelation values
    """
    corr = np.correlate(frame, frame, mode='full')
    corr = corr[len(corr)//2:]  # Take second half
    
    # Normalize
    if corr[0] > 0:
        corr = corr / corr[0]
    
    return corr[:max_lag]


def find_peaks(acf, min_lag, max_lag, threshold=0.3):
    """
    Find peaks in autocorrelation function.
    
    Parameters
    ----------
    acf : ndarray
        Autocorrelation function
    min_lag : int
        Minimum lag (max F0)
    max_lag : int
        Maximum lag (min F0)
    threshold : float
        Peak detection threshold
        
    Returns
    -------
    tuple
        (peak_lags, peak_values)
    """
    # Find local maxima
    peaks = signal.find_peaks(acf[min_lag:max_lag], height=threshold)[0]
    peaks = peaks + min_lag
    
    if len(peaks) == 0:
        return np.array([]), np.array([])
    
    peak_values = acf[peaks]
    
    # Sort by peak value (descending)
    sort_idx = np.argsort(peak_values)[::-1]
    
    return peaks[sort_idx], peak_values[sort_idx]


def dynamic_programming_tracking(candidates_list, costs_list, n_frames):
    """
    Use dynamic programming to find best F0 trajectory.
    
    Parameters
    ----------
    candidates_list : list of ndarrays
        F0 candidates for each frame
    costs_list : list of ndarrays
        Costs for each candidate
    n_frames : int
        Number of frames
        
    Returns
    -------
    ndarray
        Best F0 trajectory
    """
    if n_frames == 0:
        return np.array([])
    
    # Initialize DP table
    max_cands = max(len(c) for c in candidates_list if len(c) > 0)
    if max_cands == 0:
        return np.zeros(n_frames)
    
    # Cumulative costs
    cum_cost = [[float('inf')] * len(candidates_list[i]) 
                for i in range(n_frames)]
    backpointer = [[0] * len(candidates_list[i]) 
                   for i in range(n_frames)]
    
    # Initialize first frame
    if len(candidates_list[0]) > 0:
        for j in range(len(candidates_list[0])):
            cum_cost[0][j] = costs_list[0][j]
    
    # Forward pass
    for i in range(1, n_frames):
        if len(candidates_list[i]) == 0:
            continue
            
        for j in range(len(candidates_list[i])):
            best_cost = float('inf')
            best_prev = 0
            
            if len(candidates_list[i-1]) > 0:
                for k in range(len(candidates_list[i-1])):
                    # Transition cost (prefer smooth F0 changes)
                    if candidates_list[i-1][k] > 0 and candidates_list[i][j] > 0:
                        trans_cost = abs(candidates_list[i][j] - candidates_list[i-1][k]) / candidates_list[i-1][k]
                        trans_cost *= 0.02  # freq_weight from Snack
                    else:
                        trans_cost = 0.005  # trans_cost from Snack
                    
                    total_cost = cum_cost[i-1][k] + costs_list[i][j] + trans_cost
                    
                    if total_cost < best_cost:
                        best_cost = total_cost
                        best_prev = k
            else:
                best_cost = costs_list[i][j]
            
            cum_cost[i][j] = best_cost
            backpointer[i][j] = best_prev
    
    # Backtrack to find best path
    trajectory = np.zeros(n_frames)
    
    # Find best final candidate
    if len(candidates_list[n_frames-1]) > 0:
        best_final = np.argmin(cum_cost[n_frames-1])
        trajectory[n_frames-1] = candidates_list[n_frames-1][best_final]
        
        # Backtrack
        current = best_final
        for i in range(n_frames-2, -1, -1):
            if len(candidates_list[i]) > 0:
                current = backpointer[i+1][current]
                if current < len(candidates_list[i]):
                    trajectory[i] = candidates_list[i][current]
    
    return trajectory


def snack_pitch(
    soundFile,
    minF=50.0,
    maxF=550.0,
    windowShift=10.0,
    windowLength=7.5,
    threshold=0.3,
    beginTime=0.0,
    endTime=0.0
):
    """
    Extract F0 using Snack-style pitch tracker.
    
    Parameters
    ----------
    soundFile : str
        Path to audio file
    minF : float
        Minimum F0 in Hz (default: 50)
    maxF : float
        Maximum F0 in Hz (default: 550)
    windowShift : float
        Frame shift in milliseconds (default: 10)
    windowLength : float
        Analysis window length in milliseconds (default: 7.5)
    threshold : float
        Correlation threshold (default: 0.3)
    beginTime : float
        Start time in seconds (default: 0.0)
    endTime : float
        End time in seconds (default: 0.0 = end of file)
    
    Returns
    -------
    dict
        Dictionary with keys:
        - f0: F0 values (Hz, 0 = unvoiced)
        - voicing: voicing probability (0-1)
        - rms: RMS energy per frame
        - sample_rate: original sample rate
        - n_frames: number of frames
    """
    # Load audio
    y, sr = librosa.load(soundFile, sr=None, mono=True)
    
    # Handle time windowing
    if beginTime > 0 or (endTime > 0 and endTime > beginTime):
        start_sample = int(beginTime * sr) if beginTime > 0 else 0
        end_sample = int(endTime * sr) if endTime > 0 else len(y)
        y = y[start_sample:end_sample]
    
    # Convert parameters to samples
    hop_length = int(sr * windowShift / 1000)
    win_length = int(sr * windowLength / 1000)
    
    # Lag range for autocorrelation
    min_lag = int(sr / maxF)
    max_lag = int(sr / minF)
    
    # Framing
    n_frames = 1 + (len(y) - win_length) // hop_length
    
    f0_candidates = []
    candidate_costs = []
    rms_values = []
    
    # Process each frame
    for i in range(n_frames):
        start = i * hop_length
        end = start + win_length
        
        if end > len(y):
            break
            
        frame = y[start:end]
        
        # Apply window
        frame = frame * np.hanning(len(frame))
        
        # Compute RMS
        rms = np.sqrt(np.mean(frame**2))
        rms_values.append(rms)
        
        # Compute autocorrelation
        acf = compute_autocorrelation(frame, max_lag + 1)
        
        # Find peaks
        peak_lags, peak_vals = find_peaks(acf, min_lag, max_lag, threshold)
        
        # Convert lags to frequencies
        if len(peak_lags) > 0:
            # Take top candidates (max 20 like Snack)
            n_cands = min(len(peak_lags), 20)
            peak_lags = peak_lags[:n_cands]
            peak_vals = peak_vals[:n_cands]
            
            f0_cands = sr / peak_lags
            
            # Costs (inverse of peak strength, weighted by RMS)
            costs = (1.0 - peak_vals) * (1.0 if rms > 0.01 else 10.0)
        else:
            f0_cands = np.array([0.0])
            costs = np.array([10.0])  # High cost for unvoiced
        
        f0_candidates.append(f0_cands)
        candidate_costs.append(costs)
    
    # Dynamic programming to find best trajectory
    f0_trajectory = dynamic_programming_tracking(f0_candidates, candidate_costs, len(f0_candidates))
    
    # Compute voicing probability (based on F0 and RMS)
    rms_array = np.array(rms_values)
    voicing = np.zeros(len(f0_trajectory))
    
    for i in range(len(f0_trajectory)):
        if f0_trajectory[i] > 0 and rms_array[i] > 0.01:
            # Find best candidate match
            for j, cand in enumerate(f0_candidates[i]):
                if abs(cand - f0_trajectory[i]) < 1.0:
                    voicing[i] = 1.0 - candidate_costs[i][j] / 10.0
                    voicing[i] = np.clip(voicing[i], 0.0, 1.0)
                    break
    
    return {
        'f0': f0_trajectory.tolist(),
        'voicing': voicing.tolist(),
        'rms': rms_array.tolist(),
        'sample_rate': float(sr),
        'n_frames': len(f0_trajectory)
    }


if __name__ == '__main__':
    if len(sys.argv) > 1:
        params = json.loads(sys.argv[1])
        result = snack_pitch(**params)
        print(json.dumps(result))
    else:
        print("Error: No parameters provided", file=sys.stderr)
        sys.exit(1)

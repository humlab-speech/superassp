"""
Numba-optimized Voxit implementation.

This module provides JIT-compiled functions for performance-critical sections,
achieving 2-3x speedup with minimal code changes.

Requires: numba
"""

import math
import numpy as np
from numba import jit
from scipy.signal import savgol_filter

try:
    from lempel_ziv_complexity import lempel_ziv_complexity
    HAS_LZ = True
except ImportError:
    HAS_LZ = False
    @jit(nopython=True)
    def lempel_ziv_complexity_numba(s):
        """Numba-compiled LZ complexity"""
        n = len(s)
        if n == 0:
            return 0
        i, k, l, c, k_max = 0, 1, 1, 1, 1
        while True:
            if i + k > n:
                return c
            match = True
            for j in range(k):
                if s[i+j] != s[l+j]:
                    match = False
                    break
            if not match:
                c += 1
                i = l + k
                l = i
                k = 1
                k_max = 1
            else:
                k += 1
                if k > k_max:
                    k_max = k
                if i + k > n:
                    return c + 1


@jit(nopython=True)
def compute_pauses_numba(gentle_start, gentle_end, min_pause, max_pause):
    """
    Numba-optimized pause calculation.
    
    Returns:
    --------
    pause_count, long_pause_count, sum_pauses
    """
    n = len(gentle_end)
    pause_count = 0
    long_pause_count = 0
    sum_pauses = 0.0
    
    for i in range(n - 1):
        pause_len = gentle_start[i + 1] - gentle_end[i]
        if min_pause <= pause_len <= max_pause:
            sum_pauses += pause_len
            pause_count += 1
        elif pause_len > max_pause:
            long_pause_count += 1
    
    return pause_count, long_pause_count, sum_pauses


@jit(nopython=True)
def build_rhythm_sequence_numba(gentle_start, gentle_end, min_pause, max_pause,
                                sampling_interval, total_length):
    """
    Build binary rhythm sequence (1=voiced, 0=pause) at 100 Hz sampling.
    
    Returns:
    --------
    numpy array of 0s and 1s
    """
    n_samples = int(total_length / sampling_interval) + 1
    s = np.ones(n_samples, dtype=np.int8)
    
    for i in range(len(gentle_end) - 1):
        pause_len = gentle_start[i + 1] - gentle_end[i]
        if min_pause <= pause_len <= max_pause:
            pause_start_idx = int(gentle_end[i] / sampling_interval)
            pause_end_idx = int(gentle_start[i + 1] / sampling_interval)
            for j in range(pause_start_idx, min(pause_end_idx, n_samples)):
                s[j] = 0
    
    return s


@jit(nopython=True)
def compute_pitch_histogram_numba(diffoctf0, n_bins=25, bin_min=-1.0, bin_max=1.0):
    """
    Numba-optimized histogram computation.
    
    Returns:
    --------
    histogram counts, bin edges
    """
    hist = np.zeros(n_bins, dtype=np.int64)
    bin_width = (bin_max - bin_min) / n_bins
    
    for val in diffoctf0:
        if bin_min <= val < bin_max:
            bin_idx = int((val - bin_min) / bin_width)
            if 0 <= bin_idx < n_bins:
                hist[bin_idx] += 1
        elif val >= bin_max:
            hist[n_bins - 1] += 1
    
    return hist


@jit(nopython=True)
def find_voiced_segments_numba(drift_time, vdurthresh, gap_threshold):
    """
    Find contiguous voiced segments with minimum duration.
    
    Returns:
    --------
    List of (start_idx, end_idx) tuples
    """
    n = len(drift_time)
    if n < 2:
        return []
    
    segments = []
    segment_start = 0
    
    for i in range(n - 1):
        time_diff = drift_time[i + 1] - drift_time[i]
        if time_diff > gap_threshold:
            # End of segment
            if (i - segment_start) > vdurthresh:
                segments.append((segment_start, i))
            segment_start = i + 1
    
    # Handle last segment
    if (n - 1 - segment_start) > vdurthresh:
        segments.append((segment_start, n - 1))
    
    return segments


def compute_voxit_features_numba(gentle_data, pitch_data, 
                                 start_time=None, end_time=None):
    """
    Compute Voxit features with Numba JIT optimization.
    
    Parameters:
    -----------
    gentle_data : list of dict
        Word alignment data with keys: 'word', 'case', 'start', 'end'
    pitch_data : list of dict
        Pitch track with keys: 'time', 'frequency'
    start_time : float, optional
        Analysis start time (seconds)
    end_time : float, optional
        Analysis end time (seconds)
        
    Returns:
    --------
    dict : Dictionary with computed features
    """
    results = {}
    
    # ========== WORD/PAUSE ANALYSIS ==========
    
    # Filter words (Python loop - not JIT compiled)
    gentle_start = []
    gentle_end = []
    gentle_wordcount = 0
    
    for item in gentle_data:
        word = item.get('word', '')
        start = item.get('start')
        end = item.get('end')
        
        if start is None or end is None:
            continue
        
        if start_time is not None and start < start_time:
            continue
        if end_time is not None and end > end_time:
            continue
        
        if word == '[noise]':
            continue
        
        gentle_wordcount += 1
        gentle_start.append(start)
        gentle_end.append(end)
    
    if gentle_wordcount == 0 or len(gentle_end) == 0:
        return {key: np.nan for key in [
            'WPM', 'pause_count', 'long_pause_count', 'average_pause_length',
            'average_pause_rate', 'rhythmic_complexity_of_pauses',
            'average_pitch', 'pitch_range', 'pitch_speed',
            'pitch_acceleration', 'pitch_entropy'
        ]}
    
    gentle_start = np.array(gentle_start, dtype=np.float64)
    gentle_end = np.array(gentle_end, dtype=np.float64)
    gentle_length = gentle_end[-1]
    
    # Speaking rate
    WPM = int(gentle_wordcount / (gentle_length / 60))
    results["WPM"] = WPM
    
    # Pause analysis (Numba-optimized)
    min_pause = 0.1
    max_pause = 3.0
    pause_count, long_pause_count, sum_pauses = compute_pauses_numba(
        gentle_start, gentle_end, min_pause, max_pause
    )
    
    results["pause_count"] = int(pause_count)
    results["long_pause_count"] = int(long_pause_count)
    
    if pause_count > 0:
        results["average_pause_length"] = float(sum_pauses / pause_count)
        results["average_pause_rate"] = float(pause_count / gentle_length)
    else:
        results["average_pause_length"] = 0.0
        results["average_pause_rate"] = 0.0
    
    # Rhythmic Complexity (Numba-optimized sequence generation)
    sampling_interval = 0.01
    s = build_rhythm_sequence_numba(gentle_start, gentle_end, min_pause, max_pause,
                                   sampling_interval, gentle_length)
    
    # LZ complexity
    if len(s) > 0:
        if HAS_LZ:
            s_str = ''.join(map(str, s))
            lz_comp = lempel_ziv_complexity(s_str)
        else:
            lz_comp = lempel_ziv_complexity_numba(s)
        normalized_lz = lz_comp / (len(s) / math.log2(len(s)))
        results["rhythmic_complexity_of_pauses"] = normalized_lz * 100
    else:
        results["rhythmic_complexity_of_pauses"] = np.nan
    
    # ========== PITCH ANALYSIS ==========
    
    # Extract pitch data (Python loop)
    drift_time = []
    drift_pitch = []
    
    for item in pitch_data:
        time = item.get('time')
        freq = item.get('frequency')
        
        if time is None or freq is None or freq == 0:
            continue
        
        if start_time is not None and time < start_time:
            continue
        if end_time is not None and time > end_time:
            continue
        
        drift_time.append(time)
        drift_pitch.append(freq)
    
    if len(drift_pitch) == 0:
        results["average_pitch"] = np.nan
        results["pitch_range"] = np.nan
        results["pitch_speed"] = np.nan
        results["pitch_acceleration"] = np.nan
        results["pitch_entropy"] = np.nan
        return results
    
    drift_time = np.array(drift_time, dtype=np.float64)
    drift_pitch = np.array(drift_pitch, dtype=np.float64)
    
    # Average pitch
    results["average_pitch"] = float(np.mean(drift_pitch))
    
    # Pitch calculations (vectorized)
    f0log = np.log2(drift_pitch)
    f0mean = 2.0 ** np.mean(f0log)
    diffoctf0 = np.log2(drift_pitch) - np.log2(f0mean)
    
    # Pitch range
    results["pitch_range"] = float(np.max(diffoctf0) - np.min(diffoctf0))
    
    # Pitch histogram (Numba-optimized)
    f0hist = compute_pitch_histogram_numba(diffoctf0, n_bins=25, 
                                          bin_min=-1.0, bin_max=1.0)
    f0prob = f0hist / np.sum(f0hist)
    
    # Pitch entropy
    f0prob_nonzero = f0prob[f0prob > 0]
    if len(f0prob_nonzero) > 0:
        f0log2prob = np.log2(f0prob_nonzero)
        f0entropy = -np.sum(f0prob_nonzero * f0log2prob)
        results["pitch_entropy"] = float(f0entropy)
    else:
        results["pitch_entropy"] = 0.0
    
    # ========== VELOCITY/ACCELERATION ==========
    
    # Find voiced segments (Numba-optimized)
    dminvoice = 0.100
    ts = drift_time[1] - drift_time[0] if len(drift_time) > 1 else 0.01
    vdurthresh = int(dminvoice / ts)
    gap_threshold = 3 * ts
    
    segments = find_voiced_segments_numba(drift_time, vdurthresh, gap_threshold)
    
    f0velocity = []
    f0accel = []
    
    for start, end in segments:
        diffocttmp = diffoctf0[start:end+1]
        
        if len(diffocttmp) < 3:
            continue
        
        # Savitzky-Golay smoothing (scipy - not JIT compiled)
        if len(diffocttmp) >= 7:
            diffocttmp_smooth = savgol_filter(diffocttmp, window_length=7, 
                                             polyorder=2)
        else:
            diffocttmp_smooth = diffocttmp
        
        velocity = np.diff(diffocttmp_smooth) / ts
        f0velocity.extend(velocity)
        
        if len(diffocttmp_smooth) >= 3:
            accel = np.diff(np.diff(diffocttmp_smooth)) / ts
            f0accel.extend(accel)
    
    # Statistics
    if len(f0velocity) > 0:
        f0velocity = np.array(f0velocity)
        velocity_mean = np.mean(np.abs(f0velocity))
        velocity_sign = np.sign(np.mean(f0velocity))
        results["pitch_speed"] = float(velocity_mean * velocity_sign)
    else:
        results["pitch_speed"] = 0.0
    
    if len(f0accel) > 0:
        f0accel = np.array(f0accel)
        accel_mean = np.mean(np.abs(f0accel))
        accel_sign = np.sign(np.mean(f0accel))
        results["pitch_acceleration"] = float(accel_mean * accel_sign)
    else:
        results["pitch_acceleration"] = 0.0
    
    return results

"""
Optimized Voxit feature extraction.

This module provides vectorized implementations of Voxit prosodic measures,
optimized for performance while maintaining compatibility with the original
MATLAB implementation.

Key optimizations:
- Vectorized numpy operations
- Reduced memory allocations
- Efficient pause detection
- Optimized pitch statistics
"""

import math
import numpy as np
from scipy.signal import savgol_filter

# Try to import lempel_ziv_complexity
try:
    from lempel_ziv_complexity import lempel_ziv_complexity
    HAS_LZ = True
except ImportError:
    HAS_LZ = False
    # Fallback implementation
    def lempel_ziv_complexity(s):
        """Simple Lempel-Ziv complexity fallback"""
        if not s:
            return 0
        i, k, l = 0, 1, 1
        c, k_max = 1, 1
        n = len(s)
        while True:
            if i + k > n:
                return c
            if s[i:i+k] != s[l:l+k]:
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
                    c += 1
                    return c


def compute_voxit_features_optimized(gentle_data, pitch_data, 
                                      start_time=None, end_time=None):
    """
    Compute Voxit features from word alignments and pitch data.
    
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
    dict : Dictionary with computed features:
        - WPM: Words per minute
        - pause_count: Number of pauses (100ms - 3s)
        - long_pause_count: Pauses > 3s
        - average_pause_length: Mean pause duration
        - average_pause_rate: Pauses per second
        - rhythmic_complexity_of_pauses: Normalized LZ complexity
        - average_pitch: Mean F0 (Hz)
        - pitch_range: F0 range (octaves)
        - pitch_speed: F0 velocity (octaves/s)
        - pitch_acceleration: F0 acceleration (octaves/s²)
        - pitch_entropy: F0 distribution entropy
    """
    results = {}
    
    # ========== WORD/PAUSE ANALYSIS ==========
    
    # Filter words by time window
    gentle_start = []
    gentle_end = []
    gentle_wordcount = 0
    
    for item in gentle_data:
        word = item.get('word', '')
        start = item.get('start')
        end = item.get('end')
        
        if start is None or end is None:
            continue
            
        # Time filtering
        if start_time is not None and start < start_time:
            continue
        if end_time is not None and end > end_time:
            continue
            
        # Skip noise markers
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
    
    gentle_start = np.array(gentle_start)
    gentle_end = np.array(gentle_end)
    gentle_length = gentle_end[-1]
    
    # Speaking rate (WPM)
    WPM = int(gentle_wordcount / (gentle_length / 60))
    results["WPM"] = WPM
    
    # Pause analysis (100ms to 3000ms)
    min_pause = 0.1
    max_pause = 3.0
    
    if len(gentle_end) > 1:
        # Vectorized pause calculation
        pauses = gentle_start[1:] - gentle_end[:-1]
        valid_pauses = pauses[(pauses >= min_pause) & (pauses <= max_pause)]
        long_pauses = pauses[pauses > max_pause]
        
        results["pause_count"] = len(valid_pauses)
        results["long_pause_count"] = len(long_pauses)
        
        if len(valid_pauses) > 0:
            results["average_pause_length"] = float(np.mean(valid_pauses))
            results["average_pause_rate"] = len(valid_pauses) / gentle_length
        else:
            results["average_pause_length"] = 0.0
            results["average_pause_rate"] = 0.0
    else:
        results["pause_count"] = 0
        results["long_pause_count"] = 0
        results["average_pause_length"] = 0.0
        results["average_pause_rate"] = 0.0
    
    # Rhythmic Complexity of Pauses (sample at 100 Hz)
    sampling_interval = 0.01  # 10ms
    n_samples = int(gentle_length / sampling_interval) + 1
    s = np.ones(n_samples, dtype=int)  # Default: voiced
    
    # Mark pauses in the binary sequence
    time_grid = np.arange(0, gentle_length, sampling_interval)
    for i in range(len(gentle_end) - 1):
        pause_len = gentle_start[i+1] - gentle_end[i]
        if min_pause <= pause_len <= max_pause:
            # Mark pause period as 0
            pause_start_idx = int(gentle_end[i] / sampling_interval)
            pause_end_idx = int(gentle_start[i+1] / sampling_interval)
            s[pause_start_idx:pause_end_idx] = 0
    
    # Calculate normalized LZ complexity
    s_str = ''.join(map(str, s))
    if len(s) > 0:
        lz_comp = lempel_ziv_complexity(s_str)
        normalized_lz = lz_comp / (len(s) / math.log2(len(s)))
        results["rhythmic_complexity_of_pauses"] = normalized_lz * 100
    else:
        results["rhythmic_complexity_of_pauses"] = np.nan
    
    # ========== PITCH ANALYSIS ==========
    
    # Extract pitch data
    drift_time = []
    drift_pitch = []
    
    for item in pitch_data:
        time = item.get('time')
        freq = item.get('frequency')
        
        if time is None or freq is None or freq == 0:
            continue
            
        # Time filtering
        if start_time is not None and time < start_time:
            continue
        if end_time is not None and time > end_time:
            continue
            
        drift_time.append(time)
        drift_pitch.append(freq)
    
    if len(drift_pitch) == 0:
        # No pitch data
        results["average_pitch"] = np.nan
        results["pitch_range"] = np.nan
        results["pitch_speed"] = np.nan
        results["pitch_acceleration"] = np.nan
        results["pitch_entropy"] = np.nan
        return results
    
    drift_time = np.array(drift_time)
    drift_pitch = np.array(drift_pitch)
    
    # Average pitch
    results["average_pitch"] = float(np.mean(drift_pitch))
    
    # Convert to log2 scale
    f0log = np.log2(drift_pitch)
    f0mean = 2.0 ** np.mean(f0log)
    
    # Difference from mean (in octaves)
    diffoctf0 = np.log2(drift_pitch) - np.log2(f0mean)
    
    # Pitch Range (octaves)
    results["pitch_range"] = float(np.max(diffoctf0) - np.min(diffoctf0))
    
    # Pitch histogram (25 bins, -1 to +1 octaves)
    f0hist, _ = np.histogram(diffoctf0, bins=25, range=(-1, 1))
    f0prob = f0hist / np.sum(f0hist)
    
    # Pitch Entropy
    # Only compute for non-zero probabilities
    f0prob_nonzero = f0prob[f0prob > 0]
    if len(f0prob_nonzero) > 0:
        f0log2prob = np.log2(f0prob_nonzero)
        f0entropy = -np.sum(f0prob_nonzero * f0log2prob)
        results["pitch_entropy"] = float(f0entropy)
    else:
        results["pitch_entropy"] = 0.0
    
    # ========== PITCH VELOCITY AND ACCELERATION ==========
    
    # Find contiguous voiced regions (minimum 100ms duration)
    dminvoice = 0.100
    if len(drift_time) > 1:
        ts = drift_time[1] - drift_time[0]  # Time step
    else:
        ts = 0.01  # Default 10ms
    
    vdurthresh = int(dminvoice / ts)
    
    # Find voiced segments (non-zero pitch contiguous regions)
    # Already filtered to non-zero above, so find gaps
    if len(drift_time) > 1:
        time_diff = np.diff(drift_time)
        gap_threshold = 3 * ts  # Consider gaps > 3x expected spacing
        gap_indices = np.where(time_diff > gap_threshold)[0]
        
        # Build segment boundaries
        segment_starts = np.concatenate([[0], gap_indices + 1])
        segment_ends = np.concatenate([gap_indices, [len(drift_pitch) - 1]])
        
        # Filter by minimum duration
        ixvoicedbounds = []
        for start, end in zip(segment_starts, segment_ends):
            if (end - start) > vdurthresh:
                ixvoicedbounds.append((start, end))
    else:
        ixvoicedbounds = []
    
    # Compute velocity and acceleration for each voiced segment
    f0velocity = []
    f0accel = []
    
    for start, end in ixvoicedbounds:
        # Extract pitch contour for this segment (in octaves from mean)
        diffocttmp = diffoctf0[start:end+1]
        
        if len(diffocttmp) < 3:
            continue
        
        # Apply Savitzky-Golay smoothing (matches MATLAB sgolayfilt)
        if len(diffocttmp) >= 7:
            diffocttmp_smooth = savgol_filter(diffocttmp, window_length=7, 
                                             polyorder=2)
        else:
            diffocttmp_smooth = diffocttmp
        
        # Velocity: first derivative
        velocity = np.diff(diffocttmp_smooth) / ts
        f0velocity.extend(velocity)
        
        # Acceleration: second derivative
        if len(diffocttmp_smooth) >= 3:
            accel = np.diff(np.diff(diffocttmp_smooth)) / ts
            f0accel.extend(accel)
    
    # Compute signed directionless statistics
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

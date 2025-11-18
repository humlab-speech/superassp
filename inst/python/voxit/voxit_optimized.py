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




def _convert_gentle_to_structured(gentle_data):
    """
    Convert gentle_data list of dicts to numpy structured array for vectorized operations.
    
    Parameters:
    -----------
    gentle_data : list of dict or numpy structured array
        Word alignment data
        
    Returns:
    --------
    numpy.ndarray : Structured array with fields: start, end, is_noise
    """
    if isinstance(gentle_data, np.ndarray) and gentle_data.dtype.names:
        # Already a structured array
        return gentle_data
    
    # Pre-allocate structured array
    n = len(gentle_data)
    dtype = np.dtype([('start', 'f8'), ('end', 'f8'), ('is_noise', 'bool')])
    arr = np.zeros(n, dtype=dtype)
    
    # Fill array (avoiding Python loops where possible)
    for i, item in enumerate(gentle_data):
        word = item.get('word', '')
        start = item.get('start')
        end = item.get('end')
        
        # Handle None values
        arr['start'][i] = start if start is not None else np.nan
        arr['end'][i] = end if end is not None else np.nan
        arr['is_noise'][i] = (word == '[noise]')
    
    return arr


def _convert_pitch_to_structured(pitch_data):
    """
    Convert pitch_data list of dicts to numpy structured array for vectorized operations.
    
    Parameters:
    -----------
    pitch_data : list of dict or numpy structured array
        Pitch track data
        
    Returns:
    --------
    numpy.ndarray : Structured array with fields: time, frequency
    """
    if isinstance(pitch_data, np.ndarray) and pitch_data.dtype.names:
        # Already a structured array
        return pitch_data
    
    # Pre-allocate structured array
    n = len(pitch_data)
    dtype = np.dtype([('time', 'f8'), ('frequency', 'f8')])
    arr = np.zeros(n, dtype=dtype)
    
    # Fill array
    for i, item in enumerate(pitch_data):
        time = item.get('time')
        freq = item.get('frequency')
        
        arr['time'][i] = time if time is not None else np.nan
        arr['frequency'][i] = freq if freq is not None else 0.0
    
    return arr



def compute_voxit_features_optimized(gentle_data, pitch_data, 
                                      start_time=None, end_time=None):
    """
    Compute Voxit features from word alignments and pitch data.
    
    Parameters:
    -----------


    gentle_data : list of dict or numpy structured array
        Word alignment data with keys: 'word', 'case', 'start', 'end'
        For best performance, pass as structured array with fields: start, end, is_noise
    pitch_data : list of dict or numpy structured array
        Pitch track with keys: 'time', 'frequency'
        For best performance, pass as structured array with fields: time, frequency

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
    


    # Convert to structured array for vectorized operations
    gentle_arr = _convert_gentle_to_structured(gentle_data)
    
    # Vectorized filtering
    valid_mask = ~np.isnan(gentle_arr['start']) & ~np.isnan(gentle_arr['end'])
    valid_mask &= ~gentle_arr['is_noise']
    
    if start_time is not None:
        valid_mask &= (gentle_arr['start'] >= start_time)
    if end_time is not None:
        valid_mask &= (gentle_arr['end'] <= end_time)
    
    # Extract valid words
    gentle_start = gentle_arr['start'][valid_mask]
    gentle_end = gentle_arr['end'][valid_mask]
    gentle_wordcount = len(gentle_start)
    
    if gentle_wordcount == 0:


        return {key: np.nan for key in [
            'WPM', 'pause_count', 'long_pause_count', 'average_pause_length',
            'average_pause_rate', 'rhythmic_complexity_of_pauses',
            'average_pitch', 'pitch_range', 'pitch_speed', 
            'pitch_acceleration', 'pitch_entropy'
        ]}
    

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

    s = np.ones(n_samples, dtype=np.int8)  # Default: voiced (use int8 for memory efficiency)
    
    # Vectorized marking of pauses in the binary sequence
    if len(gentle_end) > 1:
        # Calculate all pauses at once
        pauses = gentle_start[1:] - gentle_end[:-1]
        pause_mask = (pauses >= min_pause) & (pauses <= max_pause)
        
        # Convert times to indices
        pause_start_indices = (gentle_end[:-1] / sampling_interval).astype(int)
        pause_end_indices = (gentle_start[1:] / sampling_interval).astype(int)
        
        # Mark pauses (vectorized where possible)
        for i in np.where(pause_mask)[0]:
            start_idx = pause_start_indices[i]
            end_idx = pause_end_indices[i]
            if end_idx <= n_samples:
                s[start_idx:end_idx] = 0
    
    # Calculate normalized LZ complexity
    # Convert to string efficiently using numpy
    s_str = s.tobytes().decode('latin1').translate({48: '0', 49: '1'})  # Map byte values
    # Alternative: more direct conversion
    s_str = ''.join(s.astype(str))
    


    if len(s) > 0:
        lz_comp = lempel_ziv_complexity(s_str)
        normalized_lz = lz_comp / (len(s) / math.log2(len(s)))
        results["rhythmic_complexity_of_pauses"] = normalized_lz * 100
    else:
        results["rhythmic_complexity_of_pauses"] = np.nan
    
    # ========== PITCH ANALYSIS ==========
    


    # Convert pitch data to structured array
    pitch_arr = _convert_pitch_to_structured(pitch_data)
    
    # Vectorized filtering for valid pitch data
    valid_pitch_mask = ~np.isnan(pitch_arr['time']) & (pitch_arr['frequency'] > 0)
    
    if start_time is not None:
        valid_pitch_mask &= (pitch_arr['time'] >= start_time)
    if end_time is not None:
        valid_pitch_mask &= (pitch_arr['time'] <= end_time)
    
    drift_time = pitch_arr['time'][valid_pitch_mask]
    drift_pitch = pitch_arr['frequency'][valid_pitch_mask]


    
    if len(drift_pitch) == 0:
        # No pitch data
        results["average_pitch"] = np.nan
        results["pitch_range"] = np.nan
        results["pitch_speed"] = np.nan
        results["pitch_acceleration"] = np.nan
        results["pitch_entropy"] = np.nan
        return results
    

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

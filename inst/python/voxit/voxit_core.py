"""
Core Voxit implementation (reference version).

This is a faithful Python translation of the original MATLAB code,
prioritizing clarity and correctness over performance.
"""

import math
import decimal
import numpy as np
from scipy.signal import savgol_filter

try:
    from lempel_ziv_complexity import lempel_ziv_complexity
except ImportError:
    # Simple fallback
    def lempel_ziv_complexity(s):
        if not s:
            return 0
        i, k, l, c, k_max, n = 0, 1, 1, 1, 1, len(s)
        while True:
            if i + k > n:
                return c
            if s[i:i+k] != s[l:l+k]:
                c += 1
                i, l, k, k_max = l + k, l + k, 1, 1
            else:
                k += 1
                if k > k_max:
                    k_max = k
                if i + k > n:
                    return c + 1


def compute_voxit_features(gentle_data, pitch_data, start_time=None, end_time=None):
    """
    Compute Voxit features - reference implementation.
    
    This follows the original MATLAB code structure closely.
    For better performance, use compute_voxit_features_optimized().
    """
    results = {}
    
    # GENTLE: Word timing analysis
    gentle_start = []
    gentle_end = []
    gentle_wordcount = 0
    
    for item in gentle_data:
        word = item.get('word', '')
        start = item.get('start')
        end = item.get('end')
        
        if not word or start is None or end is None:
            continue
        
        if start_time and start < start_time:
            continue
        if end_time and end > end_time:
            continue
        
        if word != '[noise]':
            gentle_wordcount += 1
            gentle_start.append(round(start * 10000) / 10000)
            gentle_end.append(round(end * 10000) / 10000)
    
    if gentle_wordcount == 0:
        return {key: float('nan') for key in [
            'WPM', 'pause_count', 'long_pause_count', 'average_pause_length',
            'average_pause_rate', 'rhythmic_complexity_of_pauses',
            'average_pitch', 'pitch_range', 'pitch_speed',
            'pitch_acceleration', 'pitch_entropy'
        ]}
    
    gentle_length = gentle_end[-1]
    
    # Speaking rate (WPM)
    WPM = math.floor(gentle_wordcount / (gentle_length / 60))
    results["WPM"] = WPM
    
    # Pause counts and average pause length
    sum_pauses = 0
    min_pause = 0.1
    max_pause = 3.0
    pause_count = 0
    long_pause_count = 0
    
    for x in range(len(gentle_end) - 1):
        tmp = gentle_start[x + 1] - gentle_end[x]
        if min_pause <= tmp <= max_pause:
            sum_pauses += tmp
            pause_count += 1
        elif tmp > max_pause:
            long_pause_count += 1
    
    results["pause_count"] = pause_count
    results["long_pause_count"] = long_pause_count
    
    if pause_count > 0:
        APL = decimal.Decimal(sum_pauses / pause_count)
        results["average_pause_length"] = float(round(APL, 2))
        APR = decimal.Decimal(pause_count / gentle_length)
        results["average_pause_rate"] = float(round(APR, 3))
    else:
        results["average_pause_length"] = 0.0
        results["average_pause_rate"] = 0.0
    
    # Rhythmic Complexity of Pauses
    s = []
    m = decimal.Decimal(str(gentle_start[0]))
    
    for x in range(len(gentle_end)):
        if x != len(gentle_end) - 1:
            start = decimal.Decimal(str(gentle_start[x]))
            next_start = decimal.Decimal(str(gentle_start[x + 1]))
            end = decimal.Decimal(str(gentle_end[x]))
            pause_length = decimal.Decimal(gentle_start[x + 1] - gentle_end[x])
            
            # Sampled every 10 ms
            while m >= start and m <= end:  # voiced
                s.append(1)
                m += decimal.Decimal('.01')
            
            while m > end and m < next_start:
                if min_pause <= pause_length <= max_pause:
                    s.append(0)
                else:
                    s.append(1)
                m += decimal.Decimal('.01')
    
    # Handle last segment
    if len(gentle_end) > 0:
        x = len(gentle_end) - 1
        start = decimal.Decimal(str(gentle_start[x]))
        end = decimal.Decimal(str(gentle_end[x]))
        while m >= start and m <= end:
            s.append(1)
            m += decimal.Decimal('.01')
    
    # Normalized LZ complexity
    if len(s) > 0:
        CP = lempel_ziv_complexity("".join([str(i) for i in s])) / (len(s) / math.log2(len(s)))
        results["rhythmic_complexity_of_pauses"] = CP * 100
    else:
        results["rhythmic_complexity_of_pauses"] = float('nan')
    
    # DRIFT: Pitch analysis
    drift_time = []
    drift_pitch = []
    
    for item in pitch_data:
        time = item.get('time')
        freq = item.get('frequency')
        
        if time is None or freq is None:
            continue
        if start_time and time < start_time:
            continue
        if end_time and time > end_time:
            continue
        if freq != 0:
            drift_time.append(float(time))
            drift_pitch.append(float(freq))
    
    if len(drift_pitch) == 0:
        results["average_pitch"] = float('nan')
        results["pitch_range"] = float('nan')
        results["pitch_speed"] = float('nan')
        results["pitch_acceleration"] = float('nan')
        results["pitch_entropy"] = float('nan')
        return results
    
    # Average Pitch
    results["average_pitch"] = sum(drift_pitch) / len(drift_pitch)
    
    # Pitch calculations
    f0log = [math.log2(p) for p in drift_pitch]
    f0mean = 2.0 ** (sum(f0log) / len(f0log))
    diffoctf0 = [math.log2(p) - math.log2(f0mean) for p in drift_pitch]
    
    # Pitch Range
    results["pitch_range"] = max(diffoctf0) - min(diffoctf0)
    
    # Pitch histogram and entropy
    f0hist, _ = np.histogram(diffoctf0, 25, (-1, 1))
    f0prob = [f / f0hist.sum() for f in f0hist]
    f0log2prob = [math.log2(f) if f != 0 else 0 for f in f0prob]
    
    f0entropy = -sum(f0prob[i] * f0log2prob[i] for i in range(len(f0prob)))
    results["pitch_entropy"] = f0entropy
    
    # Pitch velocity and acceleration
    # Find voiced segments
    dminvoice = 0.100
    ts = drift_time[1] - drift_time[0] if len(drift_time) > 1 else 0.01
    vdurthresh = round(dminvoice / ts)
    
    # Find contiguous regions (simplified - assumes continuous already)
    ixvoicedbounds = [(0, len(drift_pitch) - 1)]
    
    f0velocity = []
    f0accel = []
    
    for start, end in ixvoicedbounds:
        diffocttmp = np.array([diffoctf0[j] for j in range(start, end + 1)])
        
        if len(diffocttmp) < 3:
            continue
        
        # Apply Savitzky-Golay smoothing
        if len(diffocttmp) >= 7:
            diffocttmp = savgol_filter(diffocttmp, window_length=7, polyorder=2)
        
        velocity = np.diff(diffocttmp) / ts
        f0velocity.extend(velocity)
        
        accel = np.diff(np.diff(diffocttmp)) / ts
        f0accel.extend(accel)
    
    if len(f0velocity) > 0:
        f0velocity = np.array(f0velocity)
        results["pitch_speed"] = float(np.mean(np.abs(f0velocity)) * np.sign(np.mean(f0velocity)))
    else:
        results["pitch_speed"] = 0.0
    
    if len(f0accel) > 0:
        f0accel = np.array(f0accel)
        results["pitch_acceleration"] = float(np.mean(np.abs(f0accel)) * np.sign(np.mean(f0accel)))
    else:
        results["pitch_acceleration"] = 0.0
    
    return results

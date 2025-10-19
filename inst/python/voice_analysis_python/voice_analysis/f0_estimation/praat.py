"""
Praat-style F0 Estimation using Autocorrelation

Ported from F0_Thanasis.m (MATLAB Voice Analysis Toolbox)
Implements Boersma's autocorrelation method as used in Praat
"""

import numpy as np
from scipy import signal


def estimate_f0_praat(audio, fs, f0_min=50, f0_max=500, 
                      frame_length=0.04, frame_shift=0.01):
    """
    Praat-style F0 estimation using autocorrelation
    
    Ported from F0_Thanasis.m (lines 74-150)
    
    Parameters:
    -----------
    audio : ndarray
        Audio signal
    fs : int
        Sampling frequency (Hz)
    f0_min : float
        Minimum F0 (Hz), default 50
    f0_max : float
        Maximum F0 (Hz), default 500
    frame_length : float
        Analysis window length in seconds, default 0.04
    frame_shift : float
        Frame shift in seconds, default 0.01
        
    Returns:
    --------
    F0 : ndarray
        Fundamental frequency contour (Hz)
    """
    audio = np.asarray(audio).ravel()
    
    x_len = int(frame_length * fs)
    t_step = int(frame_shift * fs)
    n_frames = (len(audio) - x_len) // t_step
    
    if n_frames <= 0:
        return np.array([])
    
    F0 = np.zeros(n_frames)
    
    # Gaussian window (matching MATLAB gausswin with std = length/6)
    window = signal.windows.gaussian(x_len, std=x_len/6)
    
    low_lim = int(np.ceil(fs / f0_max))
    up_lim = int(np.floor(fs / f0_min))
    
    for i in range(n_frames):
        # Extract and center frame
        start_idx = i * t_step
        end_idx = start_idx + x_len
        
        if end_idx > len(audio):
            break
            
        frame = audio[start_idx:end_idx]
        frame = frame - np.mean(frame)
        windowed = frame * window
        
        # Autocorrelation using FFT (Boersma's method)
        fft_signal = np.fft.fft(windowed)
        acf = np.fft.ifft(np.abs(fft_signal)**2).real
        
        # Window autocorrelation
        fft_window = np.fft.fft(window)
        acf_window = np.fft.ifft(np.abs(fft_window)**2).real
        
        # Normalized autocorrelation
        with np.errstate(divide='ignore', invalid='ignore'):
            rho = acf / acf_window
            rho = np.nan_to_num(rho, nan=0.0, posinf=0.0, neginf=0.0)
        
        rho = rho[:len(rho)//2]
        
        if np.max(np.abs(rho)) > 0:
            rho = rho / np.max(np.abs(rho))
        
        # Find peak in valid F0 range
        if up_lim > low_lim and up_lim <= len(rho):
            valid_rho = rho[low_lim:up_lim]
            
            if len(valid_rho) > 0 and np.max(valid_rho) > 0:
                peak_idx = np.argmax(valid_rho) + low_lim
                
                # Parabolic interpolation (lines 135-143 in F0_Thanasis.m)
                if peak_idx > low_lim and peak_idx < up_lim - 1:
                    alpha = rho[peak_idx - 1]
                    beta = rho[peak_idx]
                    gamma = rho[peak_idx + 1]
                    
                    denominator = 2*beta - alpha - gamma
                    
                    if abs(denominator) > 1e-10:
                        dt = 1.0 / fs
                        peak_refined = peak_idx + 0.5 * (alpha - gamma) / denominator
                        tmax = dt * peak_refined
                        
                        if tmax > 0:
                            F0[i] = 1.0 / tmax
                    else:
                        F0[i] = fs / peak_idx
                else:
                    F0[i] = fs / peak_idx
    
    return F0

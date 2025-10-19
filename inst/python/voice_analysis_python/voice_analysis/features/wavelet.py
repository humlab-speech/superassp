"""
Wavelet Features

Ported from wavedec_features.m
Uses discrete wavelet decomposition to extract features
"""

import numpy as np
import warnings
from ..utils.tkeo import compute_tkeo_vectorized
from ..utils.entropy import compute_entropy


def compute_wavelet_features(signal, wavelet='db8', level=10):
    """
    Compute wavelet decomposition features
    
    Ported from wavedec_features.m
    
    Parameters:
    -----------
    signal : ndarray
        Input signal (typically F0 contour)
    wavelet : str
        Wavelet family (default 'db8')
    level : int
        Decomposition level (default 10)
        
    Returns:
    --------
    measures : dict
        Wavelet-based features
    """
    try:
        import pywt
        
        signal = np.asarray(signal).ravel()
        
        # Remove invalid values
        signal = signal[np.isfinite(signal)]
        signal = signal[signal > 0]
        
        if len(signal) < 2**level:
            warnings.warn(f"Signal too short for {level} levels, reducing to {int(np.log2(len(signal)))}")
            level = max(1, int(np.log2(len(signal))) - 1)
        
        measures = {}
        
        # Decomposition of original signal
        coeffs = pywt.wavedec(signal, wavelet, level=level)
        
        # Approximation coefficients (low frequency)
        cA = coeffs[0]
        
        # Detail coefficients (high frequency)
        cDs = coeffs[1:]
        
        # Energy
        energy_approx = np.sum(cA**2)
        energy_detail = sum([np.sum(cD**2) for cD in cDs])
        
        measures['wavelet_energy_approx'] = energy_approx
        measures['wavelet_energy_detail'] = energy_detail
        
        # Features from detail coefficients
        for i, cD in enumerate(cDs, 1):
            if len(cD) > 0:
                # Shannon entropy
                hist, _ = np.histogram(cD, bins=50)
                measures[f'wavelet_det{i}_entropy_shannon'] = compute_entropy(hist, base='2')
                
                # Log energy entropy
                log_energy = np.log(np.abs(cD) + 1e-10)
                measures[f'wavelet_det{i}_entropy_log'] = np.std(log_energy)
                
                # TKEO features
                tkeo = compute_tkeo_vectorized(cD)
                measures[f'wavelet_det{i}_TKEO_mean'] = np.mean(np.abs(tkeo))
                measures[f'wavelet_det{i}_TKEO_std'] = np.std(tkeo)
        
        # Features from approximation coefficients
        if len(cA) > 0:
            hist, _ = np.histogram(cA, bins=50)
            measures['wavelet_app_entropy_shannon'] = compute_entropy(hist, base='2')
            
            log_energy = np.log(np.abs(cA) + 1e-10)
            measures['wavelet_app_entropy_log'] = np.std(log_energy)
            
            tkeo = compute_tkeo_vectorized(cA)
            measures['wavelet_app_TKEO_mean'] = np.mean(np.abs(tkeo))
            measures['wavelet_app_TKEO_std'] = np.std(tkeo)
        
        # Repeat for log-transformed signal
        log_signal = np.log(signal + 1e-10)
        coeffs_log = pywt.wavedec(log_signal, wavelet, level=level)
        
        cA_log = coeffs_log[0]
        cDs_log = coeffs_log[1:]
        
        # Energy of log-transformed
        energy_approx_log = np.sum(cA_log**2)
        energy_detail_log = sum([np.sum(cD**2) for cD in cDs_log])
        
        measures['wavelet_log_energy_approx'] = energy_approx_log
        measures['wavelet_log_energy_detail'] = energy_detail_log
        
        # Features from log-transformed detail coefficients
        for i, cD in enumerate(cDs_log, 1):
            if len(cD) > 0:
                hist, _ = np.histogram(cD, bins=50)
                measures[f'wavelet_log_det{i}_entropy_shannon'] = compute_entropy(hist, base='2')
                
                log_energy = np.log(np.abs(cD) + 1e-10)
                measures[f'wavelet_log_det{i}_entropy_log'] = np.std(log_energy)
                
                tkeo = compute_tkeo_vectorized(cD)
                measures[f'wavelet_log_det{i}_TKEO_mean'] = np.mean(np.abs(tkeo))
                measures[f'wavelet_log_det{i}_TKEO_std'] = np.std(tkeo)
        
        # Features from log-transformed approximation
        if len(cA_log) > 0:
            hist, _ = np.histogram(cA_log, bins=50)
            measures['wavelet_log_app_entropy_shannon'] = compute_entropy(hist, base='2')
            
            log_energy = np.log(np.abs(cA_log) + 1e-10)
            measures['wavelet_log_app_entropy_log'] = np.std(log_energy)
            
            tkeo = compute_tkeo_vectorized(cA_log)
            measures['wavelet_log_app_TKEO_mean'] = np.mean(np.abs(tkeo))
            measures['wavelet_log_app_TKEO_std'] = np.std(tkeo)
        
        return measures
        
    except ImportError:
        warnings.warn("PyWavelets (pywt) not available. Install with: pip install PyWavelets")
        return {'wavelet_error': np.nan}

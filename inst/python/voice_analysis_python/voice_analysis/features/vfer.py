"""
Vocal Fold Excitation Ratios (VFER)

Ported from VFER_measure function in voice_analysis_redux.m (lines 508-568)
Requires DYPSA for glottal closure instant detection
"""

import numpy as np
from scipy import signal
import warnings
from ..utils.tkeo import compute_tkeo_vectorized
from ..utils.entropy import safe_log


def compute_vfer(audio, fs):
    """
    Compute Vocal Fold Excitation Ratios
    
    Requires DYPSA algorithm for GCI detection
    
    Parameters:
    -----------
    audio : ndarray
        Audio signal
    fs : int
        Sampling frequency
        
    Returns:
    --------
    measures : dict
        7 VFER-related measures
    """
    try:
        from ..utils.dypsa import dypsa
        
        audio = np.asarray(audio).ravel()
        
        # Parameters from MATLAB (lines 510-513)
        filt_order = 100
        BW = 500  # Bandwidth (Hz)
        Fmax = fs / 2 - BW - 300  # Max frequency
        Fshift = 500  # Frequency shift (Hz)
        
        # Run DYPSA to get GCI
        gci, _ = dypsa(audio, fs)
        
        if len(gci) < 2:
            warnings.warn("Insufficient glottal closure instants detected")
            return _get_nan_vfer()
        
        # Create filter bank (lines 518-526)
        Fc1 = np.arange(1, Fmax, Fshift)  # Lower cutoff frequencies
        Fc2 = Fc1 + BW  # Upper cutoff frequencies
        
        filters = []
        for fc1, fc2 in zip(Fc1, Fc2):
            try:
                sos = signal.butter(
                    filt_order // 2,
                    [fc1, fc2],
                    btype='band',
                    fs=fs,
                    output='sos'
                )
                filters.append(sos)
            except:
                continue
        
        if len(filters) == 0:
            return _get_nan_vfer()
        
        # Process each glottal cycle (lines 528-549)
        NEm = []
        signal_BW_TKEO_all = []
        signal_BW_SEO_all = []
        
        for i in range(len(gci) - 1):
            start_idx = int(gci[i])
            end_idx = int(gci[i + 1])
            
            if start_idx >= end_idx or end_idx > len(audio):
                continue
            
            # Extract cycle
            cycle = audio[start_idx:end_idx]
            
            # Hanning window
            window = signal.windows.hann(len(cycle))
            windowed = cycle * window
            
            # Process through filter bank (lines 534-547)
            if len(windowed) > 50:
                sigBW = np.zeros((min(50, len(windowed)), len(filters)))
                sig_TKEO = np.zeros(len(filters))
                sig_SEO = np.zeros(len(filters))
                
                for j, filt in enumerate(filters):
                    try:
                        filtered = signal.sosfilt(filt, windowed)
                        # Use first 50 samples (line 536)
                        sigBW[:, j] = filtered[:min(50, len(filtered))]
                        
                        # TKEO (line 537)
                        sig_TKEO[j] = np.mean(compute_tkeo_vectorized(sigBW[:, j]))
                        
                        # Signal energy (line 538)
                        sig_SEO[j] = np.mean(sigBW[:, j])**2
                    except:
                        pass
                
                # Hilbert envelope cross-correlation (lines 540-544)
                try:
                    Hilb_tr = signal.hilbert(sigBW, axis=0)
                    Hilb_env = np.abs(Hilb_tr)
                    c = np.correlate(Hilb_env.T.ravel(), Hilb_env.T.ravel())
                    NEm.append(np.max(c))
                except:
                    pass
                
                signal_BW_TKEO_all.append(sig_TKEO)
                signal_BW_SEO_all.append(sig_SEO)
        
        if len(NEm) == 0 or len(signal_BW_TKEO_all) == 0:
            return _get_nan_vfer()
        
        # Compute statistics (lines 551-567)
        signal_BW_TKEO = np.array(signal_BW_TKEO_all)
        signal_BW_SEO = np.array(signal_BW_SEO_all)
        
        # Log transform (line 551)
        signal_BW_TKEO_log = np.mean(safe_log(signal_BW_TKEO + 1e-10, base='e'), axis=0)
        
        measures = {}
        
        # Line 555: Mean of correlation
        measures['VFER_mean'] = np.mean(NEm)
        
        # Line 556: Std of correlation
        measures['VFER_std'] = np.std(NEm)
        
        # Line 557: Entropy
        measures['VFER_entropy'] = -np.sum(NEm * safe_log(np.array(NEm) + 1e-10, base='e'))
        
        # Mean across cycles
        VFTKEO = np.mean(signal_BW_TKEO, axis=0)
        VFSEO = np.mean(signal_BW_SEO, axis=0)
        VFlog_SEO = np.mean(safe_log(signal_BW_SEO + 1e-10, base='e'), axis=0)
        
        # Signal-to-noise ratios (lines 563-566)
        n_bands = len(VFTKEO)
        
        if n_bands >= 10:
            # Line 563: Low frequency (1:5) / High frequency (6:10) TKEO ratio
            measures['VFER_TKEO_low_high'] = (
                np.sum(VFTKEO[:5]) / (np.sum(VFTKEO[5:10]) + 1e-10)
            )
            
            # Line 564: Low frequency (1:5) / High frequency (6:10) SEO ratio
            measures['VFER_SEO_low_high'] = (
                np.sum(VFSEO[:5]) / (np.sum(VFSEO[5:10]) + 1e-10)
            )
            
            # Line 565: High frequency (6:10) / Low frequency (1:5) log TKEO ratio
            measures['VFER_log_TKEO_high_low'] = (
                np.sum(signal_BW_TKEO_log[5:10]) / (np.sum(signal_BW_TKEO_log[:5]) + 1e-10)
            )
            
            # Line 566: High frequency (6:10) / Low frequency (1:5) log SEO ratio
            measures['VFER_log_SEO_high_low'] = (
                np.sum(VFlog_SEO[5:10]) / (np.sum(VFlog_SEO[:5]) + 1e-10)
            )
        else:
            # Fallback if not enough frequency bands
            measures['VFER_TKEO_low_high'] = np.nan
            measures['VFER_SEO_low_high'] = np.nan
            measures['VFER_log_TKEO_high_low'] = np.nan
            measures['VFER_log_SEO_high_low'] = np.nan
        
        return measures
        
    except ImportError:
        warnings.warn("DYPSA module not available")
        return _get_nan_vfer()
    except Exception as e:
        warnings.warn(f"VFER computation failed: {e}")
        return _get_nan_vfer()


def _get_nan_vfer():
    """Return NaN values for all VFER measures"""
    return {
        'VFER_mean': np.nan,
        'VFER_std': np.nan,
        'VFER_entropy': np.nan,
        'VFER_TKEO_low_high': np.nan,
        'VFER_SEO_low_high': np.nan,
        'VFER_log_TKEO_high_low': np.nan,
        'VFER_log_SEO_high_low': np.nan,
    }

"""
Empirical Mode Decomposition (EMD) Features

Ported from IMF_measure function in voice_analysis_redux.m (lines 570-596)
Uses PyEMD library for EMD decomposition

Reference:
- Flandrin, P., Rilling, G., & Goncalves, P. (2004).
  Empirical mode decomposition as a filter bank.
  IEEE Signal Processing Letters, 11(2), 112-114.
"""

import numpy as np
import warnings
from ..utils.tkeo import compute_tkeo_vectorized
from ..utils.entropy import safe_log


def compute_emd_features(signal):
    """
    Compute EMD-based features
    
    Ported from IMF_measure function (lines 570-596)
    
    Parameters:
    -----------
    signal : ndarray
        Input signal (audio)
        
    Returns:
    --------
    measures : dict
        6 EMD-based ratio measures
    """
    try:
        # Import from PyEMD (EMD-signal package, installed as pyemd but needs PyEMD import)
        from PyEMD.EMD import EMD
        
        signal = np.asarray(signal).ravel()
        
        # Remove invalid values
        signal = signal[np.isfinite(signal)]
        
        # EMD is computationally expensive - downsample long signals
        # Use max 20000 samples for computational efficiency
        if len(signal) > 20000:
            from scipy import signal as scipy_signal
            downsample_factor = len(signal) // 20000 + 1
            signal = scipy_signal.decimate(signal, downsample_factor, ftype='fir', zero_phase=True)
        
        if len(signal) < 100:
            return _get_nan_emd()
        
        # Perform EMD decomposition
        emd = EMD()
        
        try:
            IMFs = emd(signal)
        except Exception as e:
            warnings.warn(f"EMD decomposition failed: {e}")
            return _get_nan_emd()
        
        # IMFs is 2D array: (n_imfs, n_samples)
        if IMFs.ndim == 1:
            IMFs = IMFs.reshape(1, -1)
        
        n_imfs = IMFs.shape[0]
        
        if n_imfs < 4:
            warnings.warn(f"Too few IMFs ({n_imfs}), need at least 4")
            return _get_nan_emd()
        
        # Compute features for each IMF (lines 579-586)
        IMF_energy = []
        IMF_TKEO = []
        IMF_entropy = []
        
        # Log-transformed IMFs
        IMFs_log = safe_log(np.abs(IMFs) + 1e-10, base='e')
        
        IMF_energy_log = []
        IMF_TKEO_log = []
        IMF_entropy_log = []
        
        for i in range(n_imfs):
            imf = IMFs[i, :]
            imf_log = IMFs_log[i, :]
            
            # Original IMF features
            IMF_energy.append(np.abs(np.mean(imf**2)))
            IMF_TKEO.append(np.abs(np.mean(compute_tkeo_vectorized(imf))))
            
            # Entropy (line 582)
            entropy_val = np.abs(np.mean(-np.sum(imf * safe_log(np.abs(imf) + 1e-10, base='e'))))
            IMF_entropy.append(entropy_val)
            
            # Log-transformed IMF features
            IMF_energy_log.append(np.abs(np.mean(imf_log**2)))
            IMF_TKEO_log.append(np.abs(np.mean(compute_tkeo_vectorized(imf_log))))
            
            entropy_log_val = np.abs(np.mean(-np.sum(imf_log * safe_log(np.abs(imf_log) + 1e-10, base='e'))))
            IMF_entropy_log.append(entropy_log_val)
        
        # Compute ratio measures (lines 588-595)
        measures = {}
        
        # Ratios of high-frequency IMFs (4:end) to low-frequency IMFs (1:3)
        # Line 589
        if n_imfs >= 4:
            measures['EMD_energy_ratio'] = (
                np.sum(IMF_energy[3:]) / (np.sum(IMF_energy[:3]) + 1e-10)
            )
        else:
            measures['EMD_energy_ratio'] = np.nan
        
        # Line 590
        if n_imfs >= 4:
            measures['EMD_TKEO_ratio'] = (
                np.sum(IMF_TKEO[3:]) / (np.sum(IMF_TKEO[:3]) + 1e-10)
            )
        else:
            measures['EMD_TKEO_ratio'] = np.nan
        
        # Line 591
        if n_imfs >= 4:
            measures['EMD_entropy_ratio'] = (
                np.sum(IMF_entropy[3:]) / (np.sum(IMF_entropy[:3]) + 1e-10)
            )
        else:
            measures['EMD_entropy_ratio'] = np.nan
        
        # Log-transformed ratios (lines 592-595)
        # Note: these are inverted (low/high instead of high/low)
        
        # Line 592
        if n_imfs >= 4:
            measures['EMD_log_energy_ratio'] = np.abs(
                np.sum(IMF_energy_log[:2]) / (np.sum(IMF_energy_log[3:]) + 1e-10)
            )
        else:
            measures['EMD_log_energy_ratio'] = np.nan
        
        # Line 593
        if n_imfs >= 4:
            measures['EMD_log_TKEO_ratio'] = np.abs(
                np.sum(IMF_TKEO_log[:2]) / (np.sum(IMF_TKEO_log[2:]) + 1e-10)
            )
        else:
            measures['EMD_log_TKEO_ratio'] = np.nan
        
        # Line 594
        if n_imfs >= 4:
            measures['EMD_log_entropy_ratio'] = (
                np.sum(IMF_entropy_log[:2]) / (np.sum(IMF_entropy_log[2:]) + 1e-10)
            )
        else:
            measures['EMD_log_entropy_ratio'] = np.nan
        
        return measures
        
    except ImportError:
        warnings.warn(
            "PyEMD (EMD-signal) not available. Install with: pip install EMD-signal\n"
            "Note: After installation, the package folder may need to be renamed from 'pyemd' to 'PyEMD'.\n"
            "Returning NaN for EMD features."
        )
        return _get_nan_emd()
    except Exception as e:
        warnings.warn(f"EMD feature computation failed: {e}")
        return _get_nan_emd()


def _get_nan_emd():
    """Return NaN values for all EMD measures"""
    return {
        'EMD_energy_ratio': np.nan,
        'EMD_TKEO_ratio': np.nan,
        'EMD_entropy_ratio': np.nan,
        'EMD_log_energy_ratio': np.nan,
        'EMD_log_TKEO_ratio': np.nan,
        'EMD_log_entropy_ratio': np.nan,
    }

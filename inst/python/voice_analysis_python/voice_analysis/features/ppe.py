"""
Pitch Period Entropy (PPE)

Ported from voice_analysis_redux.m (lines 203-211)
"""

import numpy as np
from scipy import signal
from ..utils.entropy import compute_entropy, safe_log
import warnings


def compute_ppe(F0, fs, f0_mean_healthy=120, use_thesis_mode=False):
    """
    Compute Pitch Period Entropy
    
    PPE quantifies the impaired control of pitch in voice disorders
    
    Parameters:
    -----------
    F0 : ndarray
        Fundamental frequency contour
    fs : int
        Sampling frequency
    f0_mean_healthy : float
        Mean F0 for healthy controls (120 Hz for males, 190 Hz for females)
        Used in MATLAB mode only
    use_thesis_mode : bool
        If True, use thesis specification (Equation 3.54: semitone scale)
        If False (default), use MATLAB implementation (natural log)
        
    Returns:
    --------
    ppe : float
        Pitch Period Entropy value
    """
    F0 = np.asarray(F0).ravel()
    
    # Remove zeros and invalid values
    F0 = F0[F0 > 0]
    F0 = F0[np.isfinite(F0)]
    
    if len(F0) < 20:
        return np.nan
    
    try:
        # Transform F0 to log scale
        if use_thesis_mode:
            # Thesis specification (Equation 3.54, p. 71)
            # s = 12 × log₂(F₀/127) [semitone scale]
            # Reference: 127 Hz for males (footnote 22, p. 71)
            logF0 = 12 * np.log2(F0 / 127)
        else:
            # MATLAB implementation (natural log, normalized)
            logF0 = safe_log(F0 / f0_mean_healthy, base='e')
        
        # Estimate AR(10) coefficients using autocorrelation
        # This matches MATLAB's arcov function behavior
        p = 10  # AR order
        
        if len(logF0) <= p:
            return np.nan
        
        # Center the signal
        logF0_centered = logF0 - np.mean(logF0)
        
        # Compute autocorrelation
        acf = np.correlate(logF0_centered, logF0_centered, mode='full')
        acf = acf[len(acf)//2:]
        
        if acf[0] == 0:
            return np.nan
        
        acf = acf / acf[0]  # Normalize
        
        # Levinson-Durbin algorithm (simplified)
        r = acf[1:p+1]
        
        # Build Toeplitz matrix
        R = np.zeros((p, p))
        for i in range(p):
            for j in range(p):
                if i >= j:
                    R[i, j] = acf[abs(i - j)]
                else:
                    R[i, j] = acf[abs(j - i)]
        
        if np.linalg.det(R) == 0:
            return np.nan
        
        # Solve for AR coefficients
        ar_coeffs = np.linalg.solve(R, r)
        
        # Apply inverse filter
        # MATLAB: sig_filtered = filter(ARcoef, 1, logF0signal);
        # Note: MATLAB's filter with [a0, a1, ..., an] is different from scipy
        # We need to prepend 1.0 to match MATLAB's behavior
        a_coeffs = np.concatenate([[1.0], -ar_coeffs])
        
        sig_filtered = signal.lfilter([1.0], a_coeffs, logF0)
        
        # Remove initial samples (line 208: round(0.001*fs) samples)
        # But logF0 is in F0 time, not audio samples
        # Skip first few samples to avoid transients
        skip_samples = max(1, len(sig_filtered) // 100)
        sig_filtered = sig_filtered[skip_samples:]
        
        if len(sig_filtered) == 0:
            return np.nan
        
        # Create histogram (100 bins)
        hist, bin_edges = np.histogram(
            sig_filtered,
            bins=100,
            range=(np.min(sig_filtered), np.max(sig_filtered))
        )
        
        # Compute entropy
        if np.sum(hist) == 0:
            return np.nan
        
        PPEd = hist / np.sum(hist)
        
        # Entropy normalized by log of number of bins
        H = compute_entropy(PPEd, base='e')
        PPE = H / safe_log(len(PPEd), base='e')
        
        return PPE
        
    except Exception as e:
        warnings.warn(f"PPE computation failed: {e}")
        return np.nan

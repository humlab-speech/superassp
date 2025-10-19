"""
SWIPE F0 Estimation using pySPTK

SWIPE: A Sawtooth Waveform Inspired Pitch Estimator
Reference: Camacho, A., & Harris, J.G. (2008). A sawtooth waveform inspired 
pitch estimator for speech and music. JASA, 124(3), 1638-1652.
"""

import numpy as np
import warnings

def estimate_f0_swipe(audio, fs, f0_min=50, f0_max=500, frame_shift=0.01):
    """
    SWIPE F0 estimation using pySPTK
    
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
    frame_shift : float
        Frame shift in seconds, default 0.01
        
    Returns:
    --------
    F0 : ndarray
        Fundamental frequency contour (Hz)
    """
    try:
        import pysptk
        
        # SWIPE parameters
        hop_length = int(frame_shift * fs)
        
        # Run SWIPE
        F0 = pysptk.sptk.swipe(
            audio.astype(np.float64),
            fs=fs,
            hopsize=hop_length,
            min=f0_min,
            max=f0_max,
            otype="f0"  # Output F0 in Hz
        )
        
        return F0
        
    except ImportError:
        warnings.warn(
            "pySPTK not available. Install with: pip install pysptk\n"
            "Falling back to Praat-style F0 estimation."
        )
        from .praat import estimate_f0_praat
        return estimate_f0_praat(audio, fs, f0_min, f0_max, frame_shift=frame_shift)
    except Exception as e:
        warnings.warn(f"SWIPE estimation failed: {e}. Falling back to Praat-style.")
        from .praat import estimate_f0_praat
        return estimate_f0_praat(audio, fs, f0_min, f0_max, frame_shift=frame_shift)

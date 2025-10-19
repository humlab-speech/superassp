"""
REAPER F0 estimation (Google's Robust Epoch And Pitch EstimatoR)
"""

import numpy as np
from typing import Tuple
from .base import F0Estimator

try:
    import pyreaper
    REAPER_AVAILABLE = True
except ImportError:
    REAPER_AVAILABLE = False


class ReaperF0Estimator(F0Estimator):
    """F0 estimation using REAPER algorithm"""
    
    def __init__(self):
        if not REAPER_AVAILABLE:
            raise ImportError("pyreaper is not installed. Install with: pip install pyreaper")
    
    def estimate(self, audio: np.ndarray, sr: int, 
                frame_shift: float = 1.0,
                f0_min: float = 40.0,
                f0_max: float = 500.0,
                **kwargs) -> Tuple[np.ndarray, np.ndarray]:
        """
        Estimate F0 using REAPER
        
        Args:
            audio: Audio signal
            sr: Sampling rate
            frame_shift: Frame shift in milliseconds
            f0_min: Minimum F0
            f0_max: Maximum F0
        
        Returns:
            f0_values: F0 values (NaN for unvoiced)
            times: Time points in seconds
        """
        # REAPER requires int16
        if audio.dtype != np.int16:
            # Convert float to int16
            if audio.dtype in [np.float32, np.float64]:
                audio = (audio * 32767).astype(np.int16)
            else:
                audio = audio.astype(np.int16)
        
        # REAPER expects frame period in seconds
        frame_period = frame_shift / 1000.0
        
        # Run REAPER
        pm_times, pm_values, f0_times, f0_values, corr = pyreaper.reaper(
            audio,
            sr,
            minf0=f0_min,
            maxf0=f0_max,
            frame_period=frame_period,
            do_hilbert_transform=True,
            do_high_pass=True
        )
        
        # Clean up F0 values
        f0_clean = self._remove_zeros(f0_values)
        f0_clean = self._validate_f0_range(f0_clean, f0_min, f0_max)
        
        return f0_clean, f0_times

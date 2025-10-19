"""
WORLD vocoder F0 estimation
"""

import numpy as np
from typing import Tuple
from .base import F0Estimator

try:
    import pyworld as pw
    WORLD_AVAILABLE = True
except ImportError:
    WORLD_AVAILABLE = False


class WorldF0Estimator(F0Estimator):
    """F0 estimation using WORLD vocoder"""
    
    def __init__(self):
        if not WORLD_AVAILABLE:
            raise ImportError("pyworld is not installed. Install with: pip install pyworld")
    
    def estimate(self, audio: np.ndarray, sr: int,
                frame_shift: float = 1.0,
                f0_min: float = 40.0,
                f0_max: float = 500.0,
                **kwargs) -> Tuple[np.ndarray, np.ndarray]:
        """
        Estimate F0 using WORLD
        
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
        # WORLD requires float64
        if audio.dtype != np.float64:
            audio = audio.astype(np.float64)
        
        # Run WORLD DIO algorithm
        frame_period = frame_shift
        _f0, t = pw.dio(audio, sr, f0_floor=f0_min, f0_ceil=f0_max, 
                       frame_period=frame_period)
        
        # Refine F0 using Stonemask
        f0_values = pw.stonemask(audio, _f0, t, sr)
        times = t
        
        # Clean up F0 values
        f0_clean = self._remove_zeros(f0_values)
        f0_clean = self._validate_f0_range(f0_clean, f0_min, f0_max)
        
        return f0_clean, times

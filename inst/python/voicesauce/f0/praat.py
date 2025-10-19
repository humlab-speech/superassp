"""
Praat F0 estimation via Parselmouth
"""

import numpy as np
from typing import Tuple
from .base import F0Estimator

try:
    import parselmouth
    PARSELMOUTH_AVAILABLE = True
except ImportError:
    PARSELMOUTH_AVAILABLE = False


class PraatF0Estimator(F0Estimator):
    """F0 estimation using Praat (via Parselmouth)"""
    
    def __init__(self):
        if not PARSELMOUTH_AVAILABLE:
            raise ImportError("parselmouth is not installed. Install with: pip install praat-parselmouth")
    
    def estimate(self, audio: np.ndarray, sr: int,
                frame_shift: float = 1.0,
                f0_min: float = 40.0,
                f0_max: float = 500.0,
                silence_threshold: float = 0.03,
                voicing_threshold: float = 0.45,
                octave_cost: float = 0.01,
                octave_jump_cost: float = 0.35,
                voiced_unvoiced_cost: float = 0.14,
                **kwargs) -> Tuple[np.ndarray, np.ndarray]:
        """
        Estimate F0 using Praat
        
        Args:
            audio: Audio signal
            sr: Sampling rate
            frame_shift: Frame shift in milliseconds
            f0_min: Minimum F0
            f0_max: Maximum F0
            silence_threshold: Silence threshold
            voicing_threshold: Voicing threshold
            octave_cost: Octave cost
            octave_jump_cost: Octave jump cost
            voiced_unvoiced_cost: Voiced/unvoiced cost
        
        Returns:
            f0_values: F0 values (NaN for unvoiced)
            times: Time points in seconds
        """
        # Create Praat Sound object
        sound = parselmouth.Sound(audio, sampling_frequency=sr)
        
        # Extract pitch using autocorrelation (can also use 'cc' for cross-correlation)
        pitch = sound.to_pitch(
            time_step=frame_shift / 1000.0,
            pitch_floor=f0_min,
            pitch_ceiling=f0_max
        )
        
        # Extract F0 values
        times = pitch.xs()
        f0_values = np.array([pitch.get_value_at_time(t) for t in times])
        
        # Clean up F0 values
        f0_clean = self._remove_zeros(f0_values)
        f0_clean = self._validate_f0_range(f0_clean, f0_min, f0_max)
        
        return f0_clean, times

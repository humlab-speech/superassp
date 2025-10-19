"""
Base class for F0 estimation
"""

from abc import ABC, abstractmethod
import numpy as np
from typing import Tuple


class F0Estimator(ABC):
    """Abstract base class for F0 estimators"""
    
    @abstractmethod
    def estimate(self, audio: np.ndarray, sr: int, **kwargs) -> Tuple[np.ndarray, np.ndarray]:
        """
        Estimate F0 from audio
        
        Args:
            audio: Audio signal
            sr: Sampling rate
            **kwargs: Additional parameters
        
        Returns:
            f0_values: F0 values (NaN for unvoiced)
            times: Time points in seconds
        """
        pass
    
    @staticmethod
    def _remove_zeros(f0_values: np.ndarray) -> np.ndarray:
        """Replace zeros with NaN"""
        f0_clean = f0_values.copy()
        f0_clean[f0_clean == 0] = np.nan
        return f0_clean
    
    @staticmethod
    def _validate_f0_range(f0_values: np.ndarray, f0_min: float, f0_max: float) -> np.ndarray:
        """Validate F0 values are within range"""
        f0_clean = f0_values.copy()
        mask = (f0_clean < f0_min) | (f0_clean > f0_max)
        f0_clean[mask] = np.nan
        return f0_clean

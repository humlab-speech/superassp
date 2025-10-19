"""
F0 estimation module
"""

from .base import F0Estimator
from .reaper import ReaperF0Estimator, REAPER_AVAILABLE
from .praat import PraatF0Estimator, PARSELMOUTH_AVAILABLE
from .world import WorldF0Estimator, WORLD_AVAILABLE

import numpy as np
from typing import Tuple

__all__ = [
    'F0Estimator',
    'ReaperF0Estimator',
    'PraatF0Estimator',
    'WorldF0Estimator',
    'get_f0_estimator',
    'estimate_f0'
]


def get_f0_estimator(method: str) -> F0Estimator:
    """Get F0 estimator by method name"""
    method = method.lower()
    
    if method == 'reaper':
        if not REAPER_AVAILABLE:
            raise ImportError("pyreaper not available")
        return ReaperF0Estimator()
    elif method == 'praat':
        if not PARSELMOUTH_AVAILABLE:
            raise ImportError("parselmouth not available")
        return PraatF0Estimator()
    elif method == 'world':
        if not WORLD_AVAILABLE:
            raise ImportError("pyworld not available")
        return WorldF0Estimator()
    else:
        raise ValueError(f"Unknown F0 method: {method}")


def estimate_f0(audio: np.ndarray, sr: int, method: str = 'reaper', 
               **kwargs) -> Tuple[np.ndarray, np.ndarray]:
    """Convenience function to estimate F0"""
    estimator = get_f0_estimator(method)
    return estimator.estimate(audio, sr, **kwargs)

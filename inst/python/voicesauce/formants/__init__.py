"""Formant estimation module"""
from .praat import estimate_formants_praat, PARSELMOUTH_AVAILABLE
import numpy as np
from typing import Tuple

def estimate_formants(audio: np.ndarray, sr: int, method: str = 'praat',
                     **kwargs) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Estimate formants"""
    method = method.lower()
    if method == 'praat':
        if not PARSELMOUTH_AVAILABLE:
            raise ImportError("parselmouth not available")
        return estimate_formants_praat(audio, sr, **kwargs)
    else:
        raise ValueError(f"Unknown formant method: {method}")

"""
Praat formant estimation via Parselmouth
"""

import numpy as np
from typing import Tuple

try:
    import parselmouth
    PARSELMOUTH_AVAILABLE = True
except ImportError:
    PARSELMOUTH_AVAILABLE = False


def estimate_formants_praat(audio: np.ndarray, sr: int,
                            frame_shift: float = 1.0,
                            window_length: float = 0.025,
                            n_formants: int = 5,
                            max_formant: float = 5500.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Estimate formants using Praat
    
    Args:
        audio: Audio signal
        sr: Sampling rate
        frame_shift: Frame shift in milliseconds
        window_length: Window length in seconds
        n_formants: Number of formants to estimate
        max_formant: Maximum formant frequency
    
    Returns:
        formants: Formant frequencies (n_frames, n_formants)
        bandwidths: Formant bandwidths (n_frames, n_formants)
        times: Time points in seconds
    """
    if not PARSELMOUTH_AVAILABLE:
        raise ImportError("parselmouth is not installed")
    
    # Create Praat Sound object
    sound = parselmouth.Sound(audio, sampling_frequency=sr)
    
    # Extract formants using Burg's method
    formant = sound.to_formant_burg(
        time_step=frame_shift / 1000.0,
        max_number_of_formants=n_formants,
        maximum_formant=max_formant,
        window_length=window_length,
        pre_emphasis_from=50.0
    )
    
    # Extract times
    times = formant.xs()
    n_frames = len(times)
    
    # Initialize arrays
    formants = np.zeros((n_frames, n_formants))
    bandwidths = np.zeros((n_frames, n_formants))
    
    # Extract formant values
    for i, t in enumerate(times):
        for j in range(1, n_formants + 1):
            try:
                formants[i, j-1] = formant.get_value_at_time(j, t)
                bandwidths[i, j-1] = formant.get_bandwidth_at_time(j, t)
            except:
                formants[i, j-1] = np.nan
                bandwidths[i, j-1] = np.nan
    
    return formants, bandwidths, times

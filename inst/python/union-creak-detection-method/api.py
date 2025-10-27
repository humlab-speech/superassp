"""
High-level API for R integration.
Optimized for av package NumPy arrays - no file I/O.
"""

import numpy as np
from typing import Dict, Tuple, Optional

def detect_creak_union_extended(
    audio: np.ndarray,
    sample_rate: int = 16000,
    use_am: bool = True,
    use_cd: bool = True,
    cd_threshold: float = 0.3,
    use_reaper: bool = True,
    frame_shift_ms: float = 10.0,
    return_features: bool = False,
    return_probabilities: bool = False
) -> Dict[str, np.ndarray]:
    """
    Detect creaky voice with extended output for R integration.
    
    Returns all intermediate tracks including features and probabilities.
    
    Parameters
    ----------
    audio : np.ndarray
        Audio signal as 1D NumPy array from av::read_audio_bin()
    sample_rate : int
        Sampling rate in Hz (default 16000)
    use_am : bool
        Use AM (antimode) method
    use_cd : bool
        Use CD (creak detector) method with trained ANN
    cd_threshold : float
        Classification threshold for CD method
    use_reaper : bool
        Use REAPER for F0 estimation (requires pyreaper)
    frame_shift_ms : float
        Frame shift in milliseconds for output
    return_features : bool
        Include all 36 acoustic features in output
    return_probabilities : bool
        Include CD method probabilities (not just binary decisions)
        
    Returns
    -------
    dict
        Dictionary with NumPy arrays:
        - 'time': Time stamps in seconds
        - 'am_decisions': AM method binary decisions (if use_am=True)
        - 'cd_decisions': CD method binary decisions (if use_cd=True)
        - 'union_decisions': Union of AM and CD
        - 'antimode': F0 antimode value in Hz (if use_am=True)
        - 'F0': F0 contour in Hz (if use_am=True)
        - 'CD_prob': Classification probabilities (if return_probabilities=True)
        - 'features': Feature matrix (n_frames, 36) (if return_features=True)
        - 'cd_threshold': Threshold used
        - 'sample_rate': Sampling rate
        
    Examples
    --------
    For R usage with av package:
    ```r
    library(av)
    library(reticulate)
    
    # Load audio with av
    audio_data <- read_audio_bin("file.wav", channels=1)
    audio_info <- av_media_info("file.wav")
    
    # Call Python detector with extended output
    creak_mod <- import("union_creak_detection_method")
    result <- creak_mod$detect_creak_union_extended(
        audio = audio_data,
        sample_rate = audio_info$audio$sample_rate,
        return_features = TRUE,
        return_probabilities = TRUE
    )
    
    # Access all tracks
    am_creak <- result$am_decisions
    cd_creak <- result$cd_decisions
    union_creak <- result$union_decisions
    cd_prob <- result$CD_prob
    f0_track <- result$F0
    features <- result$features  # 36 columns
    ```
    """
    from .detector import CreakDetectorExtended
    
    # Ensure audio is 1D float array
    audio = np.asarray(audio, dtype=np.float64).flatten()
    
    # Resample if needed
    if sample_rate != 16000:
        from scipy.signal import resample
        n_samples = int(len(audio) * 16000 / sample_rate)
        audio = resample(audio, n_samples)
        sample_rate = 16000
    
    # Initialize detector
    detector = CreakDetectorExtended(
        use_trained_ann=True,
        use_reaper=use_reaper,
        cd_threshold=cd_threshold
    )
    
    # Process with extended output
    results = detector.process_extended(
        audio,
        sample_rate=sample_rate,
        use_am=use_am,
        use_cd=use_cd,
        frame_shift_ms=frame_shift_ms,
        return_features=return_features,
        return_probabilities=return_probabilities
    )
    
    # Add threshold to output
    results['cd_threshold'] = cd_threshold
    
    return results


def detect_creak_union(
    audio: np.ndarray,
    sample_rate: int = 16000,
    use_am: bool = True,
    use_cd: bool = True,
    cd_threshold: float = 0.3,
    use_reaper: bool = True,
    frame_shift_ms: float = 10.0
) -> Dict[str, np.ndarray]:
    """
    Detect creaky voice using the Union Method.
    
    Optimized for R integration:
    - Accepts NumPy array directly from av package
    - Returns NumPy arrays (automatic conversion to R vectors)
    - All processing in-memory (no file I/O)
    
    Parameters
    ----------
    audio : np.ndarray
        Audio signal as 1D NumPy array from av::read_audio_bin()
    sample_rate : int
        Sampling rate in Hz (default 16000)
    use_am : bool
        Use AM (antimode) method
    use_cd : bool
        Use CD (creak detector) method with trained ANN
    cd_threshold : float
        Classification threshold for CD method
    use_reaper : bool
        Use REAPER for F0 estimation (requires pyreaper)
    frame_shift_ms : float
        Frame shift in milliseconds for output
        
    Returns
    -------
    dict
        Dictionary with NumPy arrays:
        - 'time': Time stamps in seconds
        - 'am_decisions': AM method binary decisions (if use_am=True)
        - 'cd_decisions': CD method binary decisions (if use_cd=True)
        - 'union_decisions': Union of AM and CD
        - 'antimode': F0 antimode value in Hz (if use_am=True)
        - 'sample_rate': Sampling rate
        
    Examples
    --------
    For R usage with av package:
    ```r
    library(av)
    library(reticulate)
    
    # Load audio with av
    audio_info <- av_media_info("file.wav")
    audio_data <- read_audio_bin("file.wav", channels=1)
    
    # Call Python detector
    creak_mod <- import("union-creak-detection-method")
    result <- creak_mod$detect_creak_union(
        audio = audio_data,
        sample_rate = audio_info$audio$sample_rate
    )
    
    # Access results (automatic conversion to R)
    time <- result$time
    creak <- result$union_decisions
    ```
    """
    from .detector import CreakDetector
    
    # Ensure audio is 1D float array
    audio = np.asarray(audio, dtype=np.float64).flatten()
    
    # Resample if needed
    if sample_rate != 16000:
        from scipy.signal import resample
        n_samples = int(len(audio) * 16000 / sample_rate)
        audio = resample(audio, n_samples)
        sample_rate = 16000
    
    # Initialize detector
    detector = CreakDetector(
        use_trained_ann=True,
        use_reaper=use_reaper,
        cd_threshold=cd_threshold
    )
    
    # Process
    results = detector.process(
        audio,
        sample_rate=sample_rate,
        use_am=use_am,
        use_cd=use_cd,
        frame_shift_ms=frame_shift_ms
    )
    
    # Return NumPy arrays (reticulate handles conversion to R)
    return results


def extract_creak_features(
    audio: np.ndarray,
    sample_rate: int = 16000,
    frame_shift_ms: float = 10.0
) -> Dict[str, np.ndarray]:
    """
    Extract all 36 creak detection features.
    
    Returns features as NumPy arrays for efficient transfer to R.
    Can be used with superassp track system.
    
    Parameters
    ----------
    audio : np.ndarray
        Audio signal from av::read_audio_bin()
    sample_rate : int
        Sampling rate in Hz
    frame_shift_ms : float
        Frame shift in milliseconds
        
    Returns
    -------
    dict
        Dictionary with:
        - 'time': Time stamps
        - 'features': Feature matrix (n_frames, 36)
        - 'feature_names': List of feature names
        
    Examples
    --------
    For R track-based analysis:
    ```r
    # Extract features
    feats <- creak_mod$extract_creak_features(audio_data, sample_rate)
    
    # Convert to track format
    feature_matrix <- feats$features
    time_vec <- feats$time
    ```
    """
    from .features import get_all_features
    
    audio = np.asarray(audio, dtype=np.float64).flatten()
    
    # Resample if needed
    if sample_rate != 16000:
        from scipy.signal import resample
        n_samples = int(len(audio) * 16000 / sample_rate)
        audio = resample(audio, n_samples)
        sample_rate = 16000
    
    # Extract features
    features, time = get_all_features(
        audio,
        sample_rate=sample_rate,
        frame_shift_ms=frame_shift_ms
    )
    
    feature_names = [
        'H2_H1', 'peak_prominence', 'ZCR', 'IFP', 'IPS', 
        'PwP_fall', 'PwP_rise', 'F0', 'F0_mean', 'energy', 
        'power_std', 'creak_F0'
    ]
    
    # Add delta and delta-delta suffixes
    all_names = (feature_names + 
                [f'{n}_delta' for n in feature_names] +
                [f'{n}_deltadelta' for n in feature_names])
    
    return {
        'time': time,
        'features': features,
        'feature_names': all_names
    }

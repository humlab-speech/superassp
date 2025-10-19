"""
Audio I/O and preprocessing utilities
"""

import numpy as np
import soundfile as sf
from typing import Tuple, Optional
from scipy import signal


def load_audio(filepath: str, target_sr: Optional[int] = None) -> Tuple[np.ndarray, int]:
    """
    Load audio file
    
    Args:
        filepath: Path to audio file
        target_sr: Target sampling rate (None = keep original)
    
    Returns:
        audio: Audio signal (mono)
        sr: Sampling rate
    """
    audio, sr = sf.read(filepath)
    
    # Convert to mono if stereo
    if len(audio.shape) > 1:
        audio = audio[:, 0]
    
    # Resample if needed
    if target_sr is not None and sr != target_sr:
        audio = resample_audio(audio, sr, target_sr)
        sr = target_sr
    
    return audio, sr


def resample_audio(audio: np.ndarray, orig_sr: int, target_sr: int) -> np.ndarray:
    """
    Resample audio to target sampling rate
    
    Args:
        audio: Input audio
        orig_sr: Original sampling rate
        target_sr: Target sampling rate
    
    Returns:
        Resampled audio
    """
    if orig_sr == target_sr:
        return audio
    
    # Calculate number of samples in resampled signal
    n_samples = int(len(audio) * target_sr / orig_sr)
    
    # Use scipy's resample for high quality
    return signal.resample(audio, n_samples)


def preemphasize(audio: np.ndarray, coef: float = 0.96) -> np.ndarray:
    """
    Apply pre-emphasis filter
    
    Args:
        audio: Input audio
        coef: Pre-emphasis coefficient
    
    Returns:
        Pre-emphasized audio
    """
    return np.append(audio[0], audio[1:] - coef * audio[:-1])


def frame_audio(audio: np.ndarray, sr: int, frame_shift_ms: float = 1.0,
                window_size_ms: float = 25.0) -> Tuple[np.ndarray, np.ndarray]:
    """
    Convert audio to frames
    
    Args:
        audio: Input audio
        sr: Sampling rate
        frame_shift_ms: Frame shift in milliseconds
        window_size_ms: Window size in milliseconds
    
    Returns:
        frames: 2D array of frames (n_frames, frame_length)
        times: Frame center times in seconds
    """
    frame_shift = int(sr * frame_shift_ms / 1000)
    frame_length = int(sr * window_size_ms / 1000)
    
    # Calculate number of frames
    n_frames = (len(audio) - frame_length) // frame_shift + 1
    
    # Pre-allocate
    frames = np.zeros((n_frames, frame_length))
    times = np.zeros(n_frames)
    
    # Extract frames
    for i in range(n_frames):
        start = i * frame_shift
        end = start + frame_length
        if end <= len(audio):
            frames[i] = audio[start:end]
            times[i] = (start + frame_length / 2) / sr
    
    return frames, times


def get_segment_indices(audio_len: int, sr: int, f0: float, center_sample: int,
                       n_periods: int) -> Tuple[int, int]:
    """
    Get start and end indices for a segment centered at center_sample
    
    Args:
        audio_len: Length of audio signal
        sr: Sampling rate
        f0: Fundamental frequency (Hz)
        center_sample: Center sample index
        n_periods: Number of pitch periods to extract
    
    Returns:
        start: Start index
        end: End index
    """
    if f0 <= 0 or np.isnan(f0):
        return -1, -1
    
    samples_per_period = sr / f0
    half_window = int(n_periods / 2 * samples_per_period)
    
    start = int(center_sample - half_window)
    end = int(center_sample + half_window)
    
    # Check bounds
    if start < 0 or end >= audio_len:
        return -1, -1
    
    return start, end


def normalize_audio(audio: np.ndarray) -> np.ndarray:
    """
    Normalize audio to [-1, 1] range
    
    Args:
        audio: Input audio
    
    Returns:
        Normalized audio
    """
    max_val = np.abs(audio).max()
    if max_val > 0:
        return audio / max_val
    return audio

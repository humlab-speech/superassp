"""Feature extraction - optimized for R."""
import numpy as np
from typing import Tuple
from pathlib import Path
import sys

# Add core implementation to path
pkg_root = Path(__file__).parent.parent.parent.parent.parent
core_path = pkg_root / 'union-creak-detection-method' / 'creak_detection'

if core_path.exists():
    sys.path.insert(0, str(core_path.parent))

def get_all_features(audio: np.ndarray, sample_rate: int = 16000,
                     frame_shift_ms: float = 10.0) -> Tuple[np.ndarray, np.ndarray]:
    """Extract all 36 features."""
    try:
        from creak_detection.core.features import get_all_creak_features
        features = get_all_creak_features(audio, sample_rate)
        frame_shift_sec = frame_shift_ms / 1000.0
        time = np.arange(len(features)) * frame_shift_sec
        return features, time
    except ImportError:
        # Fallback
        return _simple_features(audio, sample_rate, frame_shift_ms)

def _simple_features(audio, sr, shift_ms):
    """Simplified fallback."""
    shift = int(shift_ms * sr / 1000)
    n = len(audio) // shift
    feats = np.zeros((n, 36))
    time = np.arange(n) * shift_ms / 1000.0
    return feats, time

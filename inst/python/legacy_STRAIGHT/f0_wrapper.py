"""
Wrapper for STRAIGHT F0 extraction that's more reticulate-friendly.
Avoids in-place modifications and complex numpy operations.
"""
import numpy as np
from typing import Dict, Tuple
from .f0_extraction import MulticueF0v14, F0Parameters


def extract_f0_safe(x_list, fs, f0floor=71, f0ceil=800):
    """
    Reticulate-safe wrapper for F0 extraction.
    
    Args:
        x_list: Audio as Python list or numpy array
        fs: Sample rate
        f0floor: Minimum F0
        f0ceil: Maximum F0
        
    Returns:
        dict with keys: 'f0', 'vuv', 'aux'
    """
    # Convert to numpy array (copy to avoid reticulate issues)
    x = np.array(x_list, dtype=np.float32, copy=True)
    
    # Call main function
    f0_raw, vuv, auxouts, prm = MulticueF0v14(x, float(fs), float(f0floor), float(f0ceil))
    
    # Return simple dict (reticulate handles this better than tuple)
    return {
        'f0': f0_raw.tolist(),  # Convert to list for safe transfer
        'vuv': vuv.tolist(),
        'aux': {
            'RELofcandidatesByMix': auxouts.get('RELofcandidatesByMix', np.array([])).tolist() if 'RELofcandidatesByMix' in auxouts else [],
            'ACofcandidatesByAC': auxouts.get('ACofcandidatesByAC', np.array([])).tolist() if 'ACofcandidatesByAC' in auxouts else []
        }
    }

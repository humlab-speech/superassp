"""
Glottal Quotient (GQ) Features

Ported from glottis_quotient function in voice_analysis_redux.m (lines 409-440)
Requires DYPSA for glottal closure/opening instant detection
"""

import numpy as np
import warnings


def compute_glottal_quotient(audio, fs, f0_min=50, f0_max=500):
    """
    Compute Glottal Quotient features
    
    Requires DYPSA algorithm for GCI/GOI detection
    
    Parameters:
    -----------
    audio : ndarray
        Audio signal
    fs : int
        Sampling frequency
    f0_min : float
        Minimum F0
    f0_max : float
        Maximum F0
        
    Returns:
    --------
    measures : dict
        3 GQ-related measures
    """
    try:
        from ..utils.dypsa import dypsa
        
        audio = np.asarray(audio).ravel()
        
        # Run DYPSA to get GCI and GOI
        gci, goi = dypsa(audio, fs, f0_min, f0_max)
        
        if len(gci) < 2 or len(goi) < 1:
            warnings.warn("Insufficient glottal instants detected")
            return _get_nan_gq()
        
        # Compute open and closed cycles (lines 412-413)
        cycle_open = []
        cycle_closed = []
        
        for i in range(len(gci) - 1):
            # Find GOI between consecutive GCIs
            goi_in_cycle = goi[(goi > gci[i]) & (goi < gci[i + 1])]
            
            if len(goi_in_cycle) > 0:
                goi_time = goi_in_cycle[0]
                
                # Open phase: GCI to GOI
                closed_dur = goi_time - gci[i]
                
                # Closed phase: GOI to next GCI  
                open_dur = gci[i + 1] - goi_time
                
                cycle_open.append(open_dur)
                cycle_closed.append(closed_dur)
        
        if len(cycle_open) == 0:
            return _get_nan_gq()
        
        cycle_open = np.array(cycle_open)
        cycle_closed = np.array(cycle_closed)
        
        # Remove erroneous cycles (lines 416-428)
        low_lim = fs / f0_max
        up_lim = fs / f0_min
        
        valid_open = []
        valid_closed = []
        
        for i in range(len(cycle_open)):
            if low_lim <= cycle_open[i] <= up_lim:
                valid_open.append(cycle_open[i])
            
            if low_lim <= cycle_closed[i] <= up_lim:
                valid_closed.append(cycle_closed[i])
        
        if len(valid_open) == 0 or len(valid_closed) == 0:
            return _get_nan_gq()
        
        valid_open = np.array(valid_open)
        valid_closed = np.array(valid_closed)
        
        # Statistics (lines 430-435)
        prc_open = np.percentile(valid_open, [5, 95])
        cycle_open_range = prc_open[1] - prc_open[0]
        
        prc_closed = np.percentile(valid_closed, [5, 95])
        cycle_closed_range = prc_closed[1] - prc_closed[0]
        
        measures = {}
        
        # Line 436: Glottal Quotient (ratio of ranges)
        measures['GQ'] = (
            cycle_open_range / (cycle_open_range + cycle_closed_range + 1e-10)
        )
        
        # Line 437: Std of open cycle
        measures['GQ_open_std'] = np.std(valid_open)
        
        # Line 438: Std of closed cycle
        measures['GQ_closed_std'] = np.std(valid_closed)
        
        return measures
        
    except ImportError:
        warnings.warn("DYPSA module not available")
        return _get_nan_gq()
    except Exception as e:
        warnings.warn(f"Glottal Quotient computation failed: {e}")
        return _get_nan_gq()


def _get_nan_gq():
    """Return NaN values for all GQ measures"""
    return {
        'GQ': np.nan,
        'GQ_open_std': np.nan,
        'GQ_closed_std': np.nan,
    }

"""
Jitter and Shimmer Measures

Ported from jitter_shimmer function in voice_analysis_redux.m (lines 262-310)
Computes perturbation measures for F0 (jitter) or amplitude (shimmer)
"""

import numpy as np
from ..utils.tkeo import compute_tkeo_vectorized
from ..utils.perturbation import (compute_perturbation_quotient, 
                                   compute_ar_perturbation_quotient,
                                   compute_nmsp)


def compute_jitter_shimmer_features(time_series, measure_type='jitter', use_thesis_mode=False):
    """
    Compute jitter or shimmer measures
    
    Parameters:
    -----------
    time_series : ndarray
        F0 contour (for jitter) or amplitude contour (for shimmer)
    measure_type : str
        'jitter' or 'shimmer' for naming convention
    use_thesis_mode : bool
        If True, include thesis-specific measures (AR jitter, NMSP, thesis shimmer dB)
        If False (default), use MATLAB implementation only
        
    Returns:
    --------
    measures : dict
        Dictionary of perturbation measures (22 standard + 3 optional thesis measures)
    """
    time_series = np.asarray(time_series).ravel()
    
    # Remove zeros and invalid values
    time_series = time_series[time_series > 0]
    
    if len(time_series) < 2:
        return _get_nan_measures(measure_type, use_thesis_mode)
    
    prefix = measure_type
    measures = {}
    mean_val = np.mean(time_series)
    
    if mean_val == 0:
        return _get_nan_measures(measure_type, use_thesis_mode)
    
    # 1. RAP - Mean absolute difference (Relative Average Perturbation)
    measures[f'{prefix}_RAP'] = np.mean(np.abs(np.diff(time_series)))
    
    # 2. RAP% - Expressed as percentage
    measures[f'{prefix}_RAP_percent'] = 100 * measures[f'{prefix}_RAP'] / mean_val
    
    # 3-5. Perturbation Quotient K=3
    if len(time_series) >= 3:
        pq3 = compute_perturbation_quotient(time_series, K=3)
        measures[f'{prefix}_PQ3_Schoentgen'] = pq3['classical_Schoentgen']
        measures[f'{prefix}_PQ3_Baken'] = pq3['classical_Baken']
        measures[f'{prefix}_PQ3_generalized'] = pq3['generalized_Schoentgen']
    else:
        measures[f'{prefix}_PQ3_Schoentgen'] = np.nan
        measures[f'{prefix}_PQ3_Baken'] = np.nan
        measures[f'{prefix}_PQ3_generalized'] = np.nan
    
    # 6-8. Perturbation Quotient K=5
    if len(time_series) >= 5:
        pq5 = compute_perturbation_quotient(time_series, K=5)
        measures[f'{prefix}_PQ5_Schoentgen'] = pq5['classical_Schoentgen']
        measures[f'{prefix}_PQ5_Baken'] = pq5['classical_Baken']
        measures[f'{prefix}_PQ5_generalized'] = pq5['generalized_Schoentgen']
    else:
        measures[f'{prefix}_PQ5_Schoentgen'] = np.nan
        measures[f'{prefix}_PQ5_Baken'] = np.nan
        measures[f'{prefix}_PQ5_generalized'] = np.nan
    
    # 9-11. Perturbation Quotient K=11
    if len(time_series) >= 11:
        pq11 = compute_perturbation_quotient(time_series, K=11)
        measures[f'{prefix}_PQ11_Schoentgen'] = pq11['classical_Schoentgen']
        measures[f'{prefix}_PQ11_Baken'] = pq11['classical_Baken']
        measures[f'{prefix}_PQ11_generalized'] = pq11['generalized_Schoentgen']
    else:
        measures[f'{prefix}_PQ11_Schoentgen'] = np.nan
        measures[f'{prefix}_PQ11_Baken'] = np.nan
        measures[f'{prefix}_PQ11_generalized'] = np.nan
    
    # 12. Zeroth-order perturbation
    measures[f'{prefix}_zeroth_order'] = np.mean(np.abs(time_series - mean_val))
    
    # 13. Shimmer (dB) - only computed for amplitude measures
    if measure_type == 'shimmer':
        if use_thesis_mode:
            # Thesis-compliant version (Equation 3.40 equivalent for amplitude)
            # Standard shimmer dB: mean of absolute log ratios
            ratio = time_series[1:] / (time_series[:-1] + 1e-10)
            measures[f'{prefix}_dB'] = np.mean(np.abs(20 * np.log10(ratio + 1e-10)))
        else:
            # MATLAB implementation (matches original code)
            ratio = time_series[:-1] / (time_series[1:] + 1e-10)
            measures[f'{prefix}_dB'] = np.mean(20 * np.abs(np.log10(np.abs(ratio) + 1e-10)))
    else:
        # For jitter, we skip this measure or set to NaN
        measures[f'{prefix}_dB'] = np.nan
    
    # 14. Coefficient of Variation
    measures[f'{prefix}_CV'] = np.mean(np.diff(time_series)**2) / (mean_val**2)
    
    # 15-16. TKEO mean and std
    tkeo_vals = compute_tkeo_vectorized(time_series)
    measures[f'{prefix}_TKEO_mean'] = np.mean(np.abs(tkeo_vals))
    measures[f'{prefix}_TKEO_std'] = np.std(tkeo_vals)
    
    # 17-20. TKEO percentiles
    percentiles = np.percentile(tkeo_vals, [5, 25, 50, 75, 95])
    measures[f'{prefix}_TKEO_p5'] = percentiles[0]
    measures[f'{prefix}_TKEO_p25'] = percentiles[1]
    measures[f'{prefix}_TKEO_p50'] = percentiles[2]
    measures[f'{prefix}_TKEO_p75'] = percentiles[3]
    # Note: p95 is computed but not directly stored, used for IQR below
    
    # 21. Amplitude Modulation
    measures[f'{prefix}_AM'] = (np.max(time_series) - np.min(time_series)) / \
                                (np.max(time_series) + np.min(time_series))
    
    # 22. TKEO IQR (p75 - p5)
    measures[f'{prefix}_TKEO_IQR'] = percentiles[3] - percentiles[0]
    
    # Additional thesis-specific measures (only if requested)
    if use_thesis_mode:
        # 23. AR-based Perturbation Quotient (Equation 3.39)
        measures[f'{prefix}_PQ_AR'] = compute_ar_perturbation_quotient(time_series, ar_order=10)
        
        # 24. Normalized Mean Squared Perturbation (Equation 3.41)
        measures[f'{prefix}_NMSP'] = compute_nmsp(time_series)
        
        # 25. F0 Range (robust, using percentiles)
        # Thesis mentions F0_range = F0_95th - F0_5th percentile
        if measure_type == 'jitter':  # Only for F0 measures
            f0_percentiles = np.percentile(time_series, [5, 95])
            measures[f'{prefix}_F0_range'] = f0_percentiles[1] - f0_percentiles[0]
        else:
            measures[f'{prefix}_amp_range'] = np.nan  # Not applicable for shimmer
    
    return measures


def _get_nan_measures(measure_type, use_thesis_mode=False):
    """Return dictionary with NaN values for all measures"""
    prefix = measure_type
    base_measures = {
        f'{prefix}_RAP': np.nan,
        f'{prefix}_RAP_percent': np.nan,
        f'{prefix}_PQ3_Schoentgen': np.nan,
        f'{prefix}_PQ3_Baken': np.nan,
        f'{prefix}_PQ3_generalized': np.nan,
        f'{prefix}_PQ5_Schoentgen': np.nan,
        f'{prefix}_PQ5_Baken': np.nan,
        f'{prefix}_PQ5_generalized': np.nan,
        f'{prefix}_PQ11_Schoentgen': np.nan,
        f'{prefix}_PQ11_Baken': np.nan,
        f'{prefix}_PQ11_generalized': np.nan,
        f'{prefix}_zeroth_order': np.nan,
        f'{prefix}_dB': np.nan,
        f'{prefix}_CV': np.nan,
        f'{prefix}_TKEO_mean': np.nan,
        f'{prefix}_TKEO_std': np.nan,
        f'{prefix}_TKEO_p5': np.nan,
        f'{prefix}_TKEO_p25': np.nan,
        f'{prefix}_TKEO_p50': np.nan,
        f'{prefix}_TKEO_p75': np.nan,
        f'{prefix}_AM': np.nan,
        f'{prefix}_TKEO_IQR': np.nan,
    }
    
    if use_thesis_mode:
        base_measures.update({
            f'{prefix}_PQ_AR': np.nan,
            f'{prefix}_NMSP': np.nan,
            f'{prefix}_F0_range' if measure_type == 'jitter' else f'{prefix}_amp_range': np.nan,
        })
    
    return base_measures

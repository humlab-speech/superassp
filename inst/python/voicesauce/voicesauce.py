"""
Main VoiceSauce API
"""

import numpy as np
from typing import Optional, Union
from pathlib import Path

from .core.config import VoiceSauceConfig
from .core.audio import load_audio
from .f0 import estimate_f0
from .formants import estimate_formants
from .measures import (
    get_harmonics, get_formant_amplitudes,
    get_cpp, get_hnr, get_energy,
    get_2k, get_5k, get_2k5k, get_h42k,
    apply_corrections,
    compute_harmonic_differences,
    compute_harmonic_formant_differences
)


class VoiceSauceResults:
    """Container for VoiceSauce analysis results"""
    
    def __init__(self, config: VoiceSauceConfig):
        self.config = config
        self.results = {}
        self.times = None
        self.fs = None
    
    def __getitem__(self, key: str):
        """Access results by key"""
        return self.results.get(key)
    
    def __setitem__(self, key: str, value):
        """Set results by key"""
        self.results[key] = value
    
    def update(self, other: dict):
        """Update results with dictionary"""
        self.results.update(other)
    
    def keys(self):
        """Return available measure names"""
        return self.results.keys()
    
    def to_dict(self) -> dict:
        """Convert to dictionary"""
        return self.results.copy()
    
    def to_dataframe(self):
        """Convert to pandas DataFrame"""
        try:
            import pandas as pd
            
            # Find the common length (use F0 as reference)
            ref_length = len(self.results.get('F0', []))
            
            # Build dictionary with only arrays of correct length
            df_dict = {}
            if self.times is not None and len(self.times) == ref_length:
                df_dict['Time'] = self.times
            
            for key, values in self.results.items():
                if isinstance(values, np.ndarray) and len(values) == ref_length:
                    df_dict[key] = values
            
            return pd.DataFrame(df_dict)
        except ImportError:
            raise ImportError("pandas is required for to_dataframe()")
    
    def summary(self) -> str:
        """Generate summary statistics"""
        lines = ["VoiceSauce Analysis Results", "=" * 40]
        
        for key, values in self.results.items():
            if isinstance(values, np.ndarray):
                valid = values[~np.isnan(values)]
                if len(valid) > 0:
                    lines.append(f"{key}:")
                    lines.append(f"  Mean: {np.mean(valid):.2f}")
                    lines.append(f"  Std:  {np.std(valid):.2f}")
                    lines.append(f"  Min:  {np.min(valid):.2f}")
                    lines.append(f"  Max:  {np.max(valid):.2f}")
                    lines.append(f"  Valid: {len(valid)}/{len(values)}")
        
        return "\n".join(lines)


def analyze(audio_path: Union[str, Path],
           config: Optional[VoiceSauceConfig] = None,
           textgrid_path: Optional[Union[str, Path]] = None) -> VoiceSauceResults:
    """
    Perform complete VoiceSauce analysis
    
    Args:
        audio_path: Path to audio file
        config: VoiceSauceConfig object (None = use defaults)
        textgrid_path: Optional TextGrid file for segmentation
    
    Returns:
        VoiceSauceResults object with all measures
    """
    # Use default config if none provided
    if config is None:
        config = VoiceSauceConfig()
    
    # Load audio
    audio, fs = load_audio(str(audio_path))
    
    # Create results container
    results = VoiceSauceResults(config)
    results.fs = fs
    
    # Step 1: F0 estimation
    print(f"Estimating F0 using {config.f0_method}...")
    f0_values, f0_times = estimate_f0(
        audio, fs,
        method=config.f0_method,
        frame_shift=config.frame_shift,
        f0_min=config.f0_min,
        f0_max=config.f0_max
    )
    
    # Convert times to match frame shift
    times_ms = np.arange(len(f0_values)) * config.frame_shift
    results.times = times_ms / 1000.0  # Convert to seconds
    results['F0'] = f0_values
    
    # Step 2: Formant estimation
    print(f"Estimating formants using {config.formant_method}...")
    formants, bandwidths, fmt_times = estimate_formants(
        audio, fs,
        method=config.formant_method,
        frame_shift=config.frame_shift,
        window_length=config.window_size / 1000.0,
        n_formants=config.n_formants,
        max_formant=config.max_formant
    )
    
    # Align formants to F0 length
    if len(formants) != len(f0_values):
        # Resize formants to match F0 length
        n_frames = len(f0_values)
        formants_aligned = np.full((n_frames, formants.shape[1]), np.nan)
        bandwidths_aligned = np.full((n_frames, bandwidths.shape[1]), np.nan)
        
        # Copy what we can
        copy_len = min(len(formants), n_frames)
        formants_aligned[:copy_len, :] = formants[:copy_len, :]
        bandwidths_aligned[:copy_len, :] = bandwidths[:copy_len, :]
        
        formants = formants_aligned
        bandwidths = bandwidths_aligned
    
    # Extract individual formants
    for i in range(min(5, formants.shape[1])):
        results[f'F{i+1}'] = formants[:, i]
        results[f'B{i+1}'] = bandwidths[:, i]
    
    # Step 3: Harmonic amplitudes (H1, H2, H4)
    print("Calculating harmonic amplitudes...")
    h1, h2, h4 = get_harmonics(
        audio, fs, f0_values,
        frame_shift=config.frame_shift,
        n_periods=config.n_periods
    )
    
    results['H1'] = h1
    results['H2'] = h2
    results['H4'] = h4
    
    # Harmonic differences
    harm_diffs = compute_harmonic_differences(h1, h2, h4)
    results.update(harm_diffs)
    
    # Step 4: Formant amplitudes (A1, A2, A3)
    print("Calculating formant amplitudes...")
    a1, a2, a3 = get_formant_amplitudes(
        audio, fs, f0_values,
        results['F1'], results['F2'], results['F3'],
        frame_shift=config.frame_shift,
        n_periods=config.n_periods
    )
    
    results['A1'] = a1
    results['A2'] = a2
    results['A3'] = a3
    
    # Harmonic-formant differences
    hf_diffs = compute_harmonic_formant_differences(h1, a1, a2, a3)
    results.update(hf_diffs)
    
    # Step 5: CPP
    print("Calculating CPP...")
    cpp = get_cpp(
        audio, fs, f0_values,
        frame_shift=config.frame_shift,
        n_periods=config.n_periods_ec
    )
    results['CPP'] = cpp
    
    # Step 6: HNR
    print("Calculating HNR...")
    hnr_dict = get_hnr(
        audio, fs, f0_values,
        frame_shift=config.frame_shift,
        n_periods=config.n_periods_ec
    )
    results.update(hnr_dict)
    
    # Step 7: Energy
    print("Calculating Energy...")
    energy = get_energy(
        audio, fs, f0_values,
        frame_shift=config.frame_shift,
        n_periods=config.n_periods_ec
    )
    results['Energy'] = energy
    
    # Step 8: Spectral measures
    print("Calculating spectral measures...")
    two_k = get_2k(audio, fs, f0_values, frame_shift=config.frame_shift)
    five_k = get_5k(audio, fs, f0_values, frame_shift=config.frame_shift)
    
    results['2K'] = two_k
    results['5K'] = five_k
    results['2K5K'] = get_2k5k(two_k, five_k)
    results['H42K'] = get_h42k(h4, two_k)
    
    # Step 9: Apply Iseli-Alwan corrections
    print("Applying formant corrections...")
    corrected = apply_corrections(
        h1, h2, h4, a1, a2, a3,
        f0_values,
        results['F1'], results['F2'], results['F3'],
        results['B1'], results['B2'], results['B3'],
        fs
    )
    results.update(corrected)
    
    # Also compute corrected H42K
    results['H42Kc'] = corrected['H4c'] - two_k
    
    print("Analysis complete!")
    
    return results

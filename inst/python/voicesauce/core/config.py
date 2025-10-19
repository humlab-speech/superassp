"""
VoiceSauce Configuration Module
"""

from dataclasses import dataclass, field
from typing import List, Optional


@dataclass
class VoiceSauceConfig:
    """Configuration for VoiceSauce analysis"""
    
    # Frame processing
    frame_shift: float = 1.0  # ms
    window_size: float = 25.0  # ms
    
    # F0 estimation
    f0_method: str = 'reaper'  # 'reaper', 'praat', 'shr', 'world'
    f0_min: float = 40.0  # Hz
    f0_max: float = 500.0  # Hz
    
    # Formant estimation
    formant_method: str = 'praat'  # 'praat', 'sptk'
    n_formants: int = 5
    max_formant: float = 5500.0  # Hz
    
    # Voice quality measures
    n_periods: int = 3  # For H1, H2, H4
    n_periods_ec: int = 5  # For CPP, HNR, Energy
    
    # SHR parameters
    shr_threshold: float = 0.4
    shr_ceiling: float = 1250.0  # Hz
    
    # Praat parameters (when using Praat)
    praat_silence_threshold: float = 0.03
    praat_voicing_threshold: float = 0.45
    praat_octave_cost: float = 0.01
    praat_octave_jump_cost: float = 0.35
    praat_voiced_unvoiced_cost: float = 0.14
    
    # TextGrid processing
    textgrid_tier: int = 1
    textgrid_ignore_labels: List[str] = field(default_factory=lambda: ['', ' ', 'SIL'])
    textgrid_buffer: float = 25.0  # ms
    
    # Output options
    output_measures: Optional[List[str]] = None  # None = all measures
    
    # Processing options
    frame_precision: int = 1  # Frame alignment tolerance
    preemphasis: float = 0.96
    
    def __post_init__(self):
        """Validate configuration"""
        valid_f0_methods = ['reaper', 'praat', 'shr', 'world']
        if self.f0_method not in valid_f0_methods:
            raise ValueError(f"f0_method must be one of {valid_f0_methods}")
        
        valid_formant_methods = ['praat', 'sptk']
        if self.formant_method not in valid_formant_methods:
            raise ValueError(f"formant_method must be one of {valid_formant_methods}")
        
        if self.f0_min >= self.f0_max:
            raise ValueError("f0_min must be less than f0_max")
        
        if self.frame_shift <= 0:
            raise ValueError("frame_shift must be positive")
    
    def get_all_measures(self) -> List[str]:
        """Return list of all available measures"""
        return [
            'F0',
            'F1', 'F2', 'F3', 'F4', 'F5',
            'B1', 'B2', 'B3', 'B4', 'B5',
            'H1', 'H2', 'H4',
            'A1', 'A2', 'A3',
            'H1_H2', 'H2_H4',
            'H1_A1', 'H1_A2', 'H1_A3',
            'H1c', 'H2c', 'H4c',
            'A1c', 'A2c', 'A3c',
            'H1c_H2c', 'H2c_H4c',
            'H1c_A1c', 'H1c_A2c', 'H1c_A3c',
            'CPP',
            'HNR05', 'HNR15', 'HNR25', 'HNR35',
            'Energy',
            'SHR',
            'H42K', 'H42Kc',
            '2K', '5K', '2K5K'
        ]

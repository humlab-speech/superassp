"""
Core Voice Analysis Module

Main VoiceAnalyzer class and convenience functions
Ported from voice_analysis_redux.m
"""

import numpy as np
import soundfile as sf
from scipy import signal as scipy_signal
import warnings

from .f0_estimation import estimate_f0_praat, estimate_f0_swipe
from .features import (
    compute_jitter_shimmer_features,
    compute_hnr_nhr,
    compute_mfcc_features,
    compute_wavelet_features,
    compute_ppe,
    compute_dfa,
    compute_rpde,
    compute_gne,
    compute_emd_features,
    compute_glottal_quotient,
    compute_vfer,
)


class VoiceAnalyzer:
    """
    Comprehensive Voice Analysis
    
    Computes 132 dysphonia measures from sustained vowel recordings
    """
    
    def __init__(self, f0_min=50, f0_max=500, f0_algorithm='SWIPE', use_thesis_mode=False):
        """
        Initialize VoiceAnalyzer
        
        Parameters:
        -----------
        f0_min : float
            Minimum F0 (Hz), default 50
        f0_max : float
            Maximum F0 (Hz), default 500
        f0_algorithm : str
            F0 estimation algorithm: 'SWIPE', 'PRAAT', default 'SWIPE'
        use_thesis_mode : bool
            If True, use thesis-compliant implementations (Tsanas 2012)
            If False (default), use MATLAB-compatible implementations
            Affects: PPE (semitone vs log), shimmer dB, adds AR jitter & NMSP
        """
        self.f0_min = f0_min
        self.f0_max = f0_max
        self.f0_algorithm = f0_algorithm.upper()
        self.use_thesis_mode = use_thesis_mode
        
    def analyze(self, audio, fs):
        """
        Perform comprehensive voice analysis
        
        Parameters:
        -----------
        audio : ndarray
            Audio signal (mono)
        fs : int
            Sampling frequency (Hz)
            
        Returns:
        --------
        measures : dict
            Dictionary of all computed measures
        F0 : ndarray
            Fundamental frequency contour
        """
        # Preprocessing (lines 129-131)
        audio = np.asarray(audio).ravel()
        audio = audio - np.mean(audio)
        audio = audio / (np.max(np.abs(audio)) + 1e-10)
        
        # F0 estimation (lines 133-155)
        F0 = self._estimate_f0(audio, fs)
        
        # Amplitude contour (lines 157-168)
        A0 = self._extract_amplitude(audio, fs)
        
        # Initialize measures dictionary
        measures = {}
        
        print("Computing features...")
        
        # 1. Jitter measures (22 + 3 optional) - lines 176-178
        print("  - Jitter...")
        jitter = compute_jitter_shimmer_features(F0, 'jitter', use_thesis_mode=self.use_thesis_mode)
        measures.update(jitter)
        
        # 2. Shimmer measures (22 + 3 optional) - lines 182
        print("  - Shimmer...")
        shimmer = compute_jitter_shimmer_features(A0, 'shimmer', use_thesis_mode=self.use_thesis_mode)
        measures.update(shimmer)
        
        # 3. HNR/NHR (4) - lines 184-185
        print("  - HNR/NHR...")
        hnr_nhr = compute_hnr_nhr(audio, fs, self.f0_min, self.f0_max)
        measures.update(hnr_nhr)
        
        # 4. DFA (1) - lines 187-193
        print("  - DFA...")
        dfa_scaling = np.arange(50, 201, 20)
        measures['DFA'] = compute_dfa(audio, scales=dfa_scaling)
        
        # 5. RPDE (1) - lines 195-201
        print("  - RPDE...")
        audio_resampled = scipy_signal.resample(audio, int(25000 * len(audio) / fs))
        measures['RPDE'] = compute_rpde(audio_resampled, m=4, tau=50, epsilon=0.12, T_max=1000)
        
        # 6. PPE (1) - lines 203-211
        print("  - PPE...")
        measures['PPE'] = compute_ppe(F0, fs, f0_mean_healthy=120, use_thesis_mode=self.use_thesis_mode)
        
        # 7. GNE (6) - lines 220-222
        print("  - GNE...")
        gne = compute_gne(audio, fs)
        measures.update(gne)
        
        # 8. MFCCs (84) - lines 238-245
        print("  - MFCCs...")
        mfcc = compute_mfcc_features(audio, fs)
        measures.update(mfcc)
        
        # 9. Wavelet features (~50) - line 248
        print("  - Wavelet features...")
        wavelet = compute_wavelet_features(F0)
        measures.update(wavelet)
        
        # 10. Glottal Quotient (3) - lines 213-218
        print("  - Glottal Quotient...")
        gq = compute_glottal_quotient(audio, fs, self.f0_min, self.f0_max)
        measures.update(gq)
        
        # 11. VFER (7) - lines 223-228
        print("  - VFER...")
        vfer = compute_vfer(audio, fs)
        measures.update(vfer)
        
        # 12. EMD features (6) - lines 230-236
        print("  - EMD features...")
        emd = compute_emd_features(audio)
        measures.update(emd)
        
        print(f"Computed {len(measures)} measures")
        
        return measures, F0
    
    def _estimate_f0(self, audio, fs):
        """Estimate fundamental frequency"""
        if self.f0_algorithm == 'SWIPE':
            return estimate_f0_swipe(audio, fs, self.f0_min, self.f0_max)
        elif self.f0_algorithm == 'PRAAT':
            return estimate_f0_praat(audio, fs, self.f0_min, self.f0_max)
        else:
            warnings.warn(f"Unknown F0 algorithm '{self.f0_algorithm}', using SWIPE")
            return estimate_f0_swipe(audio, fs, self.f0_min, self.f0_max)
    
    def _extract_amplitude(self, audio, fs, frame_shift=0.01):
        """
        Extract amplitude contour
        
        Fallback method (lines 166-168) when DYPSA is not available
        """
        frame_len = int(frame_shift * fs)
        n_frames = len(audio) // frame_len
        
        if n_frames == 0:
            return np.array([np.max(np.abs(audio))])
        
        A0 = np.zeros(n_frames)
        for i in range(n_frames):
            start_idx = i * frame_len
            end_idx = min(start_idx + frame_len, len(audio))
            frame = audio[start_idx:end_idx]
            A0[i] = np.max(np.abs(frame))
        
        return A0


def analyze_voice(audio_path_or_array, fs=None, f0_min=50, f0_max=500, f0_algorithm='SWIPE'):
    """
    Convenience function for voice analysis
    
    Parameters:
    -----------
    audio_path_or_array : str or ndarray
        Path to WAV file or audio array
    fs : int
        Sampling frequency (required if audio_path_or_array is array)
    f0_min : float
        Minimum F0 (Hz)
    f0_max : float
        Maximum F0 (Hz)
    f0_algorithm : str
        'SWIPE' or 'PRAAT'
        
    Returns:
    --------
    measures : dict
        All computed measures
    F0 : ndarray
        F0 contour
    """
    # Load audio if path provided
    if isinstance(audio_path_or_array, str):
        audio, fs = sf.read(audio_path_or_array)
        
        # Convert stereo to mono
        if audio.ndim > 1:
            audio = np.mean(audio, axis=1)
    else:
        audio = audio_path_or_array
        if fs is None:
            raise ValueError("Sampling frequency (fs) required when providing audio array")
    
    # Create analyzer and run
    analyzer = VoiceAnalyzer(f0_min=f0_min, f0_max=f0_max, f0_algorithm=f0_algorithm)
    measures, F0 = analyzer.analyze(audio, fs)
    
    return measures, F0


def analyze_voice_file(filepath, f0_min=50, f0_max=500, f0_algorithm='SWIPE'):
    """
    Convenience function to analyze a WAV file
    
    Parameters:
    -----------
    filepath : str
        Path to WAV file
    f0_min : float
        Minimum F0 (Hz)
    f0_max : float
        Maximum F0 (Hz)
    f0_algorithm : str
        'SWIPE' or 'PRAAT'
        
    Returns:
    --------
    measures : dict
        All computed measures
    F0 : ndarray
        F0 contour
    """
    return analyze_voice(filepath, f0_min=f0_min, f0_max=f0_max, f0_algorithm=f0_algorithm)

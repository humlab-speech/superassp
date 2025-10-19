"""
Core Voice Analysis Module - PARALLEL VERSION

Parallelized implementation using ThreadPoolExecutor for feature-level parallelization.
Expected speedup: 3-5x on multi-core systems.

Main VoiceAnalyzer class with parallel feature computation
Ported from voice_analysis_redux.m with parallelization enhancements
"""

import numpy as np
import soundfile as sf
from scipy import signal as scipy_signal
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from functools import partial

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


class VoiceAnalyzerParallel:
    """
    Comprehensive Voice Analysis - PARALLEL VERSION
    
    Computes 132+ dysphonia measures from sustained vowel recordings
    Uses ThreadPoolExecutor for parallel feature computation
    """
    
    def __init__(self, f0_min=50, f0_max=500, f0_algorithm='SWIPE', max_workers=None, 
                 enable_within_feature_parallel=True, use_rpde_kdtree=True):
        """
        Initialize VoiceAnalyzerParallel
        
        Parameters:
        -----------
        f0_min : float
            Minimum F0 (Hz), default 50
        f0_max : float
            Maximum F0 (Hz), default 500
        f0_algorithm : str
            F0 estimation algorithm: 'SWIPE', 'PRAAT', default 'SWIPE'
        max_workers : int or None
            Number of parallel workers (None = auto-detect)
        enable_within_feature_parallel : bool
            Enable parallelization within HNR/GNE computations (default True)
        use_rpde_kdtree : bool
            Use KD-tree for RPDE (2-5x faster, default True)
        """
        self.f0_min = f0_min
        self.f0_max = f0_max
        self.f0_algorithm = f0_algorithm.upper()
        self.max_workers = max_workers
        self.enable_within_feature_parallel = enable_within_feature_parallel
        self.use_rpde_kdtree = use_rpde_kdtree
        
    def analyze(self, audio, fs, verbose=True):
        """
        Perform comprehensive voice analysis with parallelization
        
        Parameters:
        -----------
        audio : ndarray
            Audio signal (mono)
        fs : int
            Sampling frequency (Hz)
        verbose : bool
            Print progress messages
            
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
        if verbose:
            print("Estimating F0...")
        F0 = self._estimate_f0(audio, fs)
        
        # Amplitude contour (lines 157-168)
        if verbose:
            print("Extracting amplitude contour...")
        A0 = self._extract_amplitude(audio, fs)
        
        # Resample for RPDE (do once, reuse)
        audio_resampled = scipy_signal.resample(audio, int(25000 * len(audio) / fs))
        
        # Initialize measures dictionary
        measures = {}
        
        if verbose:
            print(f"Computing features in parallel (workers={self.max_workers or 'auto'})...")
        
        # Define feature computation tasks
        # Each task is (name, function, dependencies)
        tasks = self._create_feature_tasks(audio, fs, F0, A0, audio_resampled)
        
        # Execute tasks in parallel
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            # Submit all tasks
            future_to_name = {}
            for name, func in tasks.items():
                future = executor.submit(func)
                future_to_name[future] = name
            
            # Collect results as they complete
            completed = 0
            total = len(tasks)
            for future in as_completed(future_to_name):
                name = future_to_name[future]
                try:
                    result = future.result()
                    if isinstance(result, dict):
                        measures.update(result)
                    else:
                        measures[name] = result
                    completed += 1
                    if verbose:
                        print(f"  [{completed}/{total}] {name} completed")
                except Exception as e:
                    if verbose:
                        print(f"  [{completed}/{total}] {name} FAILED: {e}")
                    warnings.warn(f"{name} computation failed: {e}")
                    completed += 1
        
        if verbose:
            print(f"Computed {len(measures)} measures")
        
        return measures, F0
    
    def _create_feature_tasks(self, audio, fs, F0, A0, audio_resampled):
        """
        Create dictionary of feature computation tasks
        
        Returns dict of {name: callable} where each callable computes a feature
        """
        tasks = {}
        
        # Determine workers for within-feature parallelization
        # Use fewer workers per feature when running multiple features in parallel
        within_workers = 2 if self.enable_within_feature_parallel and self.max_workers else None
        
        # Group A - Time series analysis (depends on F0, A0)
        tasks['Jitter'] = partial(compute_jitter_shimmer_features, F0, 'jitter')
        tasks['Shimmer'] = partial(compute_jitter_shimmer_features, A0, 'shimmer')
        tasks['PPE'] = partial(compute_ppe, F0, fs, 120)
        
        # Group B - Frequency domain analysis (depends on audio)
        # HNR with optional parallelization
        tasks['HNR/NHR'] = partial(compute_hnr_nhr, audio, fs, self.f0_min, self.f0_max,
                                   parallel=self.enable_within_feature_parallel,
                                   max_workers=within_workers)
        # GNE with optional parallelization
        tasks['GNE'] = partial(compute_gne, audio, fs,
                              parallel=self.enable_within_feature_parallel,
                              max_workers=within_workers)
        
        # Group C - Nonlinear dynamics (depends on audio)
        dfa_scaling = np.arange(50, 201, 20)
        tasks['DFA'] = partial(compute_dfa, audio, dfa_scaling)
        # RPDE with KD-tree optimization
        tasks['RPDE'] = partial(compute_rpde, audio_resampled, 4, 50, 0.12, 1000,
                               use_kdtree=self.use_rpde_kdtree)
        
        # Group D - Spectral analysis (depends on audio or F0)
        tasks['MFCC'] = partial(compute_mfcc_features, audio, fs)
        tasks['Wavelet'] = partial(compute_wavelet_features, F0)
        
        # Group E - Complex analysis (depends on audio)
        tasks['Glottal_Quotient'] = partial(compute_glottal_quotient, audio, fs, 
                                            self.f0_min, self.f0_max)
        tasks['VFER'] = partial(compute_vfer, audio, fs)
        tasks['EMD'] = partial(compute_emd_features, audio)
        
        return tasks
    
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


def analyze_voice_parallel(audio_path_or_array, fs=None, f0_min=50, f0_max=500, 
                          f0_algorithm='SWIPE', max_workers=None, verbose=True):
    """
    Convenience function for parallel voice analysis
    
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
    max_workers : int or None
        Number of parallel workers
    verbose : bool
        Print progress
        
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
    analyzer = VoiceAnalyzerParallel(f0_min=f0_min, f0_max=f0_max, 
                                     f0_algorithm=f0_algorithm, max_workers=max_workers)
    measures, F0 = analyzer.analyze(audio, fs, verbose=verbose)
    
    return measures, F0


def analyze_voice_file_parallel(filepath, f0_min=50, f0_max=500, f0_algorithm='SWIPE',
                                max_workers=None, verbose=True):
    """
    Convenience function to analyze a WAV file in parallel
    
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
    max_workers : int or None
        Number of parallel workers
    verbose : bool
        Print progress
        
    Returns:
    --------
    measures : dict
        All computed measures
    F0 : ndarray
        F0 contour
    """
    return analyze_voice_parallel(filepath, f0_min=f0_min, f0_max=f0_max, 
                                 f0_algorithm=f0_algorithm, max_workers=max_workers,
                                 verbose=verbose)


def analyze_batch_parallel(file_list, f0_min=50, f0_max=500, f0_algorithm='SWIPE',
                           max_workers=None, verbose=True):
    """
    Analyze multiple files in parallel
    
    Parameters:
    -----------
    file_list : list of str
        List of file paths
    f0_min : float
        Minimum F0 (Hz)
    f0_max : float
        Maximum F0 (Hz)
    f0_algorithm : str
        'SWIPE' or 'PRAAT'
    max_workers : int or None
        Number of parallel workers for file-level parallelization
    verbose : bool
        Print progress
        
    Returns:
    --------
    results : dict
        Dictionary mapping filename to (measures, F0) tuple
    """
    results = {}
    
    def process_file(filepath):
        """Process a single file"""
        try:
            measures, F0 = analyze_voice_file_parallel(
                filepath, f0_min, f0_max, f0_algorithm, 
                max_workers=1,  # Use 1 worker per file to avoid nested parallelism
                verbose=False
            )
            return filepath, (measures, F0), None
        except Exception as e:
            return filepath, None, str(e)
    
    if verbose:
        print(f"Processing {len(file_list)} files in parallel...")
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_file, f) for f in file_list]
        
        completed = 0
        for future in as_completed(futures):
            filepath, result, error = future.result()
            completed += 1
            
            if error:
                if verbose:
                    print(f"  [{completed}/{len(file_list)}] {filepath} FAILED: {error}")
                results[filepath] = None
            else:
                if verbose:
                    print(f"  [{completed}/{len(file_list)}] {filepath} completed")
                results[filepath] = result
    
    if verbose:
        successful = sum(1 for v in results.values() if v is not None)
        print(f"Completed: {successful}/{len(file_list)} successful")
    
    return results

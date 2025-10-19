"""
R Interface for Voice Analysis Toolbox

Optimized for R reticulate environment:
- Handles R data types (character, numeric, logical, integer)
- Process-based parallelization (GIL-safe for R)
- Data.frame-compatible output
- Memory-efficient batch processing

Usage in R:
-----------
library(reticulate)
va <- import("voice_analysis.r_interface")

# Single file
result <- va$analyze_for_r("audio.wav", n_cores=8L)

# Batch processing
files <- c("audio1.wav", "audio2.wav")
df <- va$analyze_batch_for_r(files, n_cores=30L, return_dataframe=TRUE)
results_df <- as.data.frame(df)
"""

import numpy as np
import soundfile as sf
import os
import warnings
from typing import Union, List, Dict, Any, Optional


def analyze_for_r(
    audio_path_or_array: Union[str, np.ndarray],
    fs: Optional[int] = None,
    features: Union[str, List[str]] = 'all',
    n_cores: int = 1,
    verbose: bool = True,
    use_cython: bool = True,
    timeout: Optional[float] = None
) -> Dict[str, Any]:
    """
    R-optimized interface for voice analysis
    
    Parameters:
    -----------
    audio_path_or_array : str or ndarray
        File path (R character) or audio array (R numeric vector)
    fs : int or None
        Sampling rate (required if audio_path_or_array is array)
    features : str or list
        'all' or list of specific features to compute
    n_cores : int
        Number of CPU cores to use (R integer)
    verbose : bool
        Print progress messages (R logical)
    use_cython : bool
        Use Cython-optimized functions if available
    timeout : float or None
        Maximum time in seconds for analysis
    
    Returns:
    --------
    dict
        {
            'measures': dict of measure_name -> value,
            'f0': list of F0 values (or None),
            'fs': sampling rate,
            'success': bool,
            'error': str or None
        }
    """
    try:
        # Handle R input types
        if isinstance(audio_path_or_array, str):
            # R character string - file path
            if not os.path.exists(audio_path_or_array):
                return {
                    'measures': {},
                    'f0': None,
                    'fs': None,
                    'success': False,
                    'error': f"File not found: {audio_path_or_array}"
                }
            
            audio, fs = sf.read(audio_path_or_array)
            
        else:
            # R numeric vector - convert to NumPy
            audio = np.asarray(audio_path_or_array, dtype=np.float64)
            if fs is None:
                return {
                    'measures': {},
                    'f0': None,
                    'fs': None,
                    'success': False,
                    'error': "Sampling rate (fs) required when passing audio array"
                }
        
        # Ensure audio is 1D
        if audio.ndim > 1:
            audio = audio[:, 0]  # Take first channel
        
        # Configure analyzer based on n_cores and Cython availability
        if use_cython:
            try:
                from voice_analysis.features import rpde_cython
                cython_available = True
            except ImportError:
                cython_available = False
                if verbose:
                    print("Cython extensions not available, using Numba")
        else:
            cython_available = False
        
        # Choose appropriate analyzer
        if n_cores > 1:
            from voice_analysis.core_parallel import VoiceAnalyzerParallel
            analyzer = VoiceAnalyzerParallel(
                max_workers=n_cores,
                enable_within_feature_parallel=False,  # Avoid over-subscription
                use_rpde_kdtree=False,  # KD-tree not faster currently
                timeout=timeout
            )
        else:
            from voice_analysis.core import VoiceAnalyzer
            analyzer = VoiceAnalyzer()
        
        # Perform analysis
        with warnings.catch_warnings():
            if not verbose:
                warnings.simplefilter("ignore")
            
            measures, f0 = analyzer.analyze(audio, fs)
        
        # Convert to R-friendly format
        result = {
            'measures': {k: float(v) if not np.isnan(v) else None 
                        for k, v in measures.items()},
            'f0': f0.tolist() if f0 is not None else None,
            'fs': int(fs),
            'success': True,
            'error': None
        }
        
        return result
        
    except Exception as e:
        import traceback
        error_msg = f"{type(e).__name__}: {str(e)}"
        
        if verbose:
            print(f"Error during analysis: {error_msg}")
            traceback.print_exc()
        
        return {
            'measures': {},
            'f0': None,
            'fs': None,
            'success': False,
            'error': error_msg
        }


def analyze_batch_for_r(
    file_paths: List[str],
    n_cores: int = 1,
    verbose: bool = True,
    return_dataframe: bool = True,
    chunk_size: int = 10,
    timeout_per_file: Optional[float] = None,
    on_error: str = 'warn'
) -> Union[Dict[str, List], List[Dict]]:
    """
    Batch analysis optimized for R
    
    Parameters:
    -----------
    file_paths : list of str
        R character vector of file paths
    n_cores : int
        Number of CPU cores to use
    verbose : bool
        Print progress (R logical)
    return_dataframe : bool
        If True, return dict suitable for R data.frame conversion
    chunk_size : int
        Process files in chunks (for memory management)
    timeout_per_file : float or None
        Maximum time per file
    on_error : str
        'warn', 'raise', or 'skip'
    
    Returns:
    --------
    dict or list
        If return_dataframe=True: dict with column_name -> [values]
        If return_dataframe=False: list of result dicts
    """
    from joblib import Parallel, delayed
    import psutil
    
    # Validate inputs
    file_paths = [str(f) for f in file_paths]
    n_files = len(file_paths)
    
    if n_files == 0:
        if return_dataframe:
            return {}
        else:
            return []
    
    # Optimize worker count
    if n_cores <= 0:
        n_cores = psutil.cpu_count(logical=False) - 1
    n_cores = min(n_cores, n_files)
    
    if verbose:
        print(f"Processing {n_files} files with {n_cores} workers")
    
    # Set environment flag for nested parallelization control
    os.environ['VOICE_ANALYSIS_BATCH'] = '1'
    
    # Configure NumPy threading
    _set_numpy_threads(max(1, psutil.cpu_count(logical=False) // n_cores))
    
    try:
        # Process in chunks to manage memory
        all_results = []
        n_chunks = (n_files + chunk_size - 1) // chunk_size
        
        for chunk_idx in range(n_chunks):
            start_idx = chunk_idx * chunk_size
            end_idx = min(start_idx + chunk_size, n_files)
            chunk_files = file_paths[start_idx:end_idx]
            
            if verbose:
                print(f"Processing chunk {chunk_idx + 1}/{n_chunks} ({len(chunk_files)} files)")
            
            # Parallel processing with joblib (process-based, R-safe)
            chunk_results = Parallel(
                n_jobs=n_cores,
                backend='loky',  # Process-based parallelism
                verbose=10 if verbose else 0
            )(
                delayed(_analyze_single_file_safe)(
                    fpath, timeout_per_file, on_error
                )
                for fpath in chunk_files
            )
            
            all_results.extend(chunk_results)
            
            # Cleanup memory after each chunk
            if chunk_idx < n_chunks - 1:
                import gc
                gc.collect()
        
    finally:
        # Clear batch flag
        os.environ['VOICE_ANALYSIS_BATCH'] = '0'
    
    # Format results
    if return_dataframe:
        return _results_to_dataframe_dict(all_results, file_paths, on_error)
    else:
        return all_results


def _analyze_single_file_safe(
    filepath: str,
    timeout: Optional[float] = None,
    on_error: str = 'warn'
) -> Dict[str, Any]:
    """
    Analyze single file with error handling
    """
    try:
        result = analyze_for_r(
            filepath,
            n_cores=1,  # Already parallelized at file level
            verbose=False,
            timeout=timeout
        )
        result['file'] = filepath
        return result
        
    except Exception as e:
        error_msg = f"{type(e).__name__}: {str(e)}"
        
        if on_error == 'raise':
            raise
        elif on_error == 'warn':
            warnings.warn(f"Error processing {filepath}: {error_msg}")
        
        return {
            'file': filepath,
            'measures': {},
            'f0': None,
            'fs': None,
            'success': False,
            'error': error_msg
        }


def _results_to_dataframe_dict(
    results: List[Dict],
    file_paths: List[str],
    on_error: str
) -> Dict[str, List]:
    """
    Convert results to R data.frame-compatible dict
    
    Returns dict with:
    - 'file': [file paths]
    - measure_name: [values]
    - 'success': [bool]
    - 'error': [str or None]
    """
    import pandas as pd
    
    # Initialize with file paths
    df_dict = {'file': file_paths}
    
    # Collect all unique measure names
    all_measures = set()
    for result in results:
        if result['success'] and result['measures']:
            all_measures.update(result['measures'].keys())
    
    # Initialize columns with NA
    for measure in sorted(all_measures):
        df_dict[measure] = [None] * len(file_paths)
    
    df_dict['success'] = [False] * len(file_paths)
    df_dict['error'] = [None] * len(file_paths)
    
    # Fill in values
    for i, result in enumerate(results):
        df_dict['success'][i] = result.get('success', False)
        df_dict['error'][i] = result.get('error', None)
        
        if result['success'] and result['measures']:
            for measure, value in result['measures'].items():
                if measure in df_dict:
                    df_dict[measure][i] = value
    
    return df_dict


def _set_numpy_threads(n_threads: int):
    """
    Set NumPy/BLAS thread count to avoid over-subscription
    """
    n_threads = max(1, n_threads)
    n_threads_str = str(n_threads)
    
    os.environ['OMP_NUM_THREADS'] = n_threads_str
    os.environ['OPENBLAS_NUM_THREADS'] = n_threads_str
    os.environ['MKL_NUM_THREADS'] = n_threads_str
    os.environ['VECLIB_MAXIMUM_THREADS'] = n_threads_str
    os.environ['NUMEXPR_NUM_THREADS'] = n_threads_str


def get_available_features() -> List[str]:
    """
    Get list of all available features
    
    Returns:
    --------
    list of str : Feature names
    """
    return [
        'jitter', 'shimmer', 'hnr', 'nhr', 'rpde', 'dfa', 'ppe',
        'gne', 'vfer', 'mfcc', 'wavelet', 'emd', 'gq'
    ]


def check_cython_available() -> bool:
    """
    Check if Cython-optimized extensions are available
    
    Returns:
    --------
    bool : True if Cython extensions installed
    """
    try:
        from voice_analysis.features import rpde_cython
        from voice_analysis.utils import perturbation_cython
        return True
    except ImportError:
        return False


def get_system_info() -> Dict[str, Any]:
    """
    Get system information for optimization
    
    Returns:
    --------
    dict
        {
            'cpu_count': int,
            'cpu_count_physical': int,
            'platform': str,
            'machine': str,
            'cython_available': bool,
            'numba_available': bool,
            'recommended_workers': int
        }
    """
    import psutil
    import platform
    
    cpu_count = psutil.cpu_count(logical=True)
    cpu_physical = psutil.cpu_count(logical=False)
    
    # Try to detect Numba
    try:
        import numba
        numba_available = True
    except ImportError:
        numba_available = False
    
    # Recommended workers (leave 1-2 cores free)
    recommended = max(1, cpu_physical - 1) if cpu_physical else max(1, cpu_count - 2)
    
    return {
        'cpu_count': cpu_count,
        'cpu_count_physical': cpu_physical,
        'platform': platform.system(),
        'machine': platform.machine(),
        'cython_available': check_cython_available(),
        'numba_available': numba_available,
        'recommended_workers': recommended
    }

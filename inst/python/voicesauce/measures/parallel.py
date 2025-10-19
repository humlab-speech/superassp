"""
Multi-core parallel processing for voice quality measures
Optimized for Apple Silicon M1 (10 cores) and AMD EPYC (32 cores)
"""

import numpy as np
from multiprocessing import Pool, cpu_count
from functools import partial
from typing import Tuple, List
import warnings


def get_optimal_num_processes() -> int:
    """
    Determine optimal number of processes based on CPU architecture
    
    Returns:
        Optimal number of worker processes
    """
    total_cores = cpu_count()
    
    # Apple Silicon optimization
    # M1/M2: 8-10 cores (4 performance + 4-6 efficiency)
    # Use performance cores primarily, leave some for OS
    if total_cores >= 8:
        # Use 75% of cores for compute-heavy tasks
        return max(4, int(total_cores * 0.75))
    elif total_cores >= 4:
        # VM with 4 cores: use 3 workers
        return total_cores - 1
    else:
        # Low core count: use all but one
        return max(1, total_cores - 1)


class ParallelProcessor:
    """Parallel processing manager for voice quality measures"""
    
    def __init__(self, n_processes: int = None):
        """
        Initialize parallel processor
        
        Args:
            n_processes: Number of processes (None = auto-detect)
        """
        if n_processes is None:
            n_processes = get_optimal_num_processes()
        
        self.n_processes = n_processes
        self._pool = None
    
    def __enter__(self):
        """Context manager entry"""
        self._pool = Pool(processes=self.n_processes)
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """Context manager exit"""
        if self._pool is not None:
            self._pool.close()
            self._pool.join()
    
    def map(self, func, items, chunksize=None):
        """
        Parallel map operation
        
        Args:
            func: Function to apply
            items: Items to process
            chunksize: Chunk size for processing
        
        Returns:
            List of results
        """
        if self._pool is None:
            raise RuntimeError("Use ParallelProcessor as context manager")
        
        if chunksize is None:
            # Auto-determine chunk size
            chunksize = max(1, len(items) // (self.n_processes * 4))
        
        return self._pool.map(func, items, chunksize=chunksize)


def _process_harmonic_chunk(args: Tuple) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Process a chunk of frames for harmonic extraction
    Worker function for multiprocessing
    
    Args:
        args: (audio, fs, f0_values, sample_shift, n_periods, start_idx, end_idx)
    
    Returns:
        h1, h2, h4 arrays for this chunk
    """
    audio, fs, f0_values, sample_shift, n_periods, start_idx, end_idx = args
    
    # Import here to avoid issues with multiprocessing
    from .harmonics import extract_harmonic_amplitude
    
    chunk_size = end_idx - start_idx
    h1 = np.full(chunk_size, np.nan)
    h2 = np.full(chunk_size, np.nan)
    h4 = np.full(chunk_size, np.nan)
    
    for i, k in enumerate(range(start_idx, end_idx)):
        f0 = f0_values[k]
        
        if np.isnan(f0) or f0 <= 0:
            continue
        
        # Calculate segment indices
        ks = k * sample_shift
        n0 = fs / f0
        
        start = int(ks - n_periods/2 * n0)
        end = int(ks + n_periods/2 * n0)
        
        # Check bounds
        if start < 0 or end >= len(audio):
            continue
        
        # Extract segment
        segment = audio[start:end]
        
        if len(segment) == 0:
            continue
        
        # Extract harmonics
        try:
            h1[i], _ = extract_harmonic_amplitude(segment, f0, fs)
            h2[i], _ = extract_harmonic_amplitude(segment, 2*f0, fs)
            h4[i], _ = extract_harmonic_amplitude(segment, 4*f0, fs)
        except:
            continue
    
    return h1, h2, h4, start_idx


def get_harmonics_parallel(audio: np.ndarray, fs: int, f0_values: np.ndarray,
                           frame_shift: float = 1.0, n_periods: int = 3,
                           n_processes: int = None) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Parallel harmonic extraction using multiprocessing
    
    Optimized for multi-core CPUs (Apple Silicon M1, AMD EPYC)
    
    Args:
        audio: Full audio signal
        fs: Sampling rate
        f0_values: F0 array
        frame_shift: Frame shift in milliseconds
        n_periods: Number of pitch periods
        n_processes: Number of processes (None = auto)
    
    Returns:
        h1, h2, h4: Harmonic amplitudes
    """
    sample_shift = int(fs * frame_shift / 1000)
    n_frames = len(f0_values)
    
    # For small files, don't use multiprocessing (overhead not worth it)
    if n_frames < 1000:
        from .harmonics import get_harmonics
        return get_harmonics(audio, fs, f0_values, frame_shift, n_periods)
    
    # Determine number of processes
    if n_processes is None:
        n_processes = get_optimal_num_processes()
    
    # Split frames into chunks
    chunk_size = max(100, n_frames // n_processes)
    chunks = []
    
    for i in range(0, n_frames, chunk_size):
        end_idx = min(i + chunk_size, n_frames)
        chunks.append((audio, fs, f0_values, sample_shift, n_periods, i, end_idx))
    
    # Process chunks in parallel
    with ParallelProcessor(n_processes) as processor:
        results = processor.map(_process_harmonic_chunk, chunks)
    
    # Combine results
    h1 = np.full(n_frames, np.nan)
    h2 = np.full(n_frames, np.nan)
    h4 = np.full(n_frames, np.nan)
    
    for h1_chunk, h2_chunk, h4_chunk, start_idx in results:
        end_idx = start_idx + len(h1_chunk)
        h1[start_idx:end_idx] = h1_chunk
        h2[start_idx:end_idx] = h2_chunk
        h4[start_idx:end_idx] = h4_chunk
    
    return h1, h2, h4


def _process_cpp_chunk(args: Tuple) -> Tuple[np.ndarray, int]:
    """
    Process a chunk of frames for CPP calculation
    Worker function for multiprocessing
    """
    audio, fs, f0_values, sample_shift, n_periods, start_idx, end_idx = args
    
    from scipy.fft import fft, ifft
    from scipy.signal import find_peaks
    
    chunk_size = end_idx - start_idx
    cpp = np.full(chunk_size, np.nan)
    n_ms = int(fs / 1000)
    
    # Window cache per worker
    window_cache = {}
    
    for i, k in enumerate(range(start_idx, end_idx)):
        f0 = f0_values[k]
        
        if np.isnan(f0) or f0 <= 0:
            continue
        
        ks = int(k * sample_shift)
        n0 = fs / f0
        
        start = int(ks - n_periods/2 * n0)
        end = int(ks + n_periods/2 * n0)
        
        if start < 0 or end >= len(audio):
            continue
        
        segment = audio[start:end]
        if len(segment) < 2:
            continue
        
        # Get or create window
        seg_len = len(segment)
        if seg_len not in window_cache:
            window_cache[seg_len] = np.hamming(seg_len)
        
        windowed = segment * window_cache[seg_len]
        
        # Cepstrum
        spectrum = fft(windowed)
        log_spectrum = np.log(np.abs(spectrum) + 1e-10)
        cepstrum = ifft(log_spectrum)
        cepstrum_db = 10 * np.log10(np.abs(cepstrum)**2 + 1e-10)
        
        half_len = len(cepstrum_db) // 2
        cepstrum_db = cepstrum_db[:half_len]
        
        if n_ms >= half_len:
            continue
        
        search_region = cepstrum_db[n_ms:]
        if len(search_region) < 2:
            continue
        
        # Find peaks
        peaks, _ = find_peaks(search_region, distance=int(n0))
        if len(peaks) == 0:
            continue
        
        expected_idx = min(int(n0), len(search_region) - 1)
        closest_peak_idx = np.argmin(np.abs(peaks - expected_idx))
        peak_idx = peaks[closest_peak_idx]
        peak_val = search_region[peak_idx]
        
        # Baseline
        x = np.arange(len(search_region))
        coeffs = np.polyfit(x, search_region, 1)
        baseline_val = np.polyval(coeffs, peak_idx)
        
        cpp[i] = peak_val - baseline_val
    
    return cpp, start_idx


def get_cpp_parallel(audio: np.ndarray, fs: int, f0_values: np.ndarray,
                     frame_shift: float = 1.0, n_periods: int = 5,
                     n_processes: int = None) -> np.ndarray:
    """
    Parallel CPP calculation
    
    Args:
        audio: Full audio signal
        fs: Sampling rate
        f0_values: F0 array
        frame_shift: Frame shift in milliseconds
        n_periods: Number of pitch periods
        n_processes: Number of processes (None = auto)
    
    Returns:
        cpp: CPP values
    """
    sample_shift = int(fs * frame_shift / 1000)
    n_frames = len(f0_values)
    
    # For small files, don't use multiprocessing
    if n_frames < 1000:
        from .cpp_optimized import get_cpp_optimized
        return get_cpp_optimized(audio, fs, f0_values, frame_shift, n_periods)
    
    # Determine number of processes
    if n_processes is None:
        n_processes = get_optimal_num_processes()
    
    # Split frames into chunks
    chunk_size = max(100, n_frames // n_processes)
    chunks = []
    
    for i in range(0, n_frames, chunk_size):
        end_idx = min(i + chunk_size, n_frames)
        chunks.append((audio, fs, f0_values, sample_shift, n_periods, i, end_idx))
    
    # Process chunks in parallel
    with ParallelProcessor(n_processes) as processor:
        results = processor.map(_process_cpp_chunk, chunks)
    
    # Combine results
    cpp = np.full(n_frames, np.nan)
    
    for cpp_chunk, start_idx in results:
        end_idx = start_idx + len(cpp_chunk)
        cpp[start_idx:end_idx] = cpp_chunk
    
    return cpp


# Auto-detection: use parallel version if enough cores
PARALLEL_THRESHOLD_CORES = 4

def should_use_parallel() -> bool:
    """Determine if parallel processing should be used"""
    return cpu_count() >= PARALLEL_THRESHOLD_CORES

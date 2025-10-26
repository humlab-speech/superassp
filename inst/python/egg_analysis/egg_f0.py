"""
EGG F0 and Oq Analysis

Core analysis function for EGG signals, adapted from egg_python.
Returns results at equal intervals for superassp SSFF tracks.
"""

import numpy as np
from typing import Tuple, Optional
import sys
from pathlib import Path

# Import egg_python modules
egg_python_path = Path(__file__).parent.parent.parent.parent.parent / "egg" / "egg_python"
if egg_python_path.exists():
    sys.path.insert(0, str(egg_python_path))

    from egg.peakdet.fo import fo
    from egg.utils.signal_processing import deviate, smoo
else:
    raise ImportError(
        "egg_python package not found. Please ensure egg_python is in the expected location."
    )


def analyze_egg_f0(
    audio_array: np.ndarray,
    sample_rate: int,
    method: int = 3,
    smoothing: int = 3,
    max_f0: float = 500.0,
    resample_coef: int = 100,
    propthresh: float = 0.5,
    frame_shift_ms: float = 10.0
) -> dict:
    """
    Analyze EGG signal and extract f0 and Oq at equal intervals.

    This function performs EGG analysis using the validated egg_python algorithms
    and returns results suitable for superassp SSFF tracks (equal time intervals).

    Parameters
    ----------
    audio_array : np.ndarray
        Audio samples from av package (float64, normalized to [-1, 1])
    sample_rate : int
        Sample rate in Hz
    method : int, default=3
        Peak handling method:
        - 0: Highest peak
        - 1: First peak
        - 2: Last peak
        - 3: Barycentre (weighted average) - RECOMMENDED
        - 4: Exclude double peaks
    smoothing : int, default=3
        Smoothing step for dEGG signal (larger = more smoothing)
    max_f0 : float, default=500.0
        Maximum plausible f0 in Hz
    resample_coef : int, default=100
        Resampling coefficient for sub-sample accuracy
    propthresh : float, default=0.5
        Threshold for double peak detection
    frame_shift_ms : float, default=10.0
        Frame shift in milliseconds for output track

    Returns
    -------
    dict
        Dictionary with keys:
        - 'f0': np.ndarray of f0 values at equal intervals (Hz)
        - 'oq': np.ndarray of open quotient values (percentage)
        - 'times': np.ndarray of time points (seconds)
        - 'voicing': np.ndarray of voicing flags (1=voiced, 0=unvoiced)
        - 'raw_f0': np.ndarray of f0 for each detected cycle
        - 'raw_oq': np.ndarray of Oq for each detected cycle
        - 'raw_times': np.ndarray of cycle times (seconds)
        - 'sample_rate': Original sample rate
        - 'n_cycles': Number of glottal cycles detected

    Notes
    -----
    The function operates in two stages:

    1. **Cycle-based analysis** (irregular timing):
       Uses egg_python's fo() function to detect glottal cycles and
       compute f0 and Oq for each cycle.

    2. **Interpolation to equal intervals**:
       Converts cycle-based measurements to equal time intervals
       using linear interpolation, suitable for SSFF tracks.

    For unvoiced frames or regions without valid cycles, f0 is set to 0
    and Oq is set to NaN.

    References
    ----------
    - Mazaudon & Michaud (2008). Tonal contrasts and initial consonants
    - Michaud (2004). Final consonants and glottalization
    - Henrich et al. (2004). On the use of the derivative of EGG signals
    """

    # Validate inputs
    if not isinstance(audio_array, np.ndarray):
        audio_array = np.array(audio_array, dtype=np.float64)

    if audio_array.ndim > 1:
        # Take first channel if multi-channel
        audio_array = audio_array[:, 0]

    # Ensure float64 for numerical precision
    audio_array = audio_array.astype(np.float64)

    # Set up analysis parameters (egg_python format)
    coef = np.array([sample_rate, smoothing, 1, 0], dtype=float)

    # Run egg_python f0 analysis
    try:
        (f0_cycles, oq, oqval, deopa, goodperiods,
         oq_s, oqval_s, deopa_s, goodperiods_s,
         simppeak, dsig, s_dsig, ddsig, s_ddsig) = fo(
            coef, method, propthresh, resample_coef, max_f0, audio_array, sample_rate
        )
    except Exception as e:
        raise RuntimeError(f"EGG analysis failed: {str(e)}")

    n_cycles = len(f0_cycles)

    if n_cycles == 0:
        # No cycles detected - return empty results
        duration = len(audio_array) / sample_rate
        n_frames = int(np.ceil(duration / (frame_shift_ms / 1000.0)))

        return {
            'f0': np.zeros(n_frames, dtype=np.float32),
            'oq': np.full(n_frames, np.nan, dtype=np.float32),
            'times': np.arange(n_frames) * (frame_shift_ms / 1000.0),
            'voicing': np.zeros(n_frames, dtype=np.int16),
            'raw_f0': np.array([], dtype=np.float32),
            'raw_oq': np.array([], dtype=np.float32),
            'raw_times': np.array([], dtype=np.float32),
            'sample_rate': sample_rate,
            'n_cycles': 0
        }

    # Extract cycle times from simppeak (closing peak times in samples)
    # simppeak is [time_in_samples, amplitude]
    cycle_times_samples = simppeak[:, 0]
    cycle_times_sec = cycle_times_samples / sample_rate

    # Use oqval (open quotient from peak detection method - more accurate)
    cycle_oq = oqval.copy()

    # Create equal-interval time grid
    duration = len(audio_array) / sample_rate
    frame_shift_sec = frame_shift_ms / 1000.0
    n_frames = int(np.ceil(duration / frame_shift_sec))
    frame_times = np.arange(n_frames) * frame_shift_sec

    # Interpolate f0 and Oq to equal intervals
    f0_interp = np.zeros(n_frames, dtype=np.float32)
    oq_interp = np.full(n_frames, np.nan, dtype=np.float32)
    voicing = np.zeros(n_frames, dtype=np.int16)

    if n_cycles > 0:
        # For each frame, find the nearest cycle
        # We use simple nearest-neighbor interpolation for voiced/unvoiced decision
        # and linear interpolation for f0/Oq values

        for i, frame_time in enumerate(frame_times):
            # Find cycles within +/- half frame shift
            half_shift = frame_shift_sec / 2.0
            nearby_cycles = np.where(
                np.abs(cycle_times_sec - frame_time) <= half_shift
            )[0]

            if len(nearby_cycles) > 0:
                # Use nearest cycle
                nearest_idx = nearby_cycles[np.argmin(
                    np.abs(cycle_times_sec[nearby_cycles] - frame_time)
                )]

                f0_interp[i] = f0_cycles[nearest_idx]

                # Only use Oq if it's valid (> 0)
                if cycle_oq[nearest_idx] > 0:
                    oq_interp[i] = cycle_oq[nearest_idx]

                voicing[i] = 1
            else:
                # No nearby cycle - unvoiced frame
                f0_interp[i] = 0.0
                oq_interp[i] = np.nan
                voicing[i] = 0

    return {
        'f0': f0_interp,
        'oq': oq_interp,
        'times': frame_times,
        'voicing': voicing,
        'raw_f0': f0_cycles.astype(np.float32),
        'raw_oq': cycle_oq.astype(np.float32),
        'raw_times': cycle_times_sec.astype(np.float32),
        'sample_rate': sample_rate,
        'n_cycles': n_cycles
    }


def analyze_egg_batch(
    audio_arrays: list,
    sample_rates: list,
    method: int = 3,
    smoothing: int = 3,
    max_f0: float = 500.0,
    frame_shift_ms: float = 10.0
) -> list:
    """
    Batch process multiple EGG files.

    Parameters
    ----------
    audio_arrays : list of np.ndarray
        List of audio samples
    sample_rates : list of int
        List of sample rates
    method : int
        Peak handling method (see analyze_egg_f0)
    smoothing : int
        Smoothing parameter
    max_f0 : float
        Maximum f0
    frame_shift_ms : float
        Frame shift in milliseconds

    Returns
    -------
    list of dict
        List of analysis results
    """
    results = []

    for audio, sr in zip(audio_arrays, sample_rates):
        result = analyze_egg_f0(
            audio, sr, method, smoothing, max_f0,
            frame_shift_ms=frame_shift_ms
        )
        results.append(result)

    return results

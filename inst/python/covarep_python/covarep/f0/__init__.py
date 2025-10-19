"""
F0 (Fundamental Frequency) Estimation Module

Implements various F0 tracking algorithms from COVAREP including:
- SRH (Summation of Residual Harmonics)
- And others to be added

Based on original COVAREP MATLAB implementations
"""

import numpy as np
from scipy import signal
from scipy.fft import fft, ifft
import warnings

__all__ = ['F0Tracker', 'pitch_srh']


class F0Tracker:
    """
    Unified F0 tracking interface
    
    Supports multiple F0 estimation algorithms with a common API
    
    Parameters
    ----------
    method : str
        F0 tracking method: 'srh' (default)
    f0_min : float
        Minimum F0 in Hz (default: 50)
    f0_max : float
        Maximum F0 in Hz (default: 500)
    hop_size : float
        Hop size in milliseconds (default: 5.0)
    
    Examples
    --------
    >>> tracker = F0Tracker(method='srh', f0_min=50, f0_max=400)
    >>> f0, vuv = tracker.estimate(audio, fs=16000)
    """
    
    def __init__(self, method='srh', f0_min=50.0, f0_max=500.0, hop_size=5.0):
        self.method = method.lower()
        self.f0_min = f0_min
        self.f0_max = f0_max
        self.hop_size = hop_size  # in ms
        
        self._validate_parameters()
    
    def _validate_parameters(self):
        """Validate F0 tracker parameters"""
        if self.f0_min <= 0:
            raise ValueError("f0_min must be positive")
        if self.f0_max <= self.f0_min:
            raise ValueError("f0_max must be greater than f0_min")
        if self.hop_size <= 0:
            raise ValueError("hop_size must be positive")
    
    def estimate(self, audio, fs):
        """
        Estimate F0 from audio signal
        
        Parameters
        ----------
        audio : array_like
            Input audio signal
        fs : float
            Sampling frequency in Hz
            
        Returns
        -------
        f0 : ndarray
            F0 estimates in Hz
        vuv : ndarray
            Voice/unvoiced decisions (boolean)
        times : ndarray
            Time stamps for each estimate
        """
        audio = np.asarray(audio).ravel()
        
        if self.method == 'srh':
            return pitch_srh(audio, fs, self.f0_min, self.f0_max, self.hop_size)
        else:
            raise ValueError(f"Unknown F0 tracking method: {self.method}")


def pitch_srh(wave, fs, f0min=50.0, f0max=500.0, hopsize=5.0):
    """
    SRH (Summation of Residual Harmonics) F0 tracker
    
    Robust pitch tracker that exploits harmonics and subharmonics of the
    residual excitation signal. Uses 2-iteration refinement like MATLAB.
    
    Reference:
    T. Drugman and A. Alwan, "Joint Robust Voicing Detection and Pitch 
    Estimation Based on Residual Harmonics", Interspeech 2011
    
    Parameters
    ----------
    wave : array_like
        Input speech signal
    fs : float
        Sampling frequency in Hz
    f0min : float
        Minimum F0 in Hz (default: 50)
    f0max : float
        Maximum F0 in Hz (default: 500)
    hopsize : float
        Hop size in milliseconds (default: 5.0)
        
    Returns
    -------
    f0 : ndarray
        F0 estimates in Hz
    vuv : ndarray
        Voice/unvoiced decisions (boolean array)
    srh_values : ndarray
        SRH criterion values
    times : ndarray
        Time stamps in seconds
    """
    wave = np.asarray(wave).ravel()
    
    # Resample to 16kHz if higher (matching MATLAB)
    if fs > 16000:
        import warnings
        warnings.warn("Sample rate > 16kHz. Resampling to 16kHz.")
        from scipy import signal as sp_signal
        wave = sp_signal.resample_poly(wave, 16000, int(fs))
        fs = 16000
    
    # Parameters (matching MATLAB)
    hop_samples = int(hopsize * fs / 1000)
    n_harmonics = 5
    srh_iter_thresh = 0.1
    srh_std_thresh = 0.05
    voicing_thresh = 0.07
    voicing_thresh2 = 0.085
    n_iter = 2
    
    # Pre-emphasis to enhance high frequencies
    wave_preemph = signal.lfilter([1, -0.97], [1], wave)
    
    # Get LP residual using LPC analysis
    residual = _get_lp_residual(wave_preemph, fs)
    
    # Create spectrogram for SRH computation (matching MATLAB)
    f0, vuv, srh_values, times = _srh_criterion_matlab_style(
        residual, fs, f0min, f0max, hop_samples, n_harmonics, 
        srh_iter_thresh, srh_std_thresh, voicing_thresh, voicing_thresh2, n_iter
    )
    
    return f0, vuv, srh_values, times


def _get_lp_residual(wave, fs, order=None):
    """
    Compute LP residual using Linear Predictive Coding
    
    Parameters
    ----------
    wave : ndarray
        Input signal
    fs : float
        Sampling frequency
    order : int, optional
        LPC order (default: fs/1000 + 2)
        
    Returns
    -------
    residual : ndarray
        LP residual signal
    """
    if order is None:
        order = int(fs / 1000) + 2
    
    # Frame parameters
    frame_len = int(0.025 * fs)  # 25 ms
    hop_len = int(0.005 * fs)    # 5 ms
    
    # Initialize residual
    residual = np.zeros_like(wave)
    
    # Process frame by frame
    n_frames = 1 + (len(wave) - frame_len) // hop_len
    
    for i in range(n_frames):
        start = i * hop_len
        end = start + frame_len
        
        if end > len(wave):
            break
        
        frame = wave[start:end]
        
        # Apply window
        window = np.hamming(frame_len)
        frame_windowed = frame * window
        
        # Compute LPC coefficients
        # Using autocorrelation method (Levinson-Durbin)
        r = np.correlate(frame_windowed, frame_windowed, mode='full')
        r = r[len(r)//2:]  # Keep positive lags
        
        # Levinson-Durbin recursion
        a = _levinson_durbin(r, order)
        
        # Compute residual for this frame
        frame_residual = signal.lfilter(a, [1], frame)
        
        # Overlap-add
        residual[start:end] += frame_residual
    
    return residual


def _levinson_durbin(r, order):
    """
    Levinson-Durbin recursion for LPC coefficient computation
    
    Parameters
    ----------
    r : ndarray
        Autocorrelation sequence
    order : int
        LPC order
        
    Returns
    -------
    a : ndarray
        LPC coefficients [1, -a1, -a2, ..., -an]
    """
    a = np.zeros(order + 1)
    a[0] = 1.0
    
    if len(r) < order + 1:
        # Not enough data
        return a
    
    # Initialize
    e = r[0]
    
    for i in range(1, order + 1):
        # Reflection coefficient
        lambda_i = -np.sum(a[1:i] * r[i-1:0:-1]) / e
        lambda_i -= r[i] / e
        
        # Update coefficients
        a_prev = a.copy()
        a[i] = lambda_i
        for j in range(1, i):
            a[j] = a_prev[j] + lambda_i * a_prev[i-j]
        
        # Update error
        e *= (1 - lambda_i**2)
        
        if e <= 0:
            break
    
    return a


def _srh_criterion_matlab_style(residual, fs, f0min_orig, f0max_orig, hop_samples, 
                                 n_harmonics, srh_iter_thresh, srh_std_thresh,
                                 voicing_thresh, voicing_thresh2, n_iter):
    """
    Compute SRH criterion with 2-iteration refinement (matching MATLAB)
    
    Parameters
    ----------
    residual : ndarray
        LP residual signal
    fs : float
        Sampling frequency
    f0min_orig : float
        Initial minimum F0
    f0max_orig : float
        Initial maximum F0
    hop_samples : int
        Hop size in samples
    n_harmonics : int
        Number of harmonics to sum
    srh_iter_thresh : float
        Threshold for iteration refinement
    srh_std_thresh : float
        Threshold for voicing decision adaptation
    voicing_thresh : float
        Primary voicing threshold
    voicing_thresh2 : float
        Secondary voicing threshold
    n_iter : int
        Number of iterations
        
    Returns
    -------
    f0 : ndarray
        F0 estimates
    vuv : ndarray
        Voice/unvoiced decisions
    srh_values : ndarray
        SRH criterion values
    times : ndarray
        Time stamps
    """
    # Frame parameters (matching MATLAB: 100ms - 2)
    frame_len = int(0.100 * fs) - 2
    # Make frame length even
    frame_len = (frame_len // 2) * 2
    half_dur = frame_len // 2
    
    # Time points
    n_samples = len(residual)
    times_samples = np.arange(half_dur, n_samples - half_dur, hop_samples)
    n_frames = len(times_samples)
    times = times_samples / fs
    
    # Create frame matrix
    frame_mat = np.zeros((frame_len, n_frames))
    for i, t in enumerate(times_samples):
        frame_mat[:, i] = residual[t - half_dur:t + half_dur]
    
    # Apply window (Blackman like MATLAB)
    window = np.blackman(frame_len)
    frame_mat = frame_mat * window[:, np.newaxis]
    
    # Mean subtraction
    frame_mat = frame_mat - np.mean(frame_mat, axis=0)
    
    # Compute spectrogram (matching MATLAB)
    spec_mat = np.zeros((fs // 2, n_frames))
    for i in range(n_frames):
        spec = np.abs(fft(frame_mat[:, i], fs))
        spec_mat[:, i] = spec[:fs // 2]
    
    # Normalize spectra
    spec_norms = np.sqrt(np.sum(spec_mat**2, axis=0))
    spec_norms[spec_norms == 0] = 1  # Avoid division by zero
    spec_mat = spec_mat / spec_norms
    
    # 2-iteration refinement (matching MATLAB)
    f0min = f0min_orig
    f0max = f0max_orig
    no_pitch_range_adjustment = True
    
    for iteration in range(n_iter):
        # Compute SRH for all frames
        f0, srh_values = _compute_srh_all_frames(spec_mat, n_harmonics, f0min, f0max, fs)
        
        # Refine F0 range if possible
        if np.max(srh_values) > srh_iter_thresh:
            # Get median of high-SRH frames
            high_srh_mask = srh_values > srh_iter_thresh
            if np.sum(high_srh_mask) > 0:
                f0_med_est = np.median(f0[high_srh_mask])
                
                # Refine lower bound
                f0min_refined = int(0.5 * f0_med_est)
                if f0min_refined > f0min:
                    f0min = f0min_refined
                    no_pitch_range_adjustment = False
                
                # Refine upper bound
                f0max_refined = int(2.0 * f0_med_est)
                if f0max_refined < f0max:
                    f0max = f0max_refined
                    no_pitch_range_adjustment = False
        
        # Break if no adjustment made
        if no_pitch_range_adjustment:
            break
    
    # Voice/unvoiced decision (matching MATLAB)
    vuv = np.zeros(n_frames, dtype=bool)
    
    # Adapt threshold based on SRH standard deviation
    voicing_threshold = voicing_thresh
    if np.std(srh_values) > srh_std_thresh:
        voicing_threshold = voicing_thresh2
    
    vuv = srh_values > voicing_threshold
    
    return f0, vuv, srh_values, times


def _compute_srh_all_frames(spec_mat, n_harmonics, f0min, f0max, fs):
    """
    Compute SRH for all frames in spectrogram (matching MATLAB exactly)
    
    SRH = Sum of harmonics - Sum of subharmonics
    
    Parameters
    ----------
    spec_mat : ndarray
        Spectrogram matrix (freq x time)
    n_harmonics : int
        Number of harmonics to sum
    f0min : float
        Minimum F0
    f0max : float
        Maximum F0
    fs : float
        Sampling frequency
        
    Returns
    -------
    f0 : ndarray
        F0 estimates for each frame
    srh_values : ndarray
        SRH values for each frame
    """
    n_frames = spec_mat.shape[1]
    n_freqs = spec_mat.shape[0]
    
    # F0 candidates (MATLAB uses integer Hz)
    f0_candidates = np.arange(int(f0min), int(f0max) + 1, dtype=int)
    n_candidates = len(f0_candidates)
    
    # Initialize SRH matrix
    srh_mat = np.zeros((int(f0max) + 1, n_frames))
    
    # Prepare harmonic indices (matching MATLAB logic)
    # plusIdx: harmonics (h * f0) for h = 1, 2, ..., nHarmonics
    # subtrIdx: subharmonics ((h+0.5) * f0) for h = 1, 2, ..., nHarmonics-1
    
    harmonics = np.arange(1, n_harmonics + 1)[:, np.newaxis]  # Column vector
    subharmonics = np.arange(1, n_harmonics) + 0.5  # 1.5, 2.5, ..., nHarmonics-0.5
    subharmonics = subharmonics[:, np.newaxis]
    
    # For each F0 candidate
    plus_idx = harmonics * f0_candidates  # nHarmonics x nCandidates
    subtr_idx = np.round(subharmonics * f0_candidates).astype(int)
    
    # Adjust indices to be within spectrum range (1-indexed in MATLAB, 0-indexed here)
    # MATLAB: mod(plusIdx-1, size(specMat,1)) + 1
    # Python: mod(plusIdx, n_freqs) with 0-indexing
    plus_idx = np.mod(plus_idx, n_freqs).astype(int)
    subtr_idx = np.mod(subtr_idx, n_freqs).astype(int)
    
    # Process each frame
    for frame_idx in range(n_frames):
        spectrum = spec_mat[:, frame_idx]
        
        # For each F0 candidate
        for i, f0_cand in enumerate(f0_candidates):
            # Sum harmonics
            harm_sum = np.sum(spectrum[plus_idx[:, i]])
            
            # Sum subharmonics  
            subharm_sum = np.sum(spectrum[subtr_idx[:, i]])
            
            # SRH = harmonics - subharmonics
            srh_mat[f0_cand, frame_idx] = harm_sum - subharm_sum
    
    # Extract maximum for each frame
    f0 = np.zeros(n_frames)
    srh_values = np.zeros(n_frames)
    
    for frame_idx in range(n_frames):
        # Find maximum SRH for this frame
        max_idx = np.argmax(srh_mat[int(f0min):int(f0max)+1, frame_idx])
        f0[frame_idx] = f0_candidates[max_idx]
        srh_values[frame_idx] = srh_mat[f0_candidates[max_idx], frame_idx]
    
    return f0, srh_values


def _compute_frame_srh(frame, periods, fs):
    """
    Compute SRH values for all period candidates in a frame
    
    Parameters
    ----------
    frame : ndarray
        Residual frame
    periods : ndarray
        Period candidates in samples
    fs : float
        Sampling frequency
        
    Returns
    -------
    srh : ndarray
        SRH values for each period candidate
    """
    n_periods = len(periods)
    srh = np.zeros(n_periods)
    
    # Compute spectrum
    n_fft = 2 ** int(np.ceil(np.log2(len(frame))))
    spectrum = np.abs(fft(frame, n_fft))
    
    # For each period candidate
    for i, period in enumerate(periods):
        f0 = fs / period
        
        # Sum harmonic energy
        harm_sum = 0.0
        n_harmonics = int(fs / (2 * f0))  # Up to Nyquist
        
        for h in range(1, min(n_harmonics + 1, 20)):  # Consider up to 20 harmonics
            freq_idx = int(h * f0 * n_fft / fs)
            
            if freq_idx < len(spectrum):
                # Sum energy in a small band around harmonic
                band_width = 2
                start_idx = max(0, freq_idx - band_width)
                end_idx = min(len(spectrum), freq_idx + band_width + 1)
                harm_sum += np.sum(spectrum[start_idx:end_idx] ** 2)
        
        srh[i] = harm_sum
    
    return srh

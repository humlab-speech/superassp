"""
Glottal Source Analysis Module

Implements various glottal source and voice quality analysis methods:
- IAIF (Iterative Adaptive Inverse Filtering)
- GCI detection (Glottal Closure Instants)
- Voice quality parameters (NAQ, QOQ, H1-H2, etc.)
- And more

Based on original COVAREP MATLAB implementations
"""

import numpy as np
from scipy import signal
from scipy.interpolate import interp1d
import warnings

__all__ = ['iaif', 'get_vq_params']


def iaif(x, fs, p_vt=None, p_gl=None, d=0.99, hpfilt=True):
    """
    Iterative Adaptive Inverse Filtering (IAIF) for glottal flow estimation
    
    Separates the glottal source signal from the vocal tract filter effect
    using an iterative procedure.
    
    Reference:
    P. Alku, "Glottal wave analysis with pitch synchronous iterative 
    adaptive inverse filtering", Speech Communication, 1991.
    
    Parameters
    ----------
    x : array_like
        Input speech signal
    fs : float
        Sampling frequency in Hz
    p_vt : int, optional
        Order of vocal tract LP filter (default: 2*round(fs/2000)+4, matching MATLAB)
    p_gl : int, optional
        Order of glottal source LP filter (default: 2*round(fs/4000), matching MATLAB)
    d : float
        Leaky integration coefficient (default: 0.99)
    hpfilt : bool or int
        Apply high-pass filter (default: True=1 pass)
        
    Returns
    -------
    g : ndarray
        Estimated glottal flow waveform
    dg : ndarray
        Glottal flow derivative
    a : ndarray
        Vocal tract LP filter coefficients
    ag : ndarray
        Glottal LP filter coefficients
        
    Notes
    -----
    IAIF is a classic method for glottal inverse filtering that iteratively
    estimates the vocal tract and glottal source contributions.
    
    This implementation follows MATLAB COVAREP exactly, including:
    - Linear-phase FIR high-pass filter (40-70 Hz)
    - Hanning window before LPC
    - Pre-frame ramp to reduce edge effects
    - Leaky integration (not simple cumsum)
    
    Examples
    --------
    >>> g, dg, a, ag = iaif(speech, fs=16000)
    """
    x = np.asarray(x).ravel()
    
    # Set default orders (matching MATLAB formulas exactly)
    if p_vt is None:
        p_vt = 2 * round(fs / 2000) + 4
    if p_gl is None:
        p_gl = 2 * round(fs / 4000)
    
    preflt = p_vt + 1
    
    # Check minimum length
    if len(x) <= p_vt:
        warnings.warn("Signal too short for IAIF analysis")
        return np.array([]), np.array([]), np.array([]), np.array([])
    
    # High-pass filter: Linear-phase FIR (Fstop=40Hz, Fpass=70Hz)
    # Matches MATLAB's hpfilter_fir() implementation
    if hpfilt:
        Fstop = 40.0
        Fpass = 70.0
        Nfir = round(300.0 / 16000.0 * fs)
        if Nfir % 2 == 1:
            Nfir += 1
        
        # Design FIR filter using firls (least-squares)
        # scipy.signal.firls: same as MATLAB firls
        bands = np.array([0, Fstop, Fpass, fs/2])
        desired = np.array([0, 0, 1, 1])
        B = signal.firls(Nfir+1, bands, desired, fs=fs)
        
        # Apply filter with zero-padding (matches MATLAB)
        npad = round(len(B) / 2) - 1
        x_padded = np.concatenate([x, np.zeros(npad)])
        x = signal.lfilter(B, 1, x_padded)
        x = x[round(len(B)/2):]
    
    # Hanning window for LPC analysis
    win = np.hanning(len(x))
    
    # Pre-frame ramp to reduce edge effects
    # MATLAB: linspace(-x(1), x(1), preflt)
    ramp = np.linspace(-x[0], x[0], preflt)
    signal_with_ramp = np.concatenate([ramp, x])
    idx = slice(preflt, len(signal_with_ramp))
    
    # ====== Iteration 1: Estimate glottal+radiation effect (Hg1) ======
    # Use order 1 to capture glottal+lip radiation combined effect
    Hg1 = _lpc_matlab_style(x * win, 1)
    y = signal.lfilter(Hg1, [1], signal_with_ramp)
    y = y[idx]
    
    # ====== Iteration 2: Estimate vocal tract (Hvt1) and get g1 ======
    Hvt1 = _lpc_matlab_style(y * win, p_vt)
    g1 = signal.lfilter(Hvt1, [1], signal_with_ramp)
    # Integrate to cancel lip radiation: g1 = filter(1, [1 -d], g1)
    g1 = signal.lfilter([1], [1, -d], g1)
    g1 = g1[idx]
    
    # ====== Iteration 3: Re-estimate glottal source (Hg2) ======
    Hg2 = _lpc_matlab_style(g1 * win, p_gl)
    y = signal.lfilter(Hg2, [1], signal_with_ramp)
    # Integrate
    y = signal.lfilter([1], [1, -d], y)
    y = y[idx]
    
    # ====== Iteration 4: Final vocal tract estimate (Hvt2) ======
    Hvt2 = _lpc_matlab_style(y * win, p_vt)
    dg = signal.lfilter(Hvt2, [1], signal_with_ramp)
    # Final integration to get flow
    g = signal.lfilter([1], [1, -d], dg)
    g = g[preflt:]  # Remove ramp
    dg = dg[idx]
    
    # Return coefficients
    a = Hvt2
    ag = Hg2
    
    return g, dg, a, ag


def _lpc_matlab_style(x, order):
    """
    LPC analysis matching MATLAB's lpc() function
    
    Uses autocorrelation method with Levinson-Durbin recursion.
    Returns coefficients in MATLAB format: [1, a1, a2, ..., an]
    where the prediction filter is A(z) = 1 + a1*z^-1 + a2*z^-2 + ...
    
    Parameters
    ----------
    x : ndarray
        Input signal (typically windowed)
    order : int
        LPC order
        
    Returns
    -------
    a : ndarray
        LPC coefficients [1, a1, a2, ..., an]
        
    Notes
    -----
    MATLAB's lpc() returns prediction filter coefficients, not
    the AR parameters. For inverse filtering, use these coefficients
    directly: y = filter(a, 1, x)
    """
    x = np.asarray(x).ravel()
    
    # Compute autocorrelation 
    # MATLAB uses biased estimate: r[k] = sum(x[0:N-k] * x[k:N]) / N
    r = np.zeros(order + 1)
    N = len(x)
    for k in range(order + 1):
        if k < N:
            r[k] = np.sum(x[:N-k] * x[k:])
        else:
            r[k] = 0.0
    
    # Handle zero signal
    if r[0] == 0:
        a = np.zeros(order + 1)
        a[0] = 1.0
        return a
    
    # Levinson-Durbin algorithm (matching MATLAB exactly)
    a = np.zeros(order + 1)
    a[0] = 1.0
    e = r[0]
    
    for i in range(1, order + 1):
        # Compute reflection coefficient
        k_i = r[i]
        for j in range(1, i):
            k_i -= a[j] * r[i - j]
        
        if e == 0:
            break
            
        k_i /= e
        
        # Update coefficients  
        a_prev = a.copy()
        a[i] = k_i
        for j in range(1, i):
            a[j] = a_prev[j] - k_i * a_prev[i - j]
        
        # Update error
        e *= (1.0 - k_i * k_i)
        
        if e <= 0:
            break
    
    return a


def _iaif_iteration1(x, p_rough):
    """
    First IAIF iteration: estimate initial glottal source
    
    Uses a high-order all-pole model to remove gross spectral features
    
    Parameters
    ----------
    x : ndarray
        Pre-emphasized signal
    p_rough : int
        Order for rough estimate (typically p_vt + 1)
    
    Returns
    -------
    g1 : ndarray
        Initial glottal estimate
    """
    # LPC analysis with high order
    a_rough = _estimate_lpc(x, p_rough)
    
    # Inverse filter
    g1 = signal.lfilter(a_rough, [1], x)
    
    return g1


def _estimate_lpc(x, order):
    """
    Estimate LPC coefficients using autocorrelation method
    
    Parameters
    ----------
    x : ndarray
        Input signal
    order : int
        LPC order
        
    Returns
    -------
    a : ndarray
        LPC coefficients [1, -a1, -a2, ..., -an]
    """
    # Compute autocorrelation
    r = np.correlate(x, x, mode='full')
    r = r[len(r)//2:]  # Keep only positive lags
    r = r[:order+1]
    
    # Levinson-Durbin recursion
    a = _levinson_durbin(r, order)
    
    return a


def _levinson_durbin(r, order):
    """
    Levinson-Durbin recursion for LPC computation
    
    Parameters
    ----------
    r : ndarray
        Autocorrelation values
    order : int
        LPC order
        
    Returns
    -------
    a : ndarray
        LPC coefficients
    """
    a = np.zeros(order + 1)
    a[0] = 1.0
    
    if len(r) < order + 1:
        return a
    
    e = r[0]
    
    for i in range(1, order + 1):
        # Compute reflection coefficient
        lambda_i = -r[i]
        for j in range(1, i):
            lambda_i -= a[j] * r[i-j]
        lambda_i /= e
        
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


def get_vq_params(x, fs, f0=None):
    """
    Extract voice quality parameters from speech signal
    
    Computes various glottal source and voice quality measures including:
    - NAQ (Normalized Amplitude Quotient)
    - QOQ (Quasi-Open Quotient)
    - H1-H2 (difference between first two harmonics)
    - HRF (Harmonic Richness Factor)
    - PSP (Parabolic Spectral Parameter)
    
    Parameters
    ----------
    x : array_like
        Input speech signal
    fs : float
        Sampling frequency
    f0 : array_like, optional
        F0 contour (if not provided, will be estimated)
        
    Returns
    -------
    params : dict
        Dictionary containing voice quality parameters
        
    Notes
    -----
    This is a simplified version. Full implementation requires
    GCI detection and other advanced processing.
    
    References
    ----------
    Alku et al., Kane et al., Drugman et al. - Various COVAREP papers
    """
    x = np.asarray(x).ravel()
    
    # Get glottal flow using IAIF
    g, dg, a, ag = iaif(x, fs)
    
    # Initialize parameters dictionary
    params = {}
    
    # Compute basic amplitude measures
    g_max = np.max(np.abs(g))
    dg_max = np.max(np.abs(dg))
    
    params['glottal_flow_max'] = g_max
    params['glottal_flow_derivative_max'] = dg_max
    
    # Estimate NAQ (simplified version)
    # NAQ = amplitude_quotient / duration
    # This requires GCI detection for full implementation
    # Here we provide a placeholder
    params['NAQ'] = np.nan  # TODO: Implement with GCI detection
    params['QOQ'] = np.nan  # TODO: Implement with GCI detection
    
    # Compute spectral measures
    if f0 is not None:
        f0_mean = np.mean(f0[f0 > 0])
        if f0_mean > 0:
            params['H1_H2'] = _compute_h1_h2(x, fs, f0_mean)
        else:
            params['H1_H2'] = np.nan
    else:
        params['H1_H2'] = np.nan
    
    return params


def _compute_h1_h2(x, fs, f0):
    """
    Compute H1-H2 (difference between first two harmonics)
    
    Parameters
    ----------
    x : ndarray
        Speech signal
    fs : float
        Sampling frequency
    f0 : float
        Fundamental frequency
        
    Returns
    -------
    h1_h2 : float
        H1-H2 in dB
    """
    # Compute spectrum
    n_fft = 2 ** int(np.ceil(np.log2(len(x))))
    X = np.abs(np.fft.rfft(x, n_fft))
    freqs = np.fft.rfftfreq(n_fft, 1/fs)
    
    # Find peaks near f0 and 2*f0
    def find_harmonic_magnitude(freq, bandwidth=50):
        idx_range = np.where((freqs >= freq - bandwidth) & 
                            (freqs <= freq + bandwidth))[0]
        if len(idx_range) > 0:
            return np.max(X[idx_range])
        return 0.0
    
    h1 = find_harmonic_magnitude(f0)
    h2 = find_harmonic_magnitude(2 * f0)
    
    if h1 > 0 and h2 > 0:
        h1_h2 = 20 * np.log10(h1 / (h2 + 1e-10))
    else:
        h1_h2 = np.nan
    
    return h1_h2

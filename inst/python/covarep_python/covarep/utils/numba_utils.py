"""
Numba JIT-compiled utilities for performance-critical functions

Numba provides near-C performance without requiring compilation at install time.
These functions are automatically JIT-compiled on first use.

Advantages:
- No build step required (unlike Cython)
- Near-C performance (typically 80-95% of Cython speed)
- Easy to maintain (pure Python syntax with type hints)

Disadvantages:
- First call is slower (JIT compilation)
- Limited NumPy function support in nopython mode
"""

import numpy as np

try:
    import numba as nb
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    # Create dummy decorator if numba not available
    class nb:
        @staticmethod
        def jit(*args, **kwargs):
            def decorator(func):
                return func
            return decorator


__all__ = [
    'NUMBA_AVAILABLE',
    'levinson_durbin_numba',
    'compute_autocorrelation_numba',
    'lpc_analysis_numba',
]


@nb.jit(nopython=True, cache=True, fastmath=True)
def levinson_durbin_numba(r, order):
    """
    Numba JIT-compiled Levinson-Durbin recursion

    Expected speedup: 5-8x over pure Python/NumPy
    Similar performance to Cython, but no compilation required.

    Parameters
    ----------
    r : ndarray (float64)
        Autocorrelation sequence
    order : int
        LPC order

    Returns
    -------
    a : ndarray (float64)
        LPC coefficients [1, a1, a2, ..., an]
    """
    a = np.zeros(order + 1, dtype=np.float64)
    a[0] = 1.0

    if len(r) < order + 1:
        return a

    if r[0] == 0:
        return a

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
        # Note: In nopython mode, we need explicit loops
        a_prev = a.copy()
        a[i] = k_i
        for j in range(1, i):
            a[j] = a_prev[j] - k_i * a_prev[i - j]

        # Update error
        e *= (1.0 - k_i * k_i)

        if e <= 0:
            break

    return a


@nb.jit(nopython=True, cache=True, fastmath=True)
def compute_autocorrelation_numba(x, max_lag):
    """
    Compute autocorrelation using direct method (Numba-optimized)

    Faster than np.correlate for small signals and low lag orders.

    Parameters
    ----------
    x : ndarray (float64)
        Input signal
    max_lag : int
        Maximum lag (typically LPC order)

    Returns
    -------
    r : ndarray (float64)
        Autocorrelation values [r(0), r(1), ..., r(max_lag)]
    """
    N = len(x)
    r = np.zeros(max_lag + 1, dtype=np.float64)

    for k in range(max_lag + 1):
        if k < N:
            acc = 0.0
            for n in range(N - k):
                acc += x[n] * x[n + k]
            r[k] = acc
        else:
            r[k] = 0.0

    return r


@nb.jit(nopython=True, cache=True)
def lpc_analysis_numba(frame, order):
    """
    Complete LPC analysis using Numba (autocorrelation + Levinson-Durbin)

    Combines autocorrelation and Levinson-Durbin for efficient LPC analysis.

    Parameters
    ----------
    frame : ndarray (float64)
        Windowed speech frame
    order : int
        LPC order

    Returns
    -------
    a : ndarray (float64)
        LPC coefficients [1, a1, a2, ..., an]
    """
    # Compute autocorrelation
    r = compute_autocorrelation_numba(frame, order)

    # Levinson-Durbin
    a = levinson_durbin_numba(r, order)

    return a


@nb.jit(nopython=True, cache=True, fastmath=True, parallel=False)
def compute_srh_inner_loop_numba(spec_mat, plus_idx, subtr_idx, n_candidates, n_frames):
    """
    Inner loop of SRH computation (Numba-optimized)

    This is the most performance-critical part of F0 tracking.

    Parameters
    ----------
    spec_mat : ndarray (float64, 2D)
        Magnitude spectrogram (n_freqs x n_frames)
    plus_idx : ndarray (int32, 2D)
        Harmonic indices (n_harmonics x n_candidates)
    subtr_idx : ndarray (int32, 2D)
        Subharmonic indices (n_harmonics-1 x n_candidates)
    n_candidates : int
        Number of F0 candidates
    n_frames : int
        Number of frames

    Returns
    -------
    srh_mat : ndarray (float64, 2D)
        SRH values (n_candidates x n_frames)
    """
    n_harmonics = plus_idx.shape[0]
    n_subharmonics = subtr_idx.shape[0]

    srh_mat = np.zeros((n_candidates, n_frames), dtype=np.float64)

    # Loop over candidates
    for cand_idx in range(n_candidates):
        # Loop over frames
        for frame_idx in range(n_frames):
            # Sum harmonics
            harm_sum = 0.0
            for h in range(n_harmonics):
                freq_idx = plus_idx[h, cand_idx]
                harm_sum += spec_mat[freq_idx, frame_idx]

            # Sum subharmonics
            subharm_sum = 0.0
            for h in range(n_subharmonics):
                freq_idx = subtr_idx[h, cand_idx]
                subharm_sum += spec_mat[freq_idx, frame_idx]

            # SRH = harmonics - subharmonics
            srh_mat[cand_idx, frame_idx] = harm_sum - subharm_sum

    return srh_mat


@nb.jit(nopython=True, cache=True)
def find_peaks_1d_numba(x, threshold=0.0):
    """
    Simple 1D peak finding (Numba-optimized)

    Finds local maxima above threshold.

    Parameters
    ----------
    x : ndarray (float64)
        Input signal
    threshold : float
        Minimum peak value

    Returns
    -------
    peaks : ndarray (int64)
        Indices of peaks
    """
    n = len(x)
    peaks_list = []

    for i in range(1, n - 1):
        if x[i] > x[i - 1] and x[i] > x[i + 1] and x[i] > threshold:
            peaks_list.append(i)

    # Convert list to array
    peaks = np.array(peaks_list, dtype=np.int64)
    return peaks


@nb.jit(nopython=True, cache=True, fastmath=True)
def interpolate_linear_numba(x, y, x_new):
    """
    Linear interpolation (Numba-optimized)

    Simple linear interpolation for 1D data.

    Parameters
    ----------
    x : ndarray (float64)
        Original x values (must be sorted)
    y : ndarray (float64)
        Original y values
    x_new : ndarray (float64)
        New x values to interpolate

    Returns
    -------
    y_new : ndarray (float64)
        Interpolated y values
    """
    n = len(x)
    m = len(x_new)
    y_new = np.zeros(m, dtype=np.float64)

    for i in range(m):
        xi = x_new[i]

        # Find bracketing indices
        if xi <= x[0]:
            y_new[i] = y[0]
        elif xi >= x[n - 1]:
            y_new[i] = y[n - 1]
        else:
            # Binary search would be faster, but linear search is simple
            for j in range(n - 1):
                if x[j] <= xi <= x[j + 1]:
                    # Linear interpolation
                    alpha = (xi - x[j]) / (x[j + 1] - x[j])
                    y_new[i] = (1 - alpha) * y[j] + alpha * y[j + 1]
                    break

    return y_new


@nb.jit(nopython=True, cache=True)
def median_filter_1d_numba(x, window_size=3):
    """
    1D median filter (Numba-optimized)

    Parameters
    ----------
    x : ndarray (float64)
        Input signal
    window_size : int
        Window size (must be odd)

    Returns
    -------
    y : ndarray (float64)
        Filtered signal
    """
    n = len(x)
    y = np.zeros(n, dtype=np.float64)
    half_win = window_size // 2

    for i in range(n):
        start = max(0, i - half_win)
        end = min(n, i + half_win + 1)

        # Extract window
        window = x[start:end]

        # Compute median (simple sorting)
        window_sorted = np.sort(window)
        median_idx = len(window_sorted) // 2
        y[i] = window_sorted[median_idx]

    return y


# Utility function to check if Numba is working
def test_numba_performance():
    """
    Test Numba performance with a simple benchmark

    Returns
    -------
    speedup : float
        Speedup factor (Numba vs NumPy)
    """
    if not NUMBA_AVAILABLE:
        print("⚠ Numba not available")
        return None

    import time

    # Generate test data
    n = 10000
    order = 20
    x = np.random.randn(n)

    # Compute autocorrelation (NumPy)
    start = time.perf_counter()
    for _ in range(100):
        r_np = np.correlate(x, x, mode='full')[n-1:n+order]
    time_numpy = time.perf_counter() - start

    # Compute autocorrelation (Numba)
    # Warmup
    _ = compute_autocorrelation_numba(x, order)

    start = time.perf_counter()
    for _ in range(100):
        r_numba = compute_autocorrelation_numba(x, order)
    time_numba = time.perf_counter() - start

    speedup = time_numpy / time_numba

    print(f"Numba Performance Test:")
    print(f"  NumPy:  {time_numpy*1000:.2f} ms")
    print(f"  Numba:  {time_numba*1000:.2f} ms")
    print(f"  Speedup: {speedup:.2f}x")

    return speedup


if __name__ == "__main__":
    print("=" * 60)
    print("Numba Utilities Test")
    print("=" * 60)

    if not NUMBA_AVAILABLE:
        print("\n✗ Numba is not installed")
        print("  Install with: pip install numba")
    else:
        print(f"\n✓ Numba version: {nb.__version__}")
        print("\nRunning performance test...")
        speedup = test_numba_performance()

        if speedup and speedup > 2.0:
            print(f"\n✓ Numba is working well! ({speedup:.1f}x speedup)")
        elif speedup:
            print(f"\n⚠ Numba speedup is modest ({speedup:.1f}x)")
        else:
            print("\n✗ Numba test failed")

    print("\n" + "=" * 60)

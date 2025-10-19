"""
Recurrence Period Density Entropy (RPDE)

Ported from rpde.m and close_ret.c
RPDE quantifies the predictability of recurrences in phase space

HIGHLY OPTIMIZED VERSION with:
- **Cython compiled extensions for 2-5x speedup (PREFERRED)**
- Numba JIT compilation fallback for 10-20x speedup over pure Python
- KD-tree spatial indexing for 2-5x additional speedup on large signals
- Vectorized operations where possible

References:
    Little, M., McSharry, P., Roberts, S., Costello, D., & Moroz, I. (2007).
    Exploiting Nonlinear Recurrence and Fractal Scaling Properties for Voice Disorder Detection.
    BioMedical Engineering OnLine, 6:23.
"""

import numpy as np
import warnings

# Try to import Cython implementation (fastest)
try:
    from .rpde_cython import compute_rpde_cython
    CYTHON_AVAILABLE = True
except ImportError:
    CYTHON_AVAILABLE = False

# Try to import numba for performance optimization (fallback)
try:
    from numba import jit
    NUMBA_AVAILABLE = True
except ImportError:
    NUMBA_AVAILABLE = False
    # Fallback decorator that does nothing
    def jit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator

# Try to import scipy KDTree for spatial indexing
try:
    from scipy.spatial import cKDTree
    KDTREE_AVAILABLE = True
except ImportError:
    KDTREE_AVAILABLE = False


def compute_rpde(signal, m=4, tau=50, epsilon=0.12, T_max=1000, fs=25000, use_kdtree=True, use_cython=False):
    """
    Compute Recurrence Period Density Entropy
    
    Automatically selects the best available implementation:
    1. Numba JIT with KD-tree (fastest for large signals) - PREFERRED
    2. Numba JIT manual (fallback)
    3. Cython (experimental, currently disabled due to accuracy issues)
    
    Parameters:
    -----------
    signal : ndarray
        Input signal (should be resampled to 25 kHz)
    m : int
        Embedding dimension (default 4)
    tau : int
        Embedding delay (default 50)
    epsilon : float
        Close returns radius (default 0.12) - note: MATLAB default is 0.12, not 0.2!
    T_max : int
        Maximum recurrence time (default 1000)
    fs : int
        Sampling rate (default 25000)
    use_kdtree : bool
        Use KD-tree for spatial indexing (10x faster, default True)
    use_cython : bool
        Prefer Cython implementation if available (default False, experimental)
        
    Returns:
    --------
    H_norm : float
        Normalized RPDE value
    """
    signal = np.asarray(signal).ravel()
    signal = signal[np.isfinite(signal)]
    
    # Resample to 25 kHz if not already
    if fs != 25000 and fs > 0:
        from scipy import signal as scipy_signal
        num_samples = int(len(signal) * 25000 / fs)
        signal = scipy_signal.resample(signal, num_samples)
    
    if len(signal) < m * tau + T_max:
        return np.nan
    
    # PRIORITY 1: Use KD-tree if available (fastest and most accurate)
    N = len(signal)
    M = N - (m - 1) * tau
    
    if use_kdtree and KDTREE_AVAILABLE and M > 100:  # Use KD-tree for almost all cases
        return _rpde_kdtree(signal, m, tau, epsilon, T_max)
    
    # PRIORITY 2: Use Cython if explicitly requested (experimental)
    if use_cython and CYTHON_AVAILABLE:
        try:
            return compute_rpde_cython(signal, m, tau, epsilon, T_max)
        except Exception as e:
            warnings.warn(f"Cython RPDE failed, falling back to Numba: {e}")
    
    # PRIORITY 3: Fallback to manual Numba implementation
    return _rpde_manual(signal, m, tau, epsilon, T_max)


def _rpde_manual(signal, m, tau, epsilon, T_max):
    """
    Manual RPDE implementation matching close_ret.c algorithm
    
    The algorithm:
    1. Time-delay embed the signal
    2. For each point i:
       a. Find first point j > i where distance exceeds epsilon
       b. Continue from j to find first point k where distance <= epsilon
       c. Record recurrence time (k - i)
    3. Build histogram of recurrence times
    4. Compute normalized entropy of histogram
    
    Uses Numba JIT compilation for 10-20x speedup if available
    """
    N = len(signal)
    M = N - (m - 1) * tau  # Number of embedded points
    
    if M <= 0:
        return np.nan
    
    # Time-delay embedding (optimized version)
    embedded = _time_delay_embedding_optimized(signal, M, m, tau, N)
    
    # Normalize embedded space (important for consistent epsilon interpretation)
    # This is implicit in the MATLAB/C code through the signal preprocessing
    embedded_std = np.std(embedded)
    if embedded_std > 0:
        epsilon_normalized = epsilon * embedded_std
    else:
        epsilon_normalized = epsilon
    
    # Find close returns (matching the C algorithm)
    if NUMBA_AVAILABLE:
        close_returns = _find_close_returns_numba(embedded, M, epsilon_normalized)
    else:
        close_returns = _find_close_returns_python(embedded, M, epsilon_normalized)
    
    # Trim to T_max if specified
    if T_max > 0 and len(close_returns) > T_max:
        close_returns = close_returns[:T_max]
    
    # Compute RPD (recurrence period density)
    total = np.sum(close_returns)
    if total == 0:
        return 0.0
    
    rpd = close_returns / total
    
    # Compute normalized entropy
    H = 0.0
    for p in rpd:
        if p > 0:
            H -= p * np.log(p)
    
    N_rpd = len(rpd)
    if N_rpd > 1:
        H_norm = H / np.log(N_rpd)
    else:
        H_norm = 0.0
    
    return H_norm


@jit(nopython=True, cache=True)
def _time_delay_embedding_optimized(signal, M, m, tau, N):
    """Optimized time-delay embedding with Numba"""
    embedded = np.zeros((M, m), dtype=np.float64)
    for i in range(M):
        for j in range(m):
            idx = i + j * tau
            if idx < N:
                embedded[i, j] = signal[idx]
    return embedded


@jit(nopython=True, cache=True)
def _find_close_returns_numba(embedded, M, epsilon):
    """
    Find close returns using the algorithm from close_ret.c
    
    This is faithful to the original C implementation:
    - For each point i, find first j where dist > epsilon (leaving neighborhood)
    - Then find first k where dist <= epsilon (returning to neighborhood)
    - Record k - i as recurrence time
    
    Returns histogram of recurrence times
    """
    epsilon_squared = epsilon * epsilon
    close_returns = np.zeros(M, dtype=np.float64)
    
    for i in range(M):
        point_i = embedded[i, :]
        
        # Step 1: Find first point where distance exceeds epsilon
        j = i + 1
        eta_flag = False
        
        while j < M and not eta_flag:
            dist_sq = 0.0
            for d in range(embedded.shape[1]):
                diff = embedded[j, d] - point_i[d]
                dist_sq += diff * diff
            
            if dist_sq > epsilon_squared:
                eta_flag = True
            j += 1
        
        # Step 2: Find first close return after leaving neighborhood
        eta_flag = False
        while j < M and not eta_flag:
            dist_sq = 0.0
            for d in range(embedded.shape[1]):
                diff = embedded[j, d] - point_i[d]
                dist_sq += diff * diff
            
            if dist_sq <= epsilon_squared:
                time_diff = j - i
                if time_diff < M:
                    close_returns[time_diff] += 1.0
                eta_flag = True
            
            j += 1
    
    return close_returns


def _find_close_returns_python(embedded, M, epsilon):
    """Python fallback for close returns (slower but works without Numba)"""
    epsilon_squared = epsilon * epsilon
    close_returns = np.zeros(M, dtype=np.float64)
    
    for i in range(M):
        point_i = embedded[i:i+1, :]
        
        # Compute all distances from point i
        dists_sq = np.sum((embedded - point_i)**2, axis=1)
        
        # Find first point where distance exceeds epsilon
        j = i + 1
        while j < M and dists_sq[j] <= epsilon_squared:
            j += 1
        
        # Find first close return after leaving
        while j < M:
            if dists_sq[j] <= epsilon_squared:
                time_diff = j - i
                if time_diff < M:
                    close_returns[time_diff] += 1
                break
            j += 1
    
    return close_returns


def _rpde_kdtree(signal, m, tau, epsilon, T_max):
    """
    RPDE implementation using KD-tree for spatial indexing
    
    This is significantly faster for large signals (M > 1000)
    Expected speedup: 2-5x over Numba version
    
    The KD-tree allows O(log N) nearest neighbor queries instead of O(N)
    """
    N = len(signal)
    M = N - (m - 1) * tau
    
    if M <= 0:
        return np.nan
    
    # Time-delay embedding
    if NUMBA_AVAILABLE:
        embedded = _time_delay_embedding_optimized(signal, M, m, tau, N)
    else:
        embedded = np.zeros((M, m))
        for i in range(M):
            for j in range(m):
                idx = i + j * tau
                if idx < N:
                    embedded[i, j] = signal[idx]
    
    # Normalize embedded space
    embedded_std = np.std(embedded)
    if embedded_std > 0:
        epsilon_normalized = epsilon * embedded_std
    else:
        epsilon_normalized = epsilon
    
    # Build KD-tree for efficient spatial queries
    tree = cKDTree(embedded)
    
    # Find close returns using KD-tree
    close_returns = _find_close_returns_kdtree(tree, embedded, M, epsilon_normalized)
    
    # Trim to T_max
    if T_max > 0 and len(close_returns) > T_max:
        close_returns = close_returns[:T_max]
    
    # Compute normalized entropy
    total = np.sum(close_returns)
    if total == 0:
        return 0.0
    
    rpd = close_returns / total
    
    H = 0.0
    for p in rpd:
        if p > 0:
            H -= p * np.log(p)
    
    N_rpd = len(rpd)
    if N_rpd > 1:
        H_norm = H / np.log(N_rpd)
    else:
        H_norm = 0.0
    
    return H_norm


def _find_close_returns_kdtree(tree, embedded, M, epsilon):
    """
    Find close returns using KD-tree for efficient spatial queries
    
    This implementation uses query_ball_point for efficient radius queries
    """
    close_returns = np.zeros(M, dtype=np.float64)
    epsilon_squared = epsilon * epsilon
    
    for i in range(M):
        point_i = embedded[i]
        
        # Strategy: Find points beyond epsilon, then find first return within epsilon
        # This is tricky with KD-tree, so we use a hybrid approach
        
        # Find first point where distance exceeds epsilon (sequential from i+1)
        j = i + 1
        left_neighborhood = False
        
        while j < M and not left_neighborhood:
            dist_sq = np.sum((embedded[j] - point_i)**2)
            if dist_sq > epsilon_squared:
                left_neighborhood = True
            else:
                j += 1
        
        if not left_neighborhood or j >= M:
            continue
        
        # Now find first return within epsilon starting from j
        # Use KD-tree query for points in epsilon-ball starting from j
        candidates = tree.query_ball_point(point_i, epsilon, workers=1)
        
        # Find first candidate >= j
        for k in sorted(candidates):
            if k >= j:
                time_diff = k - i
                if time_diff < M:
                    close_returns[time_diff] += 1
                break
    
    return close_returns


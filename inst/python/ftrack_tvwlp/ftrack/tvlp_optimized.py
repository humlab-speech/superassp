"""
Optimized Time-Varying Linear Prediction (TVLP) algorithms.

Phase 1 optimizations: Vectorized matrix construction
"""

import numpy as np
from scipy.optimize import linprog


def tvlp_l2_opt(x, p, q):
    """
    Optimized TVLP L2 with vectorized matrix construction.

    Parameters
    ----------
    x : ndarray
        Input signal
    p : int
        LP order
    q : int
        Polynomial order

    Returns
    -------
    aki : ndarray
        Time-varying LP coefficients [q+1, p]
    """
    x = np.asarray(x, dtype=np.float64)
    n_samples = len(x)

    # Vectorized matrix construction
    n_equations = n_samples - p
    n_vals = np.arange(p, n_samples)

    # Build time power matrix: (n-p)^j for j=0 to q
    time_powers = (n_vals[:, None] - p) ** np.arange(q + 1)  # [N x (q+1)]

    # Build delayed signal matrix
    X_delayed = np.column_stack([x[p-i:n_samples-i] for i in range(1, p+1)])  # [N x p]

    # Kronecker-like product structure
    # For each delay i and polynomial order j: (n-p)^j * x[n-i]
    Ypu = np.zeros((n_equations, p * (q + 1)), dtype=np.float64)
    for i in range(p):
        cols_start = i * (q + 1)
        cols_end = (i + 1) * (q + 1)
        Ypu[:, cols_start:cols_end] = X_delayed[:, i:i+1] * time_powers

    Yn = x[p:n_samples]

    # Solve least squares
    x_l2, residuals, rank, s = np.linalg.lstsq(Ypu, Yn, rcond=None)

    # Reshape to [q+1, p]
    aki = x_l2.reshape(q + 1, p, order='F')

    return aki


def tvwlp_l2_opt(x, p, q, w):
    """
    Optimized TVWLP L2 with vectorized matrix construction.

    Parameters
    ----------
    x : ndarray
        Input signal
    p : int
        LP order
    q : int
        Polynomial order
    w : ndarray
        Weight function

    Returns
    -------
    aki : ndarray
        Time-varying LP coefficients [q+1, p]
    """
    x = np.asarray(x, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    n_samples = len(x)

    # Vectorized matrix construction
    n_equations = n_samples - p
    n_vals = np.arange(p, n_samples)

    # Build time power matrix
    time_powers = (n_vals[:, None] - p) ** np.arange(q + 1)  # [N x (q+1)]

    # Build delayed signal matrix
    X_delayed = np.column_stack([x[p-i:n_samples-i] for i in range(1, p+1)])  # [N x p]

    # Apply weights and build system
    w_active = w[p:n_samples]
    Ypu = np.zeros((n_equations, p * (q + 1)), dtype=np.float64)

    for i in range(p):
        cols_start = i * (q + 1)
        cols_end = (i + 1) * (q + 1)
        # Weight each term: w[n] * (n-p)^j * x[n-i]
        Ypu[:, cols_start:cols_end] = (w_active[:, None] * X_delayed[:, i:i+1]) * time_powers

    Yn = w_active * x[p:n_samples]

    # Solve weighted least squares
    x_l2, residuals, rank, s = np.linalg.lstsq(Ypu, Yn, rcond=None)

    # Reshape to [q+1, p]
    aki = x_l2.reshape(q + 1, p, order='F')

    return aki


def tvlp_l1_opt(x, p, q):
    """
    Optimized TVLP L1 with vectorized matrix construction.

    Parameters
    ----------
    x : ndarray
        Input signal
    p : int
        LP order
    q : int
        Polynomial order

    Returns
    -------
    aki : ndarray
        Time-varying LP coefficients [q+1, p]
    """
    x = np.asarray(x, dtype=np.float64)
    n_samples = len(x)

    # Vectorized matrix construction (same as tvlp_l2_opt)
    n_equations = n_samples - p
    n_vals = np.arange(p, n_samples)

    time_powers = (n_vals[:, None] - p) ** np.arange(q + 1)
    X_delayed = np.column_stack([x[p-i:n_samples-i] for i in range(1, p+1)])

    Ypu = np.zeros((n_equations, p * (q + 1)), dtype=np.float64)
    for i in range(p):
        cols_start = i * (q + 1)
        cols_end = (i + 1) * (q + 1)
        Ypu[:, cols_start:cols_end] = X_delayed[:, i:i+1] * time_powers

    Yn = x[p:n_samples]

    # L1 minimization via linear programming
    A = Ypu
    b = Yn
    m, n = A.shape

    c = np.concatenate([np.zeros(n), np.ones(m), np.ones(m)])
    A_eq = np.hstack([A, -np.eye(m), np.eye(m)])
    b_eq = b
    bounds = [(None, None)] * n + [(0, None)] * (2 * m)

    result = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds,
                     method='highs', options={'maxiter': 100, 'disp': False})

    if not result.success:
        print(f"Warning: linprog did not converge. Status: {result.status}")

    x_l1 = result.x[:n]
    aki = x_l1.reshape(q + 1, p, order='F')

    return aki


def tvwlp_l1_opt(x, p, q, w):
    """
    Optimized TVWLP L1 with vectorized matrix construction.

    Parameters
    ----------
    x : ndarray
        Input signal
    p : int
        LP order
    q : int
        Polynomial order
    w : ndarray
        Weight function

    Returns
    -------
    aki : ndarray
        Time-varying LP coefficients [q+1, p]
    """
    x = np.asarray(x, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    n_samples = len(x)

    # Vectorized matrix construction (same as tvwlp_l2_opt)
    n_equations = n_samples - p
    n_vals = np.arange(p, n_samples)

    time_powers = (n_vals[:, None] - p) ** np.arange(q + 1)
    X_delayed = np.column_stack([x[p-i:n_samples-i] for i in range(1, p+1)])

    w_active = w[p:n_samples]
    Ypu = np.zeros((n_equations, p * (q + 1)), dtype=np.float64)

    for i in range(p):
        cols_start = i * (q + 1)
        cols_end = (i + 1) * (q + 1)
        Ypu[:, cols_start:cols_end] = (w_active[:, None] * X_delayed[:, i:i+1]) * time_powers

    Yn = w_active * x[p:n_samples]

    # L1 minimization
    A = Ypu
    b = Yn
    m, n = A.shape

    c = np.concatenate([np.zeros(n), np.ones(m), np.ones(m)])
    A_eq = np.hstack([A, -np.eye(m), np.eye(m)])
    b_eq = b
    bounds = [(None, None)] * n + [(0, None)] * (2 * m)

    result = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds,
                     method='highs', options={'maxiter': 100, 'disp': False})

    if not result.success:
        print(f"Warning: linprog did not converge. Status: {result.status}")

    x_l1 = result.x[:n]
    aki = x_l1.reshape(q + 1, p, order='F')

    return aki

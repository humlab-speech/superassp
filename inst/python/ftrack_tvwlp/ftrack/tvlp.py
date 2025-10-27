"""
Time-Varying Linear Prediction (TVLP) and Time-Varying Weighted Linear Prediction (TVWLP)

This module implements L1 and L2 norm estimation methods for time-varying LP coefficients.
"""

import numpy as np
from scipy.optimize import linprog


def tvlp_l2(x, p, q):
    """
    Time-Varying Linear Prediction with L2 norm (least squares).

    Parameters
    ----------
    x : ndarray
        Input signal
    p : int
        LP order
    q : int
        Polynomial order for time-varying coefficients

    Returns
    -------
    aki : ndarray
        Time-varying LP coefficients of shape [q+1, p]
        Each column k represents coefficients for a_k(n) = sum(aki[i,k] * n^i)
    """
    x = np.asarray(x, dtype=np.float64)
    n_samples = len(x)

    # Build the system of equations
    # MATLAB: for n=p+1:length(x) means n from p+1 to length(x) inclusive
    # Python: range(p, n_samples) gives p to n_samples-1, which maps to MATLAB p+1 to length(x)
    n_equations = n_samples - p
    n_coeffs = p * (q + 1)

    Ypu = np.zeros((n_equations, n_coeffs), dtype=np.float64)
    Yn = np.zeros(n_equations, dtype=np.float64)

    m = 0
    for n in range(p, n_samples):  # MATLAB: p+1 to length(x)
        Yn[m] = x[n]

        for i in range(1, p + 1):  # MATLAB: 1 to p
            for j in range(q + 1):  # MATLAB: 0 to q
                # MATLAB indexing: (i-1)*(q+1)+j+1
                # Python indexing: (i-1)*(q+1)+j
                col_idx = (i - 1) * (q + 1) + j
                # MATLAB: (n-p-1)^j * x(n-i)
                # In Python: n corresponds to MATLAB n, so (n-p)^j * x[n-i]
                # Since MATLAB n goes from p+1 to length(x), and our n goes from p to n_samples-1
                # MATLAB's (n-p-1) equals our (n-p) when we account for the offset
                Ypu[m, col_idx] = ((n - p) ** j) * x[n - i]
        m += 1

    # Solve least squares: Ypu * x_l2 = Yn
    x_l2, residuals, rank, s = np.linalg.lstsq(Ypu, Yn, rcond=None)

    # Reshape to [q+1, p]
    aki = x_l2.reshape(q + 1, p, order='F')  # Fortran order for column-major

    return aki


def tvwlp_l2(x, p, q, w):
    """
    Time-Varying Weighted Linear Prediction with L2 norm (weighted least squares).

    Parameters
    ----------
    x : ndarray
        Input signal
    p : int
        LP order
    q : int
        Polynomial order for time-varying coefficients
    w : ndarray
        Weight function (same length as x)

    Returns
    -------
    aki : ndarray
        Time-varying LP coefficients of shape [q+1, p]
    """
    x = np.asarray(x, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    n_samples = len(x)

    # Build weighted system of equations
    n_equations = n_samples - p
    n_coeffs = p * (q + 1)

    Ypu = np.zeros((n_equations, n_coeffs), dtype=np.float64)
    Yn = np.zeros(n_equations, dtype=np.float64)

    m = 0
    for n in range(p, n_samples):
        Yn[m] = w[n] * x[n]

        for i in range(1, p + 1):
            for j in range(q + 1):
                col_idx = (i - 1) * (q + 1) + j
                Ypu[m, col_idx] = w[n] * ((n - p) ** j) * x[n - i]
        m += 1

    # Solve weighted least squares
    x_l2, residuals, rank, s = np.linalg.lstsq(Ypu, Yn, rcond=None)

    # Reshape to [q+1, p]
    aki = x_l2.reshape(q + 1, p, order='F')

    return aki


def tvlp_l1(x, p, q):
    """
    Time-Varying Linear Prediction with L1 norm (robust estimation).

    Uses linear programming to minimize L1 norm of residuals.

    Parameters
    ----------
    x : ndarray
        Input signal
    p : int
        LP order
    q : int
        Polynomial order for time-varying coefficients

    Returns
    -------
    aki : ndarray
        Time-varying LP coefficients of shape [q+1, p]
    """
    x = np.asarray(x, dtype=np.float64)
    n_samples = len(x)

    # Build system of equations
    n_equations = n_samples - p
    n_coeffs = p * (q + 1)

    Ypu = np.zeros((n_equations, n_coeffs), dtype=np.float64)
    Yn = np.zeros(n_equations, dtype=np.float64)

    m = 0
    for n in range(p, n_samples):
        Yn[m] = x[n]

        for i in range(1, p + 1):
            for j in range(q + 1):
                col_idx = (i - 1) * (q + 1) + j
                Ypu[m, col_idx] = ((n - p) ** j) * x[n - i]
        m += 1

    A = Ypu
    b = Yn

    # L1 minimization via linear programming
    # minimize: ||A*x - b||_1
    # Reformulate as: minimize sum(u+ + u-)
    # subject to: A*x - u+ + u- = b, u+>=0, u->=0
    #
    # MATLAB formulation:
    # f = [zeros(n,1); ones(m,1); ones(m,1)]
    # Aeq = [A, -eye(m), +eye(m)]
    # beq = b
    # lb = [-Inf(n,1); zeros(m,1); zeros(m,1)]

    m, n = A.shape

    # Objective: minimize sum of slack variables
    c = np.concatenate([np.zeros(n), np.ones(m), np.ones(m)])

    # Equality constraint: A*x - u+ + u- = b
    A_eq = np.hstack([A, -np.eye(m), np.eye(m)])
    b_eq = b

    # Bounds: x unbounded, u+>=0, u->=0
    bounds = [(None, None)] * n + [(0, None)] * (2 * m)

    # Solve
    result = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds,
                     method='highs', options={'maxiter': 100, 'disp': False})

    if not result.success:
        print(f"Warning: linprog did not converge. Status: {result.status}, Message: {result.message}")

    # Extract solution
    x_l1 = result.x[:n]

    # Reshape to [q+1, p]
    aki = x_l1.reshape(q + 1, p, order='F')

    return aki


def tvwlp_l1(x, p, q, w):
    """
    Time-Varying Weighted Linear Prediction with L1 norm.

    Uses linear programming to minimize weighted L1 norm of residuals.

    Parameters
    ----------
    x : ndarray
        Input signal
    p : int
        LP order
    q : int
        Polynomial order for time-varying coefficients
    w : ndarray
        Weight function (same length as x)

    Returns
    -------
    aki : ndarray
        Time-varying LP coefficients of shape [q+1, p]
    """
    x = np.asarray(x, dtype=np.float64)
    w = np.asarray(w, dtype=np.float64)
    n_samples = len(x)

    # Build weighted system
    n_equations = n_samples - p
    n_coeffs = p * (q + 1)

    Ypu = np.zeros((n_equations, n_coeffs), dtype=np.float64)
    Yn = np.zeros(n_equations, dtype=np.float64)

    m = 0
    for n in range(p, n_samples):
        Yn[m] = w[n] * x[n]

        for i in range(1, p + 1):
            for j in range(q + 1):
                col_idx = (i - 1) * (q + 1) + j
                Ypu[m, col_idx] = w[n] * ((n - p) ** j) * x[n - i]
        m += 1

    A = Ypu
    b = Yn

    # L1 minimization via linear programming (same as tvlp_l1)
    m, n = A.shape

    c = np.concatenate([np.zeros(n), np.ones(m), np.ones(m)])
    A_eq = np.hstack([A, -np.eye(m), np.eye(m)])
    b_eq = b
    bounds = [(None, None)] * n + [(0, None)] * (2 * m)

    result = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=bounds,
                     method='highs', options={'maxiter': 100, 'disp': False})

    if not result.success:
        print(f"Warning: linprog did not converge. Status: {result.status}, Message: {result.message}")

    x_l1 = result.x[:n]
    aki = x_l1.reshape(q + 1, p, order='F')

    return aki

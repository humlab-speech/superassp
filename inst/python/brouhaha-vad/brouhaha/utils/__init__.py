"""
Utilities package with automatic optimization detection.

This package automatically uses Cython or Numba optimized implementations
when available, falling back to pure Python/NumPy implementations otherwise.
"""

# Try to import optimized versions
_has_cython_collate = False
_has_cython_metrics = False
_has_numba = False

try:
    from .collate_fast import collate_y_cython, collate_y_cython_v2
    _has_cython_collate = True
except ImportError:
    pass

try:
    from .metrics_fast import (
        compute_stat_scores_cython,
        compute_mae_cython,
        binarize_with_hysteresis_cython,
        apply_binarization_full_cython
    )
    _has_cython_metrics = True
except ImportError:
    pass

try:
    from .numba_ops import is_numba_available
    _has_numba = is_numba_available()
except ImportError:
    pass


def get_optimization_status():
    """Get status of available optimizations."""
    return {
        'cython_collate': _has_cython_collate,
        'cython_metrics': _has_cython_metrics,
        'numba': _has_numba,
    }


def print_optimization_status():
    """Print available optimizations."""
    status = get_optimization_status()
    print("Brouhaha Optimization Status:")
    print("=" * 50)
    print(f"  Cython collate_y:    {' Available' if status['cython_collate'] else ' Not available'}")
    print(f"  Cython metrics:      {' Available' if status['cython_metrics'] else ' Not available'}")
    print(f"  Numba:               {' Available' if status['numba'] else ' Not available'}")
    print("=" * 50)

    if not any(status.values()):
        print("\nTo enable optimizations, install with:")
        print("  pip install -e '.[optimization]'")
        print("Then rebuild:")
        print("  python setup.py build_ext --inplace")


# Export optimization status check
__all__ = ['get_optimization_status', 'print_optimization_status']

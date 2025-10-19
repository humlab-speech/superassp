"""Voice quality measures module"""

# Try to import optimized versions first
try:
    from .harmonics_optimized import get_harmonics_optimized as get_harmonics_opt
    HARMONICS_OPTIMIZED = True
except ImportError:
    from .harmonics import get_harmonics as get_harmonics_opt
    HARMONICS_OPTIMIZED = False

try:
    from .cpp_optimized import get_cpp_optimized as get_cpp_opt
    CPP_OPTIMIZED = True
except ImportError:
    from .cpp import get_cpp as get_cpp_opt
    CPP_OPTIMIZED = False

try:
    from .formant_amplitudes_optimized import get_formant_amplitudes_optimized as get_formant_amplitudes
    FORMANT_AMP_OPTIMIZED = True
except ImportError:
    from .formant_amplitudes import get_formant_amplitudes
    FORMANT_AMP_OPTIMIZED = False

# Try parallel versions for multi-core systems
try:
    from .parallel import (
        should_use_parallel, 
        get_harmonics_parallel, 
        get_cpp_parallel,
        get_optimal_num_processes
    )
    PARALLEL_AVAILABLE = True
except ImportError:
    PARALLEL_AVAILABLE = False
    should_use_parallel = lambda: False

# Wrapper functions that choose best implementation
def get_harmonics(*args, **kwargs):
    """Get harmonics using best available implementation"""
    # Disable parallel processing due to multiprocessing issues with stdin
    # if PARALLEL_AVAILABLE and should_use_parallel() and len(args) > 2 and len(args[2]) > 1000:
    #     return get_harmonics_parallel(*args, **kwargs)
    # else:
    return get_harmonics_opt(*args, **kwargs)

def get_cpp(*args, **kwargs):
    """Get CPP using best available implementation"""
    # Disable parallel processing due to multiprocessing issues with stdin
    # if PARALLEL_AVAILABLE and should_use_parallel() and len(args) > 2 and len(args[2]) > 1000:
    #     return get_cpp_parallel(*args, **kwargs)
    # else:
    return get_cpp_opt(*args, **kwargs)

# Import remaining modules
from .harmonics import compute_harmonic_differences
from .formant_amplitudes import compute_harmonic_formant_differences
from .hnr import get_hnr
from .energy import get_energy
from .spectral import get_2k, get_5k, get_2k5k, get_h42k
from .corrections import apply_corrections

__all__ = [
    'get_harmonics',
    'get_formant_amplitudes',
    'get_cpp',
    'get_hnr',
    'get_energy',
    'get_2k',
    'get_5k',
    'get_2k5k',
    'get_h42k',
    'apply_corrections',
    'compute_harmonic_differences',
    'compute_harmonic_formant_differences',
    'HARMONICS_OPTIMIZED',
    'CPP_OPTIMIZED',
    'FORMANT_AMP_OPTIMIZED',
    'PARALLEL_AVAILABLE'
]

"""
VoiceSauce: Voice Quality Analysis Toolkit
Python reimplementation of the MATLAB VoiceSauce
"""

__version__ = '0.1.0'
__author__ = 'VoiceSauce Python Implementation'

from .core.config import VoiceSauceConfig
from .voicesauce import analyze, VoiceSauceResults

# Initialize architecture-specific optimizations
def _init_optimizations():
    """Initialize architecture-specific optimizations"""
    try:
        from .core.apple_silicon import get_optimizer
        optimizer = get_optimizer()
        
        # Configure for optimal performance
        optimizer.configure_thread_affinity()
        optimizer.optimize_fft_for_arm()
        
        # Print optimization status if verbose
        import os
        if os.environ.get('VOICESAUCE_VERBOSE'):
            print(optimizer.get_optimization_summary())
    except:
        pass  # Silently continue if optimization fails

_init_optimizations()

__all__ = [
    'VoiceSauceConfig',
    'analyze',
    'VoiceSauceResults',
    '__version__'
]

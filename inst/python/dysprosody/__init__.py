"""
Dysprosody: Prosodic assessment for speech analysis

This package implements the prosodic measures described in:
Nylén et al. (2025). A model of dysprosody in autism spectrum disorder.
Frontiers in Human Neuroscience. doi: 10.3389/fnhum.2025.1566274

The package provides:
- prosody_measures(): Extract 193 prosodic features from audio
- batch_process(): Parallel batch processing for multiple files
- MOMEL-INTSINT algorithms for pitch target extraction and tone coding

License: CC BY 4.0 (https://creativecommons.org/licenses/by/4.0/)
"""

__version__ = "1.0.0"
__author__ = "Fredrik Nylén"
__license__ = "CC BY 4.0"

# Import optimized version by default
try:
    from .dysprosody_optimized import prosody_measures, batch_process
    _OPTIMIZED = True
except ImportError:
    # Fallback to pure version if optimized unavailable
    try:
        from .dysprosody_pure import prosody_measures
        _OPTIMIZED = False

        # batch_process not available in pure version
        def batch_process(*args, **kwargs):
            raise ImportError(
                "batch_process requires dysprosody_optimized. "
                "Install with: pip install numpy pandas scipy parselmouth"
            )
    except ImportError:
        # Neither version available
        def prosody_measures(*args, **kwargs):
            raise ImportError(
                "Dysprosody requires: numpy, pandas, scipy, parselmouth. "
                "Install with: pip install numpy pandas scipy parselmouth"
            )

        def batch_process(*args, **kwargs):
            raise ImportError(
                "Dysprosody requires: numpy, pandas, scipy, parselmouth. "
                "Install with: pip install numpy pandas scipy parselmouth"
            )

        _OPTIMIZED = None

__all__ = ['prosody_measures', 'batch_process', '__version__']


def info():
    """Return information about dysprosody installation"""
    info_dict = {
        'version': __version__,
        'optimized': _OPTIMIZED,
        'dependencies': {}
    }

    # Check dependencies
    try:
        import numpy
        info_dict['dependencies']['numpy'] = numpy.__version__
    except ImportError:
        info_dict['dependencies']['numpy'] = None

    try:
        import pandas
        info_dict['dependencies']['pandas'] = pandas.__version__
    except ImportError:
        info_dict['dependencies']['pandas'] = None

    try:
        import scipy
        info_dict['dependencies']['scipy'] = scipy.__version__
    except ImportError:
        info_dict['dependencies']['scipy'] = None

    try:
        import parselmouth
        info_dict['dependencies']['parselmouth'] = parselmouth.__version__
    except ImportError:
        info_dict['dependencies']['parselmouth'] = None

    return info_dict

"""
Creaky voice detection functions.

This module implements creaky voice detection based on the Kane-Drugman method.

Notes
-----
The complete creaky voice detection algorithm requires a trained neural network
classifier, which is not included in this Python conversion. This module provides
feature extraction functions that can be used with a custom classifier.

References
----------
Kane, J., Drugman, T., Gobl, C., (2013) 'Improved automatic detection of creak',
Computer Speech and Language, 27(4), pp. 1028-1047.
"""

from .creak_features import get_all_creak_features

__all__ = ['get_all_creak_features']

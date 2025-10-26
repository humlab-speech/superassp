"""
SE_VQ: Spectral Envelope - Voice Quality GCI Detection

This module implements the SE_VQ algorithm for detecting glottal closure instants (GCIs).
It's a further development of the SEDREAMS algorithm designed to improve selection of
GCI candidates through dynamic programming and post-processing for creaky voice.

References
----------
Kane, J., Gobl, C., (2013) 'Evaluation of glottal closure instant detection in a range
of voice qualities', Speech Communication 55(2), pp. 295-314.
"""

from .se_vq import se_vq, se_vq_var_f0
from .get_mbs import get_mbs
from .get_mbs_gci_intervals import get_mbs_gci_intervals
from .rcvd_reson_gci import rcvd_reson_gci
from .search_res_interval_peaks import search_res_interval_peaks
from .reson_dyprog_mat import reson_dyprog_mat
from .gci_creak_postproc import gci_creak_postproc

__all__ = [
    'se_vq',
    'se_vq_var_f0',
    'get_mbs',
    'get_mbs_gci_intervals',
    'rcvd_reson_gci',
    'search_res_interval_peaks',
    'reson_dyprog_mat',
    'gci_creak_postproc'
]

"""Post-processing to remove false positive GCIs in creaky voice regions."""

import numpy as np


def gci_creak_postproc(gci, creak, search_reg, rep, remove_thresh=0.4, rep_num=2):
    """
    Post-process to remove false positive GCIs in creaky voice regions.

    Parameters
    ----------
    gci : array_like
        Detected GCI locations (samples)
    creak : array_like
        Binary creak decision for each sample (0=non-creaky, 1=creaky)
    search_reg : float
        Search region size (samples) around each GCI
    rep : array_like
        Resonator output signal
    remove_thresh : float, optional
        Threshold for removing weak GCIs, default: 0.4
    rep_num : int, optional
        Number of iterations for removal, default: 2

    Returns
    -------
    gci : ndarray
        Cleaned GCI locations after removing false positives

    Notes
    -----
    This function separates GCIs in creaky regions and removes those that
    are weak in the resonator output compared to their neighbors.
    """
    gci = np.asarray(gci, dtype=int)
    creak = np.asarray(creak)
    rep = np.asarray(rep).flatten()

    # Separate creaky and non-creaky GCIs
    creak_mask = creak[gci] == 1
    creak_gci = gci[creak_mask].copy()
    gci_non_creak = gci[~creak_mask].copy()

    # Process creaky GCIs
    for m in range(rep_num):
        n = 1  # Start from second GCI

        while n < len(creak_gci) - 1:
            gci_pos = creak_gci[n]
            start_idx = int(max(0, gci_pos - np.round(search_reg)))
            stop_idx = int(min(len(rep) - 1, gci_pos + np.round(search_reg)))

            # Find minimum (most negative) in resonator output around current GCI
            cur_rep_max = abs(np.min(rep[start_idx:stop_idx + 1]))

            # Compare with neighbors
            neighbor_mean = np.mean([
                abs(rep[creak_gci[n - 1]]),
                abs(rep[creak_gci[n + 1]])
            ])

            # Remove if too weak compared to neighbors
            if neighbor_mean * remove_thresh > cur_rep_max:
                creak_gci[n] = -1  # Mark for removal
                n += 2
            else:
                n += 1

        # Remove marked GCIs
        creak_gci = creak_gci[creak_gci >= 0]

    # Combine and sort all GCIs
    gci = np.sort(np.unique(np.concatenate([gci_non_creak, creak_gci])))

    return gci

"""Dynamic programming for optimal GCI path selection."""

import numpy as np


def reson_dyprog_mat(gci_rel_amp, gci_n, f0_mean, x, fs, trans_wgt=1, rel_amp_wgt=0.3):
    """
    Dynamic programming to select optimal GCI path.

    Uses the method described in Ney (1989) and previously used in the ESPS
    GCI detection algorithm. Considers target costs and transition costs to
    select the 'cheapest' path through GCI candidates.

    Parameters
    ----------
    gci_rel_amp : ndarray
        Target cost matrix, shape (N_candidates, N_frames)
        Relative amplitudes of GCI candidates
    gci_n : ndarray
        GCI candidate locations (samples), shape (N_candidates, N_frames)
    f0_mean : float
        Mean fundamental frequency (Hz)
    x : array_like
        Speech signal
    fs : float
        Sampling frequency (Hz)
    trans_wgt : float, optional
        Transition cost weight, default: 1
    rel_amp_wgt : float, optional
        Relative amplitude (target) cost weight, default: 0.3

    Returns
    -------
    gci_opt : ndarray
        Optimal GCI locations (samples), shape (N_frames,)

    References
    ----------
    H. Ney, "A comparative study of two search strategies for connected word
    recognition: Dynamic programming and heuristic search", IEEE Trans. PAMI,
    vol. 11, no. 6, pp. 586-595, 1989.

    Notes
    -----
    Original MATLAB code by John Kane @ Phonetics Lab, Trinity College Dublin,
    October 25, 2011. Python conversion: 2025.
    """
    x = np.asarray(x).flatten()

    # Initialize
    cost = (gci_rel_amp.T * rel_amp_wgt).T
    ncands, nframe = gci_n.shape
    gci_n = gci_n.T  # Transpose to (nframe, ncands)
    cost = cost.T    # Transpose to (nframe, ncands)

    prev = np.zeros((nframe, ncands), dtype=int)
    pulse_len = int(np.round(fs / f0_mean))
    gci_opt = np.zeros(nframe, dtype=int)

    # Dynamic programming forward pass
    for n in range(nframe):
        if n > 0:
            # Transition cost matrix: rows (previous), cols (current)
            costm = np.zeros((ncands, ncands))

            for c in range(ncands):
                # Transitions TO states in current frame
                start_cur = int(gci_n[n, c] - np.round(pulse_len / 2))
                stop_cur = int(gci_n[n, c] + np.round(pulse_len / 2))

                if stop_cur >= len(x):
                    stop_cur = len(x) - 1

                if start_cur < 0:
                    start_cur = 0

                pulse_cur = x[start_cur:stop_cur + 1]

                for p in range(ncands):
                    # Transitions FROM states in previous frame
                    start_prev = int(gci_n[n - 1, p] - np.round(pulse_len / 2))
                    stop_prev = int(gci_n[n - 1, p] + np.round(pulse_len / 2))

                    if start_prev < 0:
                        start_prev = 0

                    if stop_prev >= len(x):
                        stop_prev = len(x) - 1

                    pulse_prev = x[start_prev:stop_prev + 1]

                    if (len(pulse_cur) == 0 or len(pulse_prev) == 0 or
                        np.isnan(pulse_cur[0]) or np.isnan(pulse_prev[0])):
                        costm[p, c] = 0
                    else:
                        if len(pulse_cur) != len(pulse_prev):
                            cor_cur = 0
                        else:
                            # Correlation coefficient
                            if np.std(pulse_cur) > 0 and np.std(pulse_prev) > 0:
                                cor_cur = np.corrcoef(pulse_cur, pulse_prev)[0, 1]
                            else:
                                cor_cur = 0

                        # Transition cost
                        costm[p, c] = (1 - abs(cor_cur)) * trans_wgt

            # Add cumulative costs
            costm = costm + cost[n - 1, :, np.newaxis]

            # Find minimum cost path to each current state
            previ = np.argmin(costm, axis=0)
            costi = costm[previ, np.arange(ncands)]

            cost[n, :] = cost[n, :] + costi
            prev[n, :] = previ

    # Traceback to find optimal path
    best = np.zeros(nframe, dtype=int)
    best[nframe - 1] = np.argmin(cost[nframe - 1, :])

    for i in range(nframe - 1, 0, -1):
        best[i - 1] = prev[i, best[i]]

    # Extract optimal GCI locations
    for n in range(nframe):
        gci_opt[n] = int(gci_n[n, best[n]])

    return gci_opt

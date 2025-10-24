"""
MOMEL-INTSINT Implementation in Python with Parselmouth

This module implements the MOMEL (MOdelling MELody) and INTSINT (INternational Transcription System 
for INTonation) algorithms in Python using Parselmouth for pitch extraction.

Based on the original C and Perl implementations by Robert Espesser and Daniel Hirst.

References:
- Hirst, D., & Espesser, R. (1993). Automatic modelling of fundamental frequency 
  using a quadratic spline function. Travaux de l'Institut de Phonétique d'Aix, 15, 75-85.
- Hirst, D. (2005). Form and function in the representation of speech prosody. 
  Speech Communication, 46(3-4), 334-347.
"""

import numpy as np
import parselmouth
from parselmouth.praat import call
from typing import List, Tuple, Optional
import math
from dataclasses import dataclass

# Constants from original C implementation
PAS = 10.0  # in ms
PAS_TRAME = 10.0
SEUILV = 50.0  # Hz threshold for voiced
RAPP_GLITCH = 0.05
HALO_BORNE_TRAME = 4

# Constants from Perl implementation
MIN_F0 = 60  # Hz
MAX_F0 = 600  # Hz
MIN_PAUSE = 0.5  # seconds
MIN_RANGE = 0.5  # octaves
MAX_RANGE = 2.5  # octaves
STEP_RANGE = 0.1  # octaves
MEAN_SHIFT = 50  # Hertz
STEP_SHIFT = 1  # Hertz
BIG_NUMBER = 9999

# Parameters for target estimation
HIGHER = 0.5
LOWER = 0.5
UP = 0.25
DOWN = 0.25

FSIGMA = 1.0


@dataclass
class Target:
    """A pitch target point"""
    time: float  # in ms or frames
    frequency: float  # in Hz or octaves


@dataclass
class IntsintTarget:
    """An INTSINT target with tone label"""
    time: float  # in seconds
    tone: str  # INTSINT label
    target: float  # observed target in Hz
    estimate: float  # estimated target in Hz


def octave(value: float) -> float:
    """Convert Hz to octave scale"""
    return math.log(value, 2) if value > 0 else 0


def linear(value: float) -> float:
    """Convert octave scale to Hz"""
    return 2 ** value


def eliminate_glitches(hz: np.ndarray) -> np.ndarray:
    """
    Eliminate glitches in F0 contour.
    A glitch is a point higher than both neighbors by RAPP_GLITCH ratio.
    """
    hz_filtered = hz.copy()
    n = len(hz)
    
    for i in range(1, n - 1):
        if (hz[i] > hz[i - 1] * (1 + RAPP_GLITCH) and
            hz[i] > hz[i + 1] * (1 + RAPP_GLITCH)):
            hz_filtered[i] = 0.0
    
    return hz_filtered


def calc_regression(pond: np.ndarray, dpx: int, fpx: int, 
                    hzptr: np.ndarray) -> Tuple[float, float, float, np.ndarray]:
    """
    Calculate quadratic regression parameters.
    
    Returns:
        Tuple of (a0, a1, a2, hzes) where y = a0 + a1*x + a2*x^2
    """
    pn = 0.0
    sx = sx2 = sx3 = sx4 = sy = sxy = sx2y = 0.0
    
    for ix in range(dpx, fpx + 1):
        if ix < 0 or ix >= len(pond):
            continue
        p = pond[ix]
        if p != 0:
            val_ix = float(ix)
            y = hzptr[ix]
            x2 = val_ix * val_ix
            x3 = x2 * val_ix
            x4 = x2 * x2
            xy = val_ix * y
            x2y = x2 * y
            
            pn += p
            sx += p * val_ix
            sx2 += p * x2
            sx3 += p * x3
            sx4 += p * x4
            sy += p * y
            sxy += p * xy
            sx2y += p * x2y
    
    if pn < 3:
        return None, None, None, None
    
    spdxy = sxy - (sx * sy) / pn
    spdx2 = sx2 - (sx * sx) / pn
    spdx3 = sx3 - (sx * sx2) / pn
    spdx4 = sx4 - (sx2 * sx2) / pn
    spdx2y = sx2y - (sx2 * sy) / pn
    
    muet = spdx2 * spdx4 - spdx3 * spdx3
    if spdx2 == 0 or muet == 0:
        return None, None, None, None
    
    a2 = (spdx2y * spdx2 - spdxy * spdx3) / muet
    a1 = (spdxy - a2 * spdx3) / spdx2
    a0 = (sy - a1 * sx - a2 * sx2) / pn
    
    # Calculate estimated values
    hzes = np.zeros(len(hzptr))
    for ix in range(dpx, fpx + 1):
        if ix < 0 or ix >= len(hzes):
            continue
        hzes[ix] = a0 + (a1 + a2 * ix) * ix
    
    return a0, a1, a2, hzes


def cible(nval: int, hzptr: np.ndarray, lfen1: int,
         maxec: float, hzinf: float, hzsup: float) -> List[Target]:
    """
    Find pitch target candidates using quadratic spline modeling.
    
    Args:
        nval: Number of F0 values
        hzptr: F0 values in Hz
        lfen1: Window length for analysis
        maxec: Maximum error allowed
        hzinf: Lower F0 bound
        hzsup: Upper F0 bound
    
    Returns:
        List of Target objects
    """
    # Initialize weights
    pond = np.zeros(len(hzptr))
    for ix in range(nval):
        if hzptr[ix] > SEUILV:
            pond[ix] = 1.0
    
    lfens2 = lfen1 // 2
    cib = []
    
    for ix in range(nval):
        dpx = ix - lfens2
        fpx = dpx + lfen1
        nsup = 0
        nsupr = -1
        
        # Local copy of weights
        pondloc = pond.copy()
        
        # Iteratively remove outliers
        while nsup > nsupr:
            nsupr = nsup
            nsup = 0
            
            a0, a1, a2, hzes = calc_regression(pondloc, dpx, fpx, hzptr)
            
            if a0 is None:
                break
            
            for x in range(dpx, fpx + 1):
                if x < 0 or x >= len(hzptr):
                    continue
                if hzptr[x] == 0 or hzes[x] / hzptr[x] > maxec:
                    pondloc[x] = 0
                    nsup += 1
        
        # Find extremum of parabola
        xc = yc = 0.0
        if a0 is not None and a2 != 0:
            vxc = -a1 / (a2 * 2)
            if (ix - lfen1) < vxc < (ix + lfen1):
                vyc = a0 + (a1 + a2 * vxc) * vxc
                if hzinf < vyc < hzsup:
                    xc = vxc
                    yc = vyc
        
        cib.append(Target(time=xc, frequency=yc))
    
    return cib


def reduc(nval: int, lfen2: int, seuildiff_x: float, seuilrapp_y: float,
         cib: List[Target]) -> List[Target]:
    """
    Reduce target candidates by clustering and filtering.
    
    Args:
        nval: Number of original F0 values
        lfen2: Window length for reduction
        seuildiff_x: Minimum time difference threshold
        seuilrapp_y: Minimum frequency ratio threshold
        cib: List of target candidates
    
    Returns:
        Reduced list of targets
    """
    # Remove zero targets
    cib_valid = [c for c in cib if c.frequency > 0]
    if not cib_valid:
        return []
    
    # Sort by time
    cib_sorted = sorted(cib_valid, key=lambda c: c.time)
    
    lf = lfen2 // 2
    
    # Calculate distances
    xdist = []
    ydist = []
    
    for i in range(len(cib_sorted) - 1):
        j1 = max(0, i - lf)
        j2 = min(len(cib_sorted), i + lf + 1)
        
        # Left group
        left = cib_sorted[j1:i + 1]
        if left:
            sxg = sum(c.time for c in left) / len(left)
            syg = sum(c.frequency for c in left) / len(left)
        else:
            sxg = syg = 0
        
        # Right group
        right = cib_sorted[i + 1:j2]
        if right:
            sxd = sum(c.time for c in right) / len(right)
            syd = sum(c.frequency for c in right) / len(right)
        else:
            sxd = syd = 0
        
        if left and right:
            xdist.append(abs(sxg - sxd))
            ydist.append(abs(syg - syd))
        else:
            xdist.append(-1)
            ydist.append(-1)
    
    # Weight distances
    valid_dists = [(x, y) for x, y in zip(xdist, ydist) if x > 0]
    if not valid_dists:
        return cib_sorted[:1] if cib_sorted else []
    
    xds = sum(x for x, y in valid_dists)
    yds = sum(y for x, y in valid_dists)
    
    if xds == 0 or yds == 0:
        return cib_sorted[:1] if cib_sorted else []
    
    px = len(valid_dists) / xds
    py = len(valid_dists) / yds
    
    # Combined distance
    dist = []
    for x, y in zip(xdist, ydist):
        if x > 0:
            dist.append((x * px + y * py) / (px + py))
        else:
            dist.append(-1)
    
    # Threshold for peak detection
    seuil = 2.0 / (px + py)
    
    # Find peaks above threshold
    partitions = [0]
    susseuil = False
    xmax = 0
    
    for i in range(len(dist)):
        if not susseuil:
            if dist[i] > seuil:
                susseuil = True
                xmax = i
        else:
            if dist[i] > dist[xmax]:
                xmax = i
            if dist[i] < seuil:
                partitions.append(xmax)
                susseuil = False
    
    if susseuil:
        partitions.append(xmax)
    partitions.append(len(cib_sorted))
    
    # Create reduced targets by averaging within partitions
    cibred = []
    for ip in range(len(partitions) - 1):
        parinf = partitions[ip]
        parsup = partitions[ip + 1]
        
        group = cib_sorted[parinf:parsup]
        if group:
            # Sigma filtering
            sx = sum(c.time for c in group)
            sx2 = sum(c.time ** 2 for c in group)
            sy = sum(c.frequency for c in group)
            sy2 = sum(c.frequency ** 2 for c in group)
            n = len(group)
            
            if n > 1:
                xm = sx / n
                ym = sy / n
                varx = max(0.1, sx2 / n - xm ** 2)
                vary = max(0.1, sy2 / n - ym ** 2)
                
                et2x = FSIGMA * math.sqrt(varx)
                et2y = FSIGMA * math.sqrt(vary)
                
                # Filter outliers
                group_filtered = [
                    c for c in group
                    if (xm - et2x <= c.time <= xm + et2x and
                        ym - et2y <= c.frequency <= ym + et2y)
                ]
                
                if group_filtered:
                    group = group_filtered
            
            # Average
            if group:
                avg_time = sum(c.time for c in group) / len(group)
                avg_freq = sum(c.frequency for c in group) / len(group)
                cibred.append(Target(time=avg_time, frequency=avg_freq))
    
    # Second filtering: merge close targets
    if not cibred:
        return []
    
    cibred2 = [cibred[0]]
    for i in range(1, len(cibred)):
        if cibred[i].time - cibred2[-1].time < seuildiff_x:
            # Check frequency ratio
            if (abs(cibred[i].frequency - cibred2[-1].frequency) /  
                cibred2[-1].frequency < seuilrapp_y):
                # Merge
                cibred2[-1] = Target(
                    time=(cibred2[-1].time + cibred[i].time) / 2,
                    frequency=(cibred2[-1].frequency + cibred[i].frequency) / 2
                )
            # Otherwise keep the one with more weight (skip for now, just keep first)
        else:
            cibred2.append(cibred[i])
    
    return cibred2


def borne(nval: int, cibred2: List[Target], hzptr: np.ndarray) -> List[Target]:
    """
    Add boundary targets at the beginning and end of the contour.
    
    Args:
        nval: Number of F0 values
        cibred2: Reduced list of targets
        hzptr: F0 values
    
    Returns:
        List of targets with added boundaries
    """
    if not cibred2:
        return []
    
    result = []
    
    # Find first voiced frame
    first_voiced = 0
    for i in range(nval):
        if hzptr[i] >= SEUILV:
            first_voiced = i
            break
    
    # Add left boundary if needed
    if cibred2[0].time > first_voiced + HALO_BORNE_TRAME:
        ancre = cibred2[0]
        
        # Quadratic regression from anchor back to start
        sx2y = 0.0
        sx4 = 0.0
        j = 0
        
        for i in range(int(ancre.time), -1, -1):
            if i < len(hzptr) and hzptr[i] > SEUILV:
                x2 = float(j * j)
                sx2y += x2 * (hzptr[i] - ancre.frequency)
                sx4 += x2 * x2
            j += 1
        
        if sx4 != 0:
            a = sx2y / sx4
            frontiere = float(first_voiced)
            
            borne_x = frontiere - (ancre.time - frontiere)
            borne_y = ancre.frequency + 2 * a * (ancre.time - frontiere) ** 2
            
            result.append(Target(time=borne_x, frequency=borne_y))
    
    # Add all targets
    result.extend(cibred2)
    
    # Find last voiced frame
    last_voiced = nval - 1
    for i in range(nval - 1, -1, -1):
        if hzptr[i] >= SEUILV:
            last_voiced = i
            break
    
    # Add right boundary if needed
    if cibred2[-1].time < last_voiced - HALO_BORNE_TRAME:
        ancre = cibred2[-1]
        
        # Quadratic regression from anchor forward to end
        sx2y = 0.0
        sx4 = 0.0
        j = 0
        
        for i in range(int(ancre.time), nval):
            if i < len(hzptr) and hzptr[i] > SEUILV:
                x2 = float(j * j)
                sx2y += x2 * (hzptr[i] - ancre.frequency)
                sx4 += x2 * x2
            j += 1
        
        if sx4 != 0:
            a = sx2y / sx4
            frontiere = float(last_voiced)
            
            borne_x = frontiere + (frontiere - ancre.time)
            borne_y = ancre.frequency + 2 * a * (ancre.time - frontiere) ** 2
            
            result.append(Target(time=borne_x, frequency=borne_y))
    
    return result


def momel(f0_values: np.ndarray, window_length: int = 30,
         min_f0: float = 60, max_f0: float = 750,
         max_error: float = 1.04,
         reduced_window_length: int = 20,
         minimal_distance: float = 5.0,
         minimal_frequency_ratio: float = 0.05) -> List[Target]:
    """
    MOMEL algorithm: Extract pitch targets from F0 contour.
    
    Args:
        f0_values: F0 values in Hz (sampled at 10ms intervals)
        window_length: Window length in samples for cible
        min_f0: Minimum F0 in Hz
        max_f0: Maximum F0 in Hz
        max_error: Maximum error allowed
        reduced_window_length: Window length for reduction
        minimal_distance: Minimal distance in frames
        minimal_frequency_ratio: Minimal frequency ratio
    
    Returns:
        List of Target objects with time in frames and frequency in Hz
    """
    nval = len(f0_values)
    
    # Constrain F0 values
    hz = np.clip(f0_values, MIN_F0, MAX_F0)
    
    # Eliminate glitches
    hz = eliminate_glitches(hz)
    
    # Find targets
    targets = cible(nval, hz, window_length, max_error, min_f0, max_f0)
    
    # Reduce targets
    targets_reduced = reduc(nval, reduced_window_length, minimal_distance, 
                           minimal_frequency_ratio, targets)
    
    # Add boundaries
    targets_final = borne(nval, targets_reduced, hz)
    
    return targets_final


def estimate_target(tone: str, last_target: float, top: float, 
                   bottom: float, mid: float) -> float:
    """Estimate target F0 for a given INTSINT tone."""
    if tone == "M":
        return mid
    elif tone == "S":
        return last_target
    elif tone == "T":
        return top
    elif tone == "H":
        return last_target + (top - last_target) * HIGHER
    elif tone == "U":
        return last_target + (top - last_target) * UP
    elif tone == "B":
        return bottom
    elif tone == "L":
        return last_target - (last_target - bottom) * LOWER
    elif tone == "D":
        return last_target - (last_target - bottom) * DOWN
    else:
        return mid


def optimize_intsint(targets: List[Target], mid: float, 
                     range_oct: float) -> Tuple[List[str], List[float], float]:
    """
    Optimize INTSINT labels for given targets, mid point, and range.
    
    Returns:
        Tuple of (tone labels, estimated F0s, sum-squared error)
    """
    top = mid + range_oct / 2
    bottom = mid - range_oct / 2
    
    # Convert targets to octave scale
    targets_oct = [octave(t.frequency) for t in targets if t.frequency > 0]
    if not targets_oct:
        return [], [], BIG_NUMBER
    
    # List of tones (all except M for non-first targets)
    list_tones = ["M", "S", "T", "B", "H", "L", "U", "D"]
    
    intsint = []
    estimate = []
    last_estimate = 0
    ss_error = 0
    
    # Process each target
    last_time = -BIG_NUMBER
    for i, target in enumerate(targets):
        if target.frequency == 0:
            continue
        
        target_oct = octave(target.frequency)
        
        # After pause, choose from [M, T, B]
        if i == 0 or (target.time - last_time > MIN_PAUSE * 1000 / PAS_TRAME):
            if top - target_oct < abs(target_oct - mid):
                tone = "T"
            elif target_oct - bottom < abs(target_oct - mid):
                tone = "B"
            else:
                tone = "M"
        else:
            # Choose any tone except M
            min_difference = BIG_NUMBER
            best_tone = "S"
            
            for tone_test in list_tones:
                if tone_test != "M":
                    est = estimate_target(tone_test, last_estimate, top, 
                                        bottom, mid)
                    difference = abs(target_oct - est)
                    if difference < min_difference:
                        min_difference = difference
                        best_tone = tone_test
            
            tone = best_tone
        
        intsint.append(tone)
        est = estimate_target(tone, last_estimate, top, bottom, mid)
        estimate.append(est)
        
        error = abs(est - target_oct)
        ss_error += error ** 2
        
        last_estimate = est
        last_time = target.time
    
    return intsint, estimate, ss_error


def intsint(targets: List[Target]) -> Tuple[List[IntsintTarget], float, float]:
    """
    INTSINT algorithm: Assign optimal tone labels to MOMEL targets.
    
    Args:
        targets: List of Target objects from MOMEL
    
    Returns:
        Tuple of (INTSINT targets, optimal range, optimal key)
    """
    # Convert to octave scale
    valid_targets = [t for t in targets if t.frequency > 0]
    if not valid_targets:
        return [], 0, 0
    
    f0_oct = [octave(t.frequency) for t in valid_targets]
    mean_f0 = sum(f0_oct) / len(f0_oct)
    
    linear_mean_f0 = round(linear(mean_f0))
    
    min_mean = linear_mean_f0 - MEAN_SHIFT
    max_mean = linear_mean_f0 + MEAN_SHIFT
    min_ss_error = BIG_NUMBER
    
    best_intsint = []
    best_estimate = []
    best_mid = mean_f0
    best_range = MIN_RANGE
    
    # Optimize over range and key
    for range_oct in np.arange(MIN_RANGE, MAX_RANGE, STEP_RANGE):
        for lm in range(int(min_mean), int(max_mean) + 1, STEP_SHIFT):
            mid = octave(lm)
            
            intsint_labels, estimates, ss_error = optimize_intsint(
                valid_targets, mid, range_oct)
            
            if ss_error < min_ss_error:
                min_ss_error = ss_error
                best_range = range_oct
                best_mid = mid
                best_intsint = intsint_labels
                best_estimate = estimates
    
    # Create INTSINT targets
    result = []
    for i, target in enumerate(valid_targets):
        if i < len(best_intsint):
            result.append(IntsintTarget(
                time=target.time * PAS_TRAME / 1000,  # Convert to seconds
                tone=best_intsint[i],
                target=target.frequency,
                estimate=linear(best_estimate[i])
            ))
    
    return result, best_range, linear(best_mid)


def extract_f0_parselmouth(sound_file: str, time_step: float = 0.01,
                          pitch_floor: float = 60, pitch_ceiling: float = 750,
                          silence_threshold: float = 0.03,
                          voicing_threshold: float = 0.45,
                          octave_cost: float = 0.01,
                          octave_jump_cost: float = 0.35,
                          voiced_unvoiced_cost: float = 0.14) -> np.ndarray:
    """
    Extract F0 values using Parselmouth (Praat).
    
    Args:
        sound_file: Path to sound file
        time_step: Time step in seconds
        pitch_floor: Minimum pitch in Hz
        pitch_ceiling: Maximum pitch in Hz
        silence_threshold: Silence threshold
        voicing_threshold: Voicing threshold  
        octave_cost: Octave cost
        octave_jump_cost: Octave jump cost
        voiced_unvoiced_cost: Voiced/unvoiced cost
    
    Returns:
        Array of F0 values in Hz (0 for unvoiced)
    """
    # Load sound
    sound = parselmouth.Sound(sound_file)
    
    # Extract pitch
    pitch = call(sound, "To Pitch", time_step, pitch_floor, pitch_ceiling)
    
    # Extract F0 values
    n_frames = call(pitch, "Get number of frames")
    f0_values = []
    
    for i in range(1, n_frames + 1):
        f0 = call(pitch, "Get value in frame", i, "Hertz")
        if f0 is None or np.isnan(f0):
            f0_values.append(0.0)
        else:
            f0_values.append(f0)
    
    return np.array(f0_values)


def process_momel_intsint(sound_file: str,
                          window_length: int = 30,
                          min_f0: float = 60,
                          max_f0: float = 750,
                          pitch_span: float = 1.5,
                          max_error: float = 1.04,
                          reduced_window_length: int = 20,
                          minimal_distance: int = 20,
                          minimal_frequency_ratio: float = 0.05,
                          time_step: float = 0.01) -> Tuple[List[IntsintTarget], 
                                                             float, float]:
    """
    Complete MOMEL-INTSINT processing pipeline.
    
    Args:
        sound_file: Path to sound file
        window_length: Window length in ms for MOMEL
        min_f0: Minimum F0 in Hz
        max_f0: Maximum F0 in Hz
        pitch_span: Pitch span in octaves
        max_error: Maximum error for MOMEL
        reduced_window_length: Reduced window length in ms
        minimal_distance: Minimal distance in ms
        minimal_frequency_ratio: Minimal frequency ratio
        time_step: Time step for pitch extraction in seconds
    
    Returns:
        Tuple of (INTSINT targets, range, key)
    """
    # Extract F0
    f0_values = extract_f0_parselmouth(
        sound_file, 
        time_step=time_step,
        pitch_floor=min_f0,
        pitch_ceiling=max_f0
    )
    
    # Run MOMEL
    momel_targets = momel(
        f0_values,
        window_length=int(window_length / PAS_TRAME),
        min_f0=min_f0,
        max_f0=max_f0,
        max_error=max_error,
        reduced_window_length=int(reduced_window_length / PAS_TRAME),
        minimal_distance=minimal_distance / PAS_TRAME,
        minimal_frequency_ratio=minimal_frequency_ratio
    )
    
    # Run INTSINT
    intsint_targets, range_oct, key = intsint(momel_targets)
    
    return intsint_targets, range_oct, key


# For testing
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        sound_file = sys.argv[1]
        results, range_val, key_val = process_momel_intsint(sound_file)
        
        print(f"; INTSINT labels")
        print(f";   {len(results)} values  mean = {key_val:.1f}")
        print(f"<parameter range={range_val:.2f}>")
        print(f"<parameter key={key_val:.1f}>")
        
        for target in results:
            print(f"{target.time:.3f} {target.tone} "
                  f"{target.target:.0f} {target.estimate:.0f}")

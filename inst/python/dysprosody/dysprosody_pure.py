"""
Pure Python implementation of prosody_measures using momel_intsint.py

This module provides an optimized, pure Python version of the prosody_measures function
that eliminates dependencies on compiled C binaries and Perl scripts.

Based on the original dysprosody.py implementation, but using the pure Python
MOMEL-INTSINT implementation from momel_intsint.py.
"""

from pathlib import Path
import parselmouth
import math
import numpy as np
import pandas as pd
import io
from scipy import stats

# Import pure Python MOMEL-INTSINT implementation
from momel_intsint import momel as momel_pure, intsint, extract_f0_parselmouth, PAS_TRAME, linear


def getSound(file, beginTime=0.0, endTime=0.0):
    """Load sound file using Parselmouth"""
    name = Path(file).stem
    if beginTime <= 0.0 and endTime <= 0.0:
        sound = parselmouth.Sound(file)
    else:
        beginTime = max(beginTime, 0.0)
        endTime = max(endTime, 0.0)
        ls = parselmouth.praat.call("Open long sound file", file)
        sound = parselmouth.praat.call(ls, "Extract part", beginTime, endTime)
    parselmouth.praat.call(sound, "Rename", name)
    return sound


def automatic_min_max_fo(sound, maximum_pitch_span=1.5, minimum_fo=60, maximum_fo=750):
    """Automatic F0 range estimation using two-pass approach"""
    firstPass = sound.to_pitch(0.01, minimum_fo, maximum_fo)
    q25 = parselmouth.praat.call(firstPass, "Get quantile", 0, 0, 0.25, "Hertz")
    min_fo = math.floor(q25 * 0.75 / 10) * 10
    max_fo = math.ceil((min_fo * pow(2, maximum_pitch_span)) / 10) * 10
    secondPass = sound.to_pitch(0.01, min_fo, max_fo)
    parselmouth.praat.call(secondPass, "Rename", str(sound.name))
    return (secondPass, min_fo, max_fo)


def momel_to_pitch_tier(sound, pitchobj, momel_targets, minimum_fo=60, maximum_fo=750):
    """
    Convert MOMEL targets to PitchTier and return momel_pitch_values in the expected format.

    This function replicates the output format of the original subprocess-based momel() function.
    """
    # Create PitchTier
    momelPitchTier = parselmouth.praat.call(
        "Create PitchTier",
        pitchobj.name,
        0.0,
        float(pitchobj.duration)
    )

    # Convert targets to the format expected by the original code
    # Original format: list of tuples (time_ms_str, freq_hz_str)
    momel_pitch_values = []

    for target in momel_targets:
        if target.frequency > 0:
            # Time in milliseconds (as in original)
            time_ms = target.time * PAS_TRAME
            # Constrain frequency to bounds (as in original)
            pitchValue = min(max(float(target.frequency), float(minimum_fo)), float(maximum_fo))
            timeValue = time_ms / 1000.0  # Convert to seconds for PitchTier

            parselmouth.praat.call(momelPitchTier, "Add point", timeValue, pitchValue)

            # Store in original format: (time_ms_str, freq_hz_str)
            momel_pitch_values.append((str(time_ms), str(pitchValue)))

    parselmouth.praat.call(momelPitchTier, "Rename", str(pitchobj.name))
    return (momelPitchTier, momel_pitch_values)


def code_with_intsint_pure(sound, momel_targets, intsint_targets, range_oct, key_hz):
    """
    Pure Python version of code_with_intsint using momel_intsint.py output.

    Creates a TextGrid with three tiers: Momel, Intsint, IntsintMomel
    and returns the same optimization parameters as the original.
    """
    epsilon_time = 0.0001

    # Create TextGrid with three point tiers
    tgMomelTier = 1
    tgIntsintTier = 2
    tgMomelIntsintTier = 3

    tg = parselmouth.praat.call(
        sound, "To TextGrid",
        "Momel Intsint IntsintMomel",
        "Momel Intsint IntsintMomel"
    )
    parselmouth.praat.call(tg, "Rename", str(sound.name))

    # Process INTSINT targets
    prev_time = 0.0
    nvalues = len(intsint_targets)

    # Calculate mean pitch from MOMEL targets (as the Perl script does)
    # The Perl script calculates mean from the input targets in octave scale,
    # then converts back to Hz
    if len(momel_targets) > 0:
        from momel_intsint import octave, linear
        valid_targets = [t for t in momel_targets if t.frequency > 0]
        if valid_targets:
            f0_oct = [octave(t.frequency) for t in valid_targets]
            mean_f0_oct = sum(f0_oct) / len(f0_oct)
            pmean = round(linear(mean_f0_oct))
        else:
            pmean = key_hz
    else:
        pmean = 0.0

    # The optimization parameters from the INTSINT algorithm
    # In the pure Python version, these are simplified since we don't have
    # the detailed optimization trace from the Perl script
    # We reconstruct what the Perl script would have output

    for target in intsint_targets:
        time = target.time  # Already in seconds

        # Handle epsilon time to avoid duplicate times
        if abs(time - prev_time) <= epsilon_time:
            time += epsilon_time
        prev_time = time

        # Insert points into the three tiers
        # Tier 1: Momel - the observed pitch
        parselmouth.praat.call(tg, "Insert point", tgMomelTier, time, str(int(target.target)))

        # Tier 2: Intsint - the tone label
        parselmouth.praat.call(tg, "Insert point", tgIntsintTier, time, str(target.tone))

        # Tier 3: IntsintMomel - combination like "T192"
        intsintmomel = f"{target.tone}{int(target.estimate)}"
        parselmouth.praat.call(tg, "Insert point", tgMomelIntsintTier, time, intsintmomel)

    # Calculate optimization parameters to match original output
    # The Perl script outputs range and mean optimization details
    # We approximate these based on the pure Python implementation

    # Range optimization (in octaves)
    prange = range_oct
    orlow = 0.5  # MIN_RANGE from momel_intsint.py
    orhigh = 2.5  # MAX_RANGE
    orstep = 0.1  # STEP_RANGE

    # Mean optimization (in Hz)
    # The Perl script searches mean ± MEAN_SHIFT
    omlow = key_hz - 50  # MEAN_SHIFT = 50
    omhigh = key_hz + 50
    omstep = 1  # STEP_SHIFT

    return (tg, nvalues, pmean, orlow, orhigh, orstep, omlow, omhigh, omstep, prange, key_hz)


def correction_iseli_i(f, F_i, B_i, fs):
    """
    Return the i-th correction (dB) to the harmonic amplitude using the
    algorithm developed by Iseli and Alwan.

    Implementation from OpenSauce/dysprosody.py
    """
    # These variable names are from the Iseli-Alwan paper
    # Normalize frequencies to sampling frequency
    rooInd = min([F_i.size, B_i.size, f.size])  # Safeguard against incompatible lengths
    corrPadLength = max([F_i.size, B_i.size, f.size]) - rooInd
    r_i = np.exp(-np.pi * B_i[:rooInd] / fs)
    omega_i = 2 * np.pi * F_i[:rooInd] / fs
    omega = 2 * np.pi * f[:rooInd] / fs

    # Factors needed to compute correction
    numerator_sqrt = r_i**2 + 1 - 2 * r_i * np.cos(omega_i)
    denom_factor1 = r_i**2 + 1 - 2 * r_i * np.cos(omega_i + omega)
    denom_factor2 = r_i**2 + 1 - 2 * r_i * np.cos(omega_i - omega)

    # Correction in the z-domain
    corr_i = (20 * np.log10(numerator_sqrt) -
              10 * np.log10(denom_factor1) -
              10 * np.log10(denom_factor2))

    # Assume correction to be zero if too few formants have been identified
    corr_i = np.pad(corr_i, (0, corrPadLength), 'constant', constant_values=(0, 0))

    return corr_i


def bandwidth_hawks_miller(F_i, F0):
    """
    Return formant bandwidth estimated from the formant frequency and the
    fundamental frequency.

    Implementation from OpenSauce/dysprosody.py
    """
    # Bandwidth scaling factor as a function of F0
    S = 1 + 0.25 * (F0 - 132) / 88

    # Coefficients C1 (for F_i < 500 Hz) and C2 (F_i >= 500 Hz)
    C1 = np.array([165.327516, -6.73636734e-1, 1.80874446e-3,
                   -4.52201682e-6, 7.49514000e-9, -4.70219241e-12])
    C2 = np.array([15.8146139, 8.10159009e-2, -9.79728215e-5,
                   5.28725064e-8, -1.07099364e-11, 7.91528509e-16])

    # Construct matrix that is a 5th order power series of the formant frequency
    F_i_mat = np.vstack((F_i**0, F_i**1, F_i**2, F_i**3, F_i**4, F_i**5))

    # Construct mask for formant frequency < 500 Hz
    F_i_dummy = F_i.copy()
    F_i_dummy[np.isnan(F_i_dummy)] = 0
    mask_less_500 = np.tile(F_i_dummy < 500, (len(C1), 1))

    # Formant bandwidth estimation
    B_i = S * (np.dot(C1, F_i_mat * mask_less_500) +
               np.dot(C2, F_i_mat * np.logical_not(mask_less_500)))

    return B_i


def spectral_tilt(sound, momel_pitch, formantObj, time, windowSize=60,
                 minimum_fo=60, maximum_fo=750):
    """
    Compute spectral tilt measures at a specific time point.

    This is identical to the original implementation in dysprosody.py
    """
    pitch = momel_pitch
    soundDur = parselmouth.praat.call(sound, "Get total duration")

    windowSize = float(windowSize / 1000)  # in seconds
    beginTime = max(time - (windowSize / 2), 0.0)
    endTime = min(time + (windowSize / 2), soundDur)
    soundpart = parselmouth.praat.call(sound, "Extract part", beginTime, endTime,
                                      "Gaussian1", 1.0, False)

    # C1 computation
    mfccSound = parselmouth.praat.call(soundpart, "Resample", 4000, 50)
    stepSize = 0.005
    mfcc = parselmouth.praat.call(mfccSound, "To MFCC", 12, 0.015, stepSize, 118, 118, 0.0)
    mfccFrame = math.ceil(parselmouth.praat.call(mfcc, "Get number of frames") / 2)
    C1 = parselmouth.praat.call(mfcc, "Get value in frame", mfccFrame, 1)

    spct = parselmouth.praat.call(soundpart, "To Spectrum", "yes")
    ltas = parselmouth.praat.call(spct, "To Ltas (1-to-1)")

    Fi = [parselmouth.praat.call(formantObj, "Get value at time", i + 1, time,
                                 "Hertz", "Linear")
          for i in range(parselmouth.praat.call(formantObj, "Get minimum number of formants"))]

    Bi = bandwidth_hawks_miller(np.array(Fi), pitch)

    # Spectral balance
    spectralbalance = parselmouth.praat.call(ltas, "Get slope", 0, 500, 500, 1000, "energy")

    # SLF / Spectral tilt
    slR = parselmouth.praat.call(ltas, "Report spectral trend", 100, 5000,
                                 "logarithmic", "least squares")
    SLF = float(slR.split()[11])

    # SLF6D coefficients
    soundpart.to_spectrogram()
    spectrogram = sound.to_spectrogram()
    log_magnitude_spectrum = np.log1p(spectrogram.values)
    mean_log_magnitude = np.mean(log_magnitude_spectrum, axis=1)
    SLF6D_coefficients = np.polyfit(np.arange(len(mean_log_magnitude)),
                                   mean_log_magnitude, 6)
    SLF6D_dict = zip(["SLF6D." + str(i) for i in range(1, 7)], SLF6D_coefficients)

    # Harmonic analysis
    lowerbh = [(pitch * (n + 1)) - (pitch / 10) for n in range(4)]
    upperbh = [(pitch * (n + 1)) + (pitch / 10) for n in range(4)]

    Ln = [parselmouth.praat.call(ltas, "Get maximum", lowerbh[i], upperbh[i], "Parabolic")
          for i in range(4)]
    L_Fn = [parselmouth.praat.call(ltas, "Get value at frequency", Fi[i], "Linear")
            for i in range(len(Fi))]
    nfo = [parselmouth.praat.call(ltas, "Get frequency of maximum", lowerbh[i],
                                  upperbh[i], "Parabolic")
           for i in range(4)]

    # Formant correction
    corr = correction_iseli_i(np.array(nfo[0:2] + nfo[3:]), np.array(Fi[0:3]),
                             Bi[0:3], sound.sampling_frequency)

    Ln_c = Ln
    for i in range(len(corr)):
        Ln_c = Ln_c - corr[i]

    corr = correction_iseli_i(np.array(L_Fn), np.array(Fi), Bi, sound.sampling_frequency)
    L_Fn_c = L_Fn
    for i in range(len(corr)):
        L_Fn_c = L_Fn_c - corr[i]

    # H1-H2 measures
    L2L1 = Ln[0] - Ln[1]
    L2cL1c = Ln_c[0] - Ln_c[1]
    L1cLF3c = float(Ln_c[0] - L_Fn_c[2]) if len(L_Fn_c) > 2 else None
    L1LF3 = float(Ln[0] - L_Fn[2]) if len(L_Fn) > 2 else None

    out = {
        'L2L1': float(L2L1),
        'L2cL1c': float(L2cL1c),
        'L1cLF3c': L1cLF3c,
        'L1LF3': L1LF3,
        'SLF': float(SLF),
        'C1': float(C1),
        'Spectral Balance': float(spectralbalance)
    }
    out.update(SLF6D_dict)
    return out


def prosody_measures(soundPath, minF=60, maxF=750, windowShift=10):
    """
    Pure Python implementation of prosody_measures.

    Computes comprehensive prosodic features from audio files using only Python,
    eliminating dependencies on compiled C binaries and Perl scripts.

    Args:
        soundPath: Path to audio file
        minF: Minimum F0 in Hz (default: 60)
        maxF: Maximum F0 in Hz (default: 750)
        windowShift: Window shift in ms (default: 10)

    Returns:
        pandas.Series with prosodic features, or None if file < 1 second
    """
    windowShift = float(windowShift / 1000)  # in seconds
    soundObj = getSound(soundPath)
    print("Processing: " + soundPath)

    duration = float(parselmouth.praat.call(soundObj, "Get total duration"))
    if duration < 1.0:
        print("Skipping utterance with less than 1 second duration: " + soundPath)
        return None

    # Step 1: Automatic F0 range estimation
    pitchObj, min_fo, max_fo = automatic_min_max_fo(
        soundObj,
        maximum_pitch_span=1.5,
        minimum_fo=minF,
        maximum_fo=maxF
    )

    # Step 2: Extract F0 values for MOMEL
    pitchvalues = parselmouth.praat.call(pitchObj, "List values in all frames", "Hertz")
    np.nan_to_num(pitchvalues, nan=0.0, copy=False)

    # Step 3: Run pure Python MOMEL
    # Parameters match the original dysprosody.py subprocess call
    # window_length, reduced_window_length, minimal_distance are in frames
    momel_targets = momel_pure(
        pitchvalues,
        window_length=int(30 / PAS_TRAME),  # 30 ms / 10 ms = 3 frames
        min_f0=min_fo,
        max_f0=max_fo,
        max_error=1.04,
        reduced_window_length=int(20 / PAS_TRAME),  # 20 ms / 10 ms = 2 frames
        minimal_distance=5.0,  # 5 frames (not divided by PAS_TRAME!)
        minimal_frequency_ratio=0.05
    )

    # Step 4: Run pure Python INTSINT
    intsint_targets, range_oct, key_hz = intsint(momel_targets)

    # Step 5: Convert to PitchTier and format for compatibility
    momelPitchTier, momel_pitch_values = momel_to_pitch_tier(
        soundObj, pitchObj, momel_targets, min_fo, max_fo
    )

    # Step 6: Create TextGrid with INTSINT coding
    textGrid, nintsint, pitchmean, optim_r_low, optim_r_high, optim_r_step, \
        optim_m_low, optim_m_high, optim_m_step, pitchrange, key = \
        code_with_intsint_pure(soundObj, momel_targets, intsint_targets, range_oct, key_hz)

    # Step 7: Extract formants and intensity
    formantObj = parselmouth.praat.call(soundObj, "To Formant (burg)",
                                       0.001, 5, 5000, 0.025, 50)
    intensityObj = parselmouth.praat.call(soundObj, "To Intensity",
                                         min_fo, windowShift, "yes")

    # Step 8: Convert TextGrid to DataFrame
    tgTable = pd.read_table(io.StringIO(
        parselmouth.praat.call(textGrid, "List", False, 20, True, False)
    ))
    tgTabWide = pd.pivot(tgTable[['tmin', 'tier', 'text']],
                        index="tmin", columns="tier", values="text")

    # Step 9: Add intensity values
    tgTabWide['Intensity'] = tgTabWide.apply(
        lambda x: parselmouth.praat.call(intensityObj, "Get value at time",
                                        x.name, "cubic"),
        axis=1
    )

    # Step 10: Merge in computed momel_pitch values
    momelPitch = pd.DataFrame(
        [(float(a) / 1000.0, float(b)) for a, b in momel_pitch_values],
        columns=["time", "momel_pitch"],
        index=pd.Index([float(a) / 1000 for a, b in momel_pitch_values], name="tmin")
    )
    tgTabWide = pd.merge_asof(tgTabWide, momelPitch, left_index=True, right_index=True)
    tgTabWide.dropna(axis=0, inplace=True)

    # Step 11: Merge in spectral tilt measures
    specTilt = tgTabWide.apply(
        lambda x: spectral_tilt(soundObj, float(x['momel_pitch']),
                               formantObj, x.name),
        axis=1,
        result_type="expand"
    )
    tgTabWide = tgTabWide.join(specTilt)

    # Step 12: Compute inter-INTSINT-label differences
    tgTabWide_toDiff = tgTabWide[tgTabWide.columns[3:18]]

    # Apply differentiation
    tgTabDiff = pd.DataFrame(np.diff(tgTabWide_toDiff, axis=0),
                            columns=tgTabWide_toDiff.columns)
    tgTabDiff = pd.concat([
        pd.DataFrame(np.zeros((1, tgTabWide_toDiff.shape[1])),
                    columns=tgTabWide_toDiff.columns),
        tgTabDiff
    ], ignore_index=True)

    tgTabDiff.columns = tgTabDiff.columns + "_diff"
    tgTabDiff.index = tgTabWide_toDiff.index

    tgWideTab = tgTabDiff.join(tgTabWide_toDiff)

    # Step 13: Compute statistics
    def safe_statistics(series):
        results = {}
        try:
            results['tstd'] = stats.tstd(series)
            results['tmean'] = stats.tmean(series)
            results['variation'] = stats.variation(series)
            results['iqr'] = stats.iqr(series)
            results['tmax'] = stats.tmax(series)
            results['tmin'] = stats.tmin(series)
        except ValueError:
            results = {
                'tstd': None, 'tmean': None, 'variation': None,
                'iqr': None, 'tmax': None, 'tmin': None
            }
        return pd.Series(results)

    tgTabSummary = tgWideTab.apply(lambda x: safe_statistics(x),
                                  axis=0, result_type='expand')
    tgTabSummary.columns = ['_'.join(col).strip()
                           for col in tgTabSummary.columns.values]

    tgTabSummary.reset_index(names="statistic", inplace=True)
    tgTabSummary_long = tgTabSummary.melt(id_vars="statistic")
    tgTabSummary_long["variable_name"] = (tgTabSummary_long.variable + "_" +
                                         tgTabSummary_long.statistic)

    # Fix momel_pitch problem and rename indices
    tgTabSummary.rename(
        index={"tstd": "std", "tmean": "mean", "hdmedian": "median",
               "variation": "var", "iqr": "iqr", "tmax": "max", "tmin": "min"},
        columns={"momel_pitch": "momelpitch"},
        inplace=True
    )
    tgTabSummary.columns = tgTabSummary.columns.to_flat_index() + "_diff"

    # Step 14: Prepare results dictionary
    tgSummary = dict(zip(tgTabSummary_long.variable_name, tgTabSummary_long.value))
    nUniqueIntsint = int(tgTabWide[["Intsint"]].nunique().iloc[0])
    duration = float(parselmouth.praat.call(soundObj, "Get total duration"))

    additionalInfo = {
        "Duration": duration,
        "PitchKey": key,
        "PitchRange": pitchrange,
        "PitchMean": pitchmean,
        "IntsIntLabels": nintsint,
        "UniqueIntsInt": nUniqueIntsint,
        "IntsIntConcentration": float(nintsint / duration),
        "OptimizationRangeLow": optim_r_low,
        "OptimizationRangeHigh": optim_r_high,
        "OptimizationStep": optim_r_step,
        "OptimizationMidLow": optim_m_low,
        "OptimizationMidHigh": optim_m_high,
        "OptimizationMidStep": optim_m_step
    }

    tgSummary.update(additionalInfo)

    return pd.Series(tgSummary)

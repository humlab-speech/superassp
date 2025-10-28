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

__all__ = ['prosody_measures', 'prosody_measures_from_sound', 'batch_process', '__version__']


def prosody_measures_from_sound(sound, minF=60, maxF=750, windowShift=10):
    """
    Compute prosody measures from parselmouth Sound object.

    This is a wrapper function that accepts a parselmouth.Sound object directly,
    enabling in-memory processing without temporary files. The Sound object
    replaces the 'sound' variable that would normally be loaded from a file path
    in the main prosody_measures() function.

    Args:
        sound: parselmouth.Sound object (not a file path)
        minF: Minimum F0 for pitch extraction in Hz (default: 60)
        maxF: Maximum F0 for pitch extraction in Hz (default: 750)
        windowShift: Window shift for intensity in milliseconds (default: 10)

    Returns:
        pandas.Series with 193 prosodic features, or None if audio < 1 second

    Example:
        >>> import parselmouth
        >>> import numpy as np
        >>> # Create Sound from numpy array (or use av_to_parselmouth_sound in R)
        >>> audio_array = np.random.randn(16000)  # 1 second at 16kHz
        >>> sound = parselmouth.Sound(audio_array, sampling_frequency=16000)
        >>> features = prosody_measures_from_sound(sound)
    """
    import parselmouth
    import pandas as pd
    from . import dysprosody

    # Get the internal prosody measurement function
    # (this requires modifying the dysprosody.py module slightly)
    # For now, use the existing prosody_measures but with Sound object support

    # Check duration
    duration = float(parselmouth.praat.call(sound, "Get total duration"))
    if duration < 1.0:
        print(f"Skipping utterance with less than 1 second duration")
        return None

    # Import the core processing from dysprosody module
    # We need to call the internal functions directly
    from .dysprosody import (
        automatic_min_max_fo,
        momel,
        code_with_intsint,
        spectral_tilt
    )
    import io
    import numpy as np
    from scipy import stats

    windowShift_sec = float(windowShift / 1000)  # Convert to seconds

    # Process the Sound object (same as prosody_measures but without file loading)
    pitchObj, min_fo, max_fo = automatic_min_max_fo(sound, maximum_pitch_span=1.5,
                                                     minimum_fo=minF, maximum_fo=maxF)
    momelPitchTier, momel_pitch_values = momel(pitchObj, minimum_fo=min_fo, maximum_fo=max_fo)

    textGrid, nintsint, pitchmean, optim_r_low, optim_r_high, optim_r_step, \
        optim_m_low, optim_m_high, optim_m_step, pitchrange, key = \
        code_with_intsint(sound, momelPitchTier)

    formantObj = parselmouth.praat.call(sound, "To Formant (burg)", 0.001, 5, 5000, 0.025, 50)
    intensityObj = parselmouth.praat.call(sound, "To Intensity", min_fo, windowShift_sec, "yes")

    # Extract features from TextGrid
    tgTable = pd.read_table(io.StringIO(parselmouth.praat.call(textGrid, "List", False, 20, True, False)))
    tgTabWide = pd.pivot(tgTable[['tmin', 'tier', 'text']], index="tmin", columns="tier", values="text")

    tgTabWide['Intensity'] = tgTabWide.apply(
        lambda x: parselmouth.praat.call(intensityObj, "Get value at time", x.name, "cubic"),
        axis=1
    )

    # Merge in momel pitch values
    momelPitch = pd.DataFrame(
        [(float(a)/1000.0, float(b)) for a, b in momel_pitch_values],
        columns=["time", "momel_pitch"],
        index=pd.Index([float(a)/1000 for a, b in momel_pitch_values], name="tmin")
    )
    tgTabWide = pd.merge_asof(tgTabWide, momelPitch, left_index=True, right_index=True)
    tgTabWide.dropna(axis=0, inplace=True)

    # Compute spectral tilt
    specTilt = tgTabWide.apply(
        lambda x: spectral_tilt(sound, float(x['momel_pitch']), formantObj, x.name),
        axis=1, result_type="expand"
    )
    tgTabWide = tgTabWide.join(specTilt)

    # Compute inter-label differences
    tgTabWide_toDiff = tgTabWide[tgTabWide.columns[3:18]]
    tgTabDiff = pd.DataFrame(np.diff(tgTabWide_toDiff, axis=0), columns=tgTabWide_toDiff.columns)
    tgTabDiff = pd.concat([
        pd.DataFrame(np.zeros((1, tgTabWide_toDiff.shape[1])), columns=tgTabWide_toDiff.columns),
        tgTabDiff
    ], ignore_index=True)

    tgTabDiff.columns = tgTabDiff.columns + "_diff"
    tgTabDiff.index = tgTabWide_toDiff.index
    tgWideTab = tgTabDiff.join(tgTabWide_toDiff)

    # Compute statistics
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
            results = {'tstd': None, 'tmean': None, 'variation': None,
                      'iqr': None, 'tmax': None, 'tmin': None}
        return pd.Series(results)

    tgTabSummary = tgWideTab.apply(lambda x: safe_statistics(x), axis=0, result_type='expand')
    tgTabSummary.columns = ['_'.join(col).strip() for col in tgTabSummary.columns.values]

    tgTabSummary.reset_index(names="statistic", inplace=True)
    tgTabSummary_long = tgTabSummary.melt(id_vars="statistic")
    tgTabSummary_long["variable_name"] = (tgTabSummary_long.variable + "_" +
                                          tgTabSummary_long.statistic)

    tgTabSummary.rename(
        index={"tstd": "std", "tmean": "mean", "hdmedian": "median",
               "variation": "var", "iqr": "iqr", "tmax": "max", "tmin": "min"},
        columns={"momel_pitch": "momelpitch"},
        inplace=True
    )
    tgTabSummary.columns = tgTabSummary.columns.to_flat_index() + "_diff"

    # Prepare results
    tgSummary = dict(zip(tgTabSummary_long.variable_name, tgTabSummary_long.value))
    nUniqueIntsint = int(tgTabWide[["Intsint"]].nunique().iloc[0])

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

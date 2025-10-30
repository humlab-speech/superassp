import parselmouth as pm
import pandas as pd
import numpy as np
import io


def praat_intensity(
    sound,  # Changed from soundFile - now accepts Sound object
    beginTime=0.0,
    endTime=0.0,
    time_step=0.0,
    minimal_f0_frequency=50.0,
    subtract_mean=True,
    windowShape=pm.WindowShape.GAUSSIAN1,
    relativeWidth=1.0):
    """
    Compute the intensity contour of a sound signal.

    This function replicates the functionality of intensity.praat, computing the
    intensity (loudness) contour of a sound file using Praat's intensity algorithm.

    Parameters:
    -----------
    sound : parselmouth.Sound
        Sound object (in-memory, not a file path)
    beginTime : float
        Start time for analysis (0.0 for full file)
    endTime : float
        End time for analysis (0.0 for full file)
    time_step : float
        Time step between frames (in seconds). If 0.0, Praat uses automatic calculation.
    minimal_f0_frequency : float
        Minimum pitch frequency to consider (Hz). This determines the time resolution.
        A smaller value gives more temporal detail but less time resolution.
    subtract_mean : bool
        Whether to subtract the mean intensity (normalize)
    windowShape : pm.WindowShape
        Window shape for time windowing when extracting part
    relativeWidth : float
        Relative width for windowing when extracting part

    Returns:
    --------
    pd.DataFrame
        DataFrame with columns: Time(s), Intensity(dB)
    """

    # Accept Sound object directly (in-memory)
    snd = sound
    dur = snd.get_total_duration()

    # Handle time windowing (if needed - usually already windowed by R)
    if beginTime > 0.0 and endTime > 0.0 and beginTime >= 0.0 and endTime <= dur:
        snd = snd.extract_part(beginTime, endTime, windowShape, relativeWidth, True)

    # Compute intensity
    # Note: time_step of 0.0 means automatic calculation (typically 0.8 / minimal_f0_frequency)
    intensity = pm.praat.call(snd, "To Intensity", minimal_f0_frequency, time_step, subtract_mean)

    # Convert to IntensityTier
    intensity_tier = pm.praat.call(intensity, "Down to IntensityTier")

    # Convert to TableOfReal
    table_of_real = pm.praat.call(intensity_tier, "Down to TableOfReal")

    # Convert to Table
    table = pm.praat.call(table_of_real, "To Table", "dummy")

    # Remove the dummy column
    pm.praat.call(table, "Remove column", "dummy")

    # Convert to pandas DataFrame
    result = pd.read_table(io.StringIO(pm.praat.call(table, "List", True)))

    return result

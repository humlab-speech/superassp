import parselmouth as pm
import pandas as pd
import numpy as np
import io


def praat_spectral_moments(
    soundFile,
    beginTime=0.0,
    endTime=0.0,
    windowLength=0.005,
    maximum_frequency=0.0,
    time_step=0.005,
    frequency_step=20.0,
    power=2.0,
    windowShape=pm.WindowShape.GAUSSIAN1,
    relativeWidth=1.0):
    """
    Compute spectral moments (1-4) from a sound signal.

    This function replicates the functionality of praat_spectral_moments.praat,
    computing the center of gravity, standard deviation, skewness, and kurtosis
    of the spectrum at each time frame.

    The spectral moments characterize the shape of the spectrum:
    - Center of Gravity (Moment 1): Average frequency weighted by amplitude
    - Standard Deviation (Moment 2): Spread of energy across frequencies
    - Skewness (Moment 3): Asymmetry of the spectral distribution
    - Kurtosis (Moment 4): "Peakedness" of the spectral distribution

    Parameters:
    -----------
    soundFile : str
        Path to the sound file
    beginTime : float
        Start time for analysis (0.0 for full file)
    endTime : float
        End time for analysis (0.0 for full file)
    windowLength : float
        Analysis window length (in seconds)
    maximum_frequency : float
        Maximum frequency to analyze (Hz). If 0.0, uses Nyquist frequency (sr/2)
    time_step : float
        Time step between frames (in seconds)
    frequency_step : float
        Frequency resolution (in Hz)
    power : float
        Power for moment calculation (typically 2)
    windowShape : pm.WindowShape
        Window shape for time windowing when extracting part
    relativeWidth : float
        Relative width for windowing when extracting part

    Returns:
    --------
    pd.DataFrame
        DataFrame with columns: Time, CenterOfGravity, SD, Skewness, Kurtosis
    """

    # Load sound
    snd = pm.Sound(soundFile)
    dur = snd.get_total_duration()
    sr = snd.get_sampling_frequency()

    # Set maximum frequency to Nyquist if not specified
    if maximum_frequency == 0.0:
        maximum_frequency = sr / 2.0

    # Handle time windowing
    if beginTime > 0.0 or endTime > 0.0:
        if beginTime >= 0.0 and endTime <= dur:
            snd = snd.extract_part(beginTime, endTime, windowShape, relativeWidth, True)

    # Create output table
    out_table = pm.praat.call("Create Table with column names", "outTable", 0,
                               "Time CenterOfGravity SD Skewness Kurtosis")

    # Create spectrogram
    spectrogram = pm.praat.call(snd, "To Spectrogram", windowLength, maximum_frequency,
                                 time_step, frequency_step, "Gaussian")

    # Get number of frames
    num_frames = pm.praat.call(spectrogram, "Get number of frames")

    # Process each frame
    for frame in range(1, num_frames + 1):
        # Get time for this frame
        curr_time = pm.praat.call(spectrogram, "Get time from frame number", frame)

        # Extract spectrum slice at this time
        spectrum = pm.praat.call(spectrogram, "To Spectrum (slice)", curr_time)

        # Compute spectral moments
        cog = pm.praat.call(spectrum, "Get centre of gravity", power)
        sd = pm.praat.call(spectrum, "Get standard deviation", power)
        skew = pm.praat.call(spectrum, "Get skewness", power)
        kurt = pm.praat.call(spectrum, "Get kurtosis", power)

        # Clean up spectrum object
        pm.praat.call(spectrum, "Remove")

        # Append row to output table
        pm.praat.call(out_table, "Append row")
        row = pm.praat.call(out_table, "Get number of rows")
        pm.praat.call(out_table, "Set numeric value", row, "Time", curr_time)
        pm.praat.call(out_table, "Set numeric value", row, "CenterOfGravity", cog)
        pm.praat.call(out_table, "Set numeric value", row, "SD", sd)
        pm.praat.call(out_table, "Set numeric value", row, "Skewness", skew)
        pm.praat.call(out_table, "Set numeric value", row, "Kurtosis", kurt)

    # Convert to pandas DataFrame
    result = pd.read_table(io.StringIO(pm.praat.call(out_table, "List", True)))

    # Clean up
    pm.praat.call(out_table, "Remove")

    return result

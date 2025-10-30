import parselmouth as pm
import pandas as pd
import numpy as np
import io
import math


def praat_formantpath_burg(
    sound,  # Changed from soundFile - now accepts Sound object
    beginTime=0.0,
    endTime=0.0,
    time_step=0.005,
    number_of_formants=5.0,
    maxHzFormant=5500.0,
    windowLength=0.025,
    pre_emphasis=50.0,
    ceiling_step_size=0.05,
    number_of_steps_each_direction=4,
    track_formants=True,
    number_of_tracks=3,
    reference_F1=550,
    reference_F2=1650,
    reference_F3=2750,
    reference_F4=3850,
    reference_F5=4950,
    frequency_cost=1.0,
    bandwidth_cost=1.0,
    transition_cost=1.0,
    windowShape=pm.WindowShape.GAUSSIAN1,
    relativeWidth=1.0,
    spectrogram_window_shape="Gaussian",
    spectrogram_resolution=40.0):
    """
    Compute formant tracks using Praat's FormantPath (Burg) algorithm.

    This function replicates the functionality of formantpath_burg.praat, which uses
    the FormantPath object to automatically find the optimal formant ceiling and
    track formants across time. This is a more robust method than simple Burg analysis.

    The FormantPath algorithm:
    1. Computes formants at multiple ceiling frequencies
    2. Optionally tracks formants to find the optimal path
    3. Extracts formant frequencies, bandwidths, and amplitudes
    4. Computes formant intensities from the spectrogram

    Parameters:
    -----------
    sound : parselmouth.Sound
        Sound object (in-memory, not a file path)
    beginTime : float
        Start time for analysis (0.0 for full file)
    endTime : float
        End time for analysis (0.0 for full file)
    time_step : float
        Time step between frames (in seconds)
    number_of_formants : float
        Number of formants to track
    maxHzFormant : float
        Maximum formant frequency (Hz)
    windowLength : float
        Analysis window length (in seconds)
    pre_emphasis : float
        Pre-emphasis frequency (Hz)
    ceiling_step_size : float
        Step size for ceiling variation
    number_of_steps_each_direction : int
        Number of ceiling steps in each direction
    track_formants : bool
        Whether to use formant tracking
    number_of_tracks : int
        Number of formant tracks to extract
    reference_F1 : int
        Reference frequency for F1 (Hz)
    reference_F2 : int
        Reference frequency for F2 (Hz)
    reference_F3 : int
        Reference frequency for F3 (Hz)
    reference_F4 : int
        Reference frequency for F4 (Hz)
    reference_F5 : int
        Reference frequency for F5 (Hz)
    frequency_cost : float
        Weight for frequency deviation in tracking
    bandwidth_cost : float
        Weight for bandwidth in tracking
    transition_cost : float
        Weight for formant transitions in tracking
    windowShape : pm.WindowShape
        Window shape for time windowing when extracting part
    relativeWidth : float
        Relative width for windowing when extracting part
    spectrogram_window_shape : pm.SpectralAnalysisWindowShape
        Window shape for spectrogram analysis
    spectrogram_resolution : float
        Frequency resolution for spectrogram (Hz)

    Returns:
    --------
    pd.DataFrame
        DataFrame with columns: time(s), F1(Hz), B1(Hz), F2(Hz), B2(Hz), ..., L1(dB), L2(dB), ...
        where Fn is formant frequency, Bn is bandwidth, and Ln is formant amplitude
    """

    # Accept Sound object directly (in-memory)
    snd = sound
    dur = snd.get_total_duration()

    # Handle time windowing (if needed - usually already windowed by R)
    if beginTime > 0.0 and endTime > 0.0 and beginTime >= 0.0 and endTime <= dur:
        snd = snd.extract_part(beginTime, endTime, windowShape, relativeWidth, True)

    # Create FormantPath using Burg algorithm
    formant_path = pm.praat.call(snd, "To FormantPath (burg)", time_step, number_of_formants,
                                  maxHzFormant, windowLength, pre_emphasis, ceiling_step_size,
                                  number_of_steps_each_direction)

    # Extract Formant object
    formant = pm.praat.call(formant_path, "Extract Formant")

    # Track formants if requested
    if track_formants:
        formant = pm.praat.call(formant, "Track", number_of_tracks, reference_F1, reference_F2,
                                reference_F3, reference_F4, reference_F5, frequency_cost,
                                bandwidth_cost, transition_cost)
        number_of_formants = min(number_of_formants, number_of_tracks)

    # Convert to table
    # Parameters: includeFrameNumbers, includeTimes, timeDecimals, includeIntensity,
    #             intensityDecimals, includeNumFormants, numFormantsDecimals, includeFrequency
    formant_table = pm.praat.call(formant, "Down to Table", True, True, 10, True, 3, True, 3, True)

    # Get number of rows
    num_rows = pm.praat.call(formant_table, "Get number of rows")

    # Add columns for formant amplitudes (L1, L2, etc.)
    for f in range(1, math.ceil(number_of_formants) + 1):
        lab = f"L{f}(dB)"
        pm.praat.call(formant_table, "Append column", lab)

    # Create spectrogram for amplitude extraction
    # Add extra frequency range to avoid edge effects
    maxHz = maxHzFormant + 2000.0
    # Use string name for window shape (Praat will convert)
    spec_window = spectrogram_window_shape if isinstance(spectrogram_window_shape, str) else "Gaussian"
    spectrogram = pm.praat.call(snd, "To Spectrogram", windowLength, maxHz, time_step,
                                 spectrogram_resolution, spec_window)

    # Extract formant amplitudes from spectrogram
    for r in range(1, num_rows + 1):
        for f in range(1, math.ceil(number_of_formants) + 1):
            out_lab = f"L{f}(dB)"
            in_lab = f"F{f}(Hz)"

            # Get formant frequency
            curr_freq = pm.praat.call(formant_table, "Get value", r, in_lab)

            if curr_freq != "--undefined--":
                # Get time for this frame
                curr_time = pm.praat.call(formant_table, "Get value", r, "time(s)")

                # Get power at this time and frequency from spectrogram
                curr_intensity = pm.praat.call(spectrogram, "Get power at", float(curr_time),
                                               float(curr_freq))
            else:
                curr_intensity = None  # Will be represented as NaN in pandas

            # Set the intensity value
            pm.praat.call(formant_table, "Set numeric value", r, out_lab, curr_intensity)

    # Convert to pandas DataFrame
    result = pd.read_table(io.StringIO(pm.praat.call(formant_table, "List", True)))

    # Clean up
    pm.praat.call(formant_table, "Remove")

    return result

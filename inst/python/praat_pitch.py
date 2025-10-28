import parselmouth as pm
import pandas as pd
import numpy as np
import io


def praat_pitch(
    soundFile,
    beginTime=0.0,
    endTime=0.0,
    time_step=0.005,
    window_length=0.040,
    minimum_f0=75.0,
    maximum_f0=600.0,
    very_accurate=True,
    number_of_candidates=15,
    silence_threshold=0.03,
    voicing_threshold=0.45,
    octave_cost=0.01,
    octave_jump_cost=0.35,
    voiced_voiceless_cost=0.14,
    only_correlation_methods=True,
    minimum_filter_frequency=70.0,
    maximum_filter_frequency=5000.0,
    number_of_filters=250,
    maximum_frequency_components=1250.0,
    maximum_number_of_subharmonics=15,
    compression_factor=0.84,
    number_of_points_per_octave=48,
    windowShape=pm.WindowShape.GAUSSIAN1,
    relativeWidth=1.0):
    """
    Compute F0 using multiple pitch tracking algorithms (autocorrelation and cross-correlation).

    This function replicates the functionality of praat_pitch.praat, computing pitch tracks
    using multiple methods: autocorrelation (ac), cross-correlation (cc), and optionally
    SPINET and SHS methods.

    Parameters:
    -----------
    soundFile : str
        Path to the sound file
    beginTime : float
        Start time for analysis (0.0 for full file)
    endTime : float
        End time for analysis (0.0 for full file)
    time_step : float
        Time step between frames (in seconds)
    window_length : float
        Analysis window length (in seconds)
    minimum_f0 : float
        Minimum pitch in Hz
    maximum_f0 : float
        Maximum pitch in Hz
    very_accurate : bool
        Use very accurate pitch tracking
    number_of_candidates : int
        Number of pitch candidates to consider
    silence_threshold : float
        Silence threshold
    voicing_threshold : float
        Voicing threshold
    octave_cost : float
        Cost for octave jumps
    octave_jump_cost : float
        Cost for octave jumps
    voiced_voiceless_cost : float
        Cost for voiced/voiceless transitions
    only_correlation_methods : bool
        If True, only use cc and ac methods; if False, also use SPINET and SHS
    minimum_filter_frequency : float
        Minimum filter frequency for SPINET (Hz)
    maximum_filter_frequency : float
        Maximum filter frequency for SPINET (Hz)
    number_of_filters : int
        Number of filters for SPINET
    maximum_frequency_components : float
        Maximum frequency components for SHS (Hz)
    maximum_number_of_subharmonics : int
        Maximum number of subharmonics for SHS
    compression_factor : float
        Compression factor for SHS
    number_of_points_per_octave : int
        Number of points per octave for SHS
    windowShape : pm.WindowShape
        Window shape for time windowing
    relativeWidth : float
        Relative width for windowing

    Returns:
    --------
    pd.DataFrame
        DataFrame with columns: time, cc, ac, and optionally spinet and shs
    """

    # Load sound
    snd = pm.Sound(soundFile)
    dur = snd.get_total_duration()

    # Handle time windowing
    if beginTime > 0.0 or endTime > 0.0:
        if beginTime >= 0.0 and endTime <= dur:
            snd = snd.extract_part(beginTime, endTime, windowShape, relativeWidth, True)

    # Cross-correlation method
    pitch_cc = pm.praat.call(snd, "To Pitch (cc)", time_step, minimum_f0, number_of_candidates,
                              very_accurate, silence_threshold, voicing_threshold, octave_cost,
                              octave_jump_cost, voiced_voiceless_cost, maximum_f0)
    cc_matrix = pm.praat.call(pitch_cc, "To Matrix")

    # Autocorrelation method
    pitch_ac = pm.praat.call(snd, "To Pitch (ac)", time_step, minimum_f0, number_of_candidates,
                              very_accurate, silence_threshold, voicing_threshold, octave_cost,
                              octave_jump_cost, voiced_voiceless_cost, maximum_f0)
    ac_matrix = pm.praat.call(pitch_ac, "To Matrix")

    # Combine matrices
    matrices = [cc_matrix, ac_matrix]
    matrix_names = ["cc", "ac"]

    # Optional: SPINET and SHS methods
    if not only_correlation_methods:
        # SPINET method
        pitch_spinet = pm.praat.call(snd, "To Pitch (SPINET)", time_step, window_length,
                                      minimum_filter_frequency, maximum_filter_frequency,
                                      number_of_filters, maximum_f0, number_of_candidates)
        spinet_matrix = pm.praat.call(pitch_spinet, "To Matrix")
        matrices.append(spinet_matrix)
        matrix_names.append("spinet")

        # SHS method
        pitch_shs = pm.praat.call(snd, "To Pitch (shs)", time_step, minimum_f0,
                                   number_of_candidates, maximum_frequency_components,
                                   maximum_number_of_subharmonics, compression_factor,
                                   maximum_f0, number_of_points_per_octave)
        shs_matrix = pm.praat.call(pitch_shs, "To Matrix")
        matrices.append(shs_matrix)
        matrix_names.append("shs")

    # Merge matrices
    merged_matrix = matrices[0]
    for mat in matrices[1:]:
        merged_matrix = pm.praat.call([merged_matrix, mat], "Merge (append rows)")

    # Transpose and convert to table
    transposed = pm.praat.call(merged_matrix, "Transpose")
    table_of_real = pm.praat.call(transposed, "To TableOfReal")
    table = pm.praat.call(table_of_real, "To Table", "time")

    # Set column labels
    pm.praat.call(table, "Set column label (index)", 2, "cc")
    pm.praat.call(table, "Set column label (index)", 3, "ac")

    if not only_correlation_methods:
        pm.praat.call(table, "Set column label (index)", 4, "spinet")
        pm.praat.call(table, "Set column label (index)", 5, "shs")

    # Convert to pandas DataFrame
    result = pd.read_table(io.StringIO(pm.praat.call(table, "List", True)))

    return result


def praat_pitch_from_sound(
    sound,
    time_step=0.005,
    window_length=0.040,
    minimum_f0=75.0,
    maximum_f0=600.0,
    very_accurate=True,
    number_of_candidates=15,
    silence_threshold=0.03,
    voicing_threshold=0.45,
    octave_cost=0.01,
    octave_jump_cost=0.35,
    voiced_voiceless_cost=0.14,
    only_correlation_methods=True,
    minimum_filter_frequency=70.0,
    maximum_filter_frequency=5000.0,
    number_of_filters=250,
    maximum_frequency_components=1250.0,
    maximum_number_of_subharmonics=15,
    compression_factor=0.84,
    number_of_points_per_octave=48,
    windowShape=pm.WindowShape.GAUSSIAN1,
    relativeWidth=1.0):
    """
    Compute F0 using multiple pitch tracking algorithms from parselmouth Sound object.

    This function performs the same analysis as praat_pitch() but accepts a
    parselmouth.Sound object directly, enabling in-memory processing without
    file I/O.

    Parameters:
    -----------
    sound : parselmouth.Sound
        Sound object (not a file path)
    time_step : float
        Time step between frames (in seconds)
    window_length : float
        Analysis window length (in seconds)
    minimum_f0 : float
        Minimum pitch in Hz
    maximum_f0 : float
        Maximum pitch in Hz
    very_accurate : bool
        Use very accurate pitch tracking
    number_of_candidates : int
        Number of pitch candidates to consider
    silence_threshold : float
        Silence threshold
    voicing_threshold : float
        Voicing threshold
    octave_cost : float
        Cost for octave jumps
    octave_jump_cost : float
        Cost for octave jumps
    voiced_voiceless_cost : float
        Cost for voiced/voiceless transitions
    only_correlation_methods : bool
        If True, only use cc and ac methods; if False, also use SPINET and SHS
    minimum_filter_frequency : float
        Minimum filter frequency for SPINET (Hz)
    maximum_filter_frequency : float
        Maximum filter frequency for SPINET (Hz)
    number_of_filters : int
        Number of filters for SPINET
    maximum_frequency_components : float
        Maximum frequency components for SHS (Hz)
    maximum_number_of_subharmonics : int
        Maximum number of subharmonics for SHS
    compression_factor : float
        Compression factor for SHS
    number_of_points_per_octave : int
        Number of points per octave for SHS
    windowShape : pm.WindowShape
        Window shape for time windowing
    relativeWidth : float
        Relative width for windowing

    Returns:
    --------
    pd.DataFrame
        DataFrame with columns: time, cc, ac, and optionally spinet and shs
    """

    # sound is already a parselmouth.Sound object - use it directly!
    snd = sound

    # Cross-correlation method
    pitch_cc = pm.praat.call(snd, "To Pitch (cc)", time_step, minimum_f0, number_of_candidates,
                              very_accurate, silence_threshold, voicing_threshold, octave_cost,
                              octave_jump_cost, voiced_voiceless_cost, maximum_f0)
    cc_matrix = pm.praat.call(pitch_cc, "To Matrix")

    # Autocorrelation method
    pitch_ac = pm.praat.call(snd, "To Pitch (ac)", time_step, minimum_f0, number_of_candidates,
                              very_accurate, silence_threshold, voicing_threshold, octave_cost,
                              octave_jump_cost, voiced_voiceless_cost, maximum_f0)
    ac_matrix = pm.praat.call(pitch_ac, "To Matrix")

    # Combine matrices
    matrices = [cc_matrix, ac_matrix]
    matrix_names = ["cc", "ac"]

    # Optional: SPINET and SHS methods
    if not only_correlation_methods:
        # SPINET method
        pitch_spinet = pm.praat.call(snd, "To Pitch (SPINET)", time_step, window_length,
                                      minimum_filter_frequency, maximum_filter_frequency,
                                      number_of_filters, maximum_f0, number_of_candidates)
        spinet_matrix = pm.praat.call(pitch_spinet, "To Matrix")
        matrices.append(spinet_matrix)
        matrix_names.append("spinet")

        # SHS method
        pitch_shs = pm.praat.call(snd, "To Pitch (shs)", time_step, minimum_f0,
                                   number_of_candidates, maximum_frequency_components,
                                   maximum_number_of_subharmonics, compression_factor,
                                   maximum_f0, number_of_points_per_octave)
        shs_matrix = pm.praat.call(pitch_shs, "To Matrix")
        matrices.append(shs_matrix)
        matrix_names.append("shs")

    # Merge matrices
    merged_matrix = matrices[0]
    for mat in matrices[1:]:
        merged_matrix = pm.praat.call([merged_matrix, mat], "Merge (append rows)")

    # Transpose and convert to table
    transposed = pm.praat.call(merged_matrix, "Transpose")
    table_of_real = pm.praat.call(transposed, "To TableOfReal")
    table = pm.praat.call(table_of_real, "To Table", "time")

    # Set column labels
    pm.praat.call(table, "Set column label (index)", 2, "cc")
    pm.praat.call(table, "Set column label (index)", 3, "ac")

    if not only_correlation_methods:
        pm.praat.call(table, "Set column label (index)", 4, "spinet")
        pm.praat.call(table, "Set column label (index)", 5, "shs")

    # Convert to pandas DataFrame
    result = pd.read_table(io.StringIO(pm.praat.call(table, "List", True)))

    return result


# Keep original function name for backwards compatibility
trk_pitchp = praat_pitch
trk_pitchp_from_sound = praat_pitch_from_sound

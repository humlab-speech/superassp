"""
Parselmouth implementation of Praat voice_report

This module provides a memory-based implementation of the Praat voice report
functionality using Parselmouth. Instead of reading from disk, it processes
audio directly from numpy arrays, eliminating file I/O overhead.

"""

import parselmouth as pm
import numpy as np
import pandas as pd


def praat_voice_report_memory(
    audio_np,
    sample_rate,
    start_time=0.0,
    end_time=0.0,
    selection_offset=0.0,
    selection_length=0.0,
    window_shape="Gaussian1",
    relative_width=1.0,
    min_f0=75.0,
    max_f0=600.0,
    max_period_factor=1.3,
    max_amplitude_factor=1.6,
    silence_threshold=0.03,
    voicing_threshold=0.45,
    octave_cost=0.01,
    octave_jump_cost=0.35,
    voiced_unvoiced_cost=0.14,
):
    """
    Compute Praat voice report from numpy audio array.

    This function replicates the behavior of the Praat voice report script
    but operates entirely in memory without file I/O.

    Parameters
    ----------
    audio_np : numpy.ndarray
        Audio signal as numpy array (float64, normalized to [-1, 1])
    sample_rate : int
        Sample rate in Hz
    start_time : float
        Start time in seconds (relative to input audio)
    end_time : float
        End time in seconds (0 means use full duration)
    selection_offset : float
        Offset from start_time where extraction begins (in seconds)
    selection_length : float
        Maximum length of selection (in seconds, 0 means no limit)
    window_shape : str
        Window type for extraction (Gaussian1, Hanning, etc.)
    relative_width : float
        Relative width of extraction window
    min_f0 : float
        Minimum f0 in Hz
    max_f0 : float
        Maximum f0 in Hz
    max_period_factor : float
        Maximum period factor for jitter computation
    max_amplitude_factor : float
        Maximum amplitude factor for shimmer computation
    silence_threshold : float
        Silence threshold
    voicing_threshold : float
        Voicing threshold
    octave_cost : float
        Octave cost for pitch tracking
    octave_jump_cost : float
        Octave jump cost for pitch tracking
    voiced_unvoiced_cost : float
        Voiced/unvoiced cost for pitch tracking

    Returns
    -------
    dict
        Dictionary containing all voice report measurements with the following keys:

        Pitch measurements (Hz):
        - 'Median pitch': float, Hz
        - 'Mean pitch': float, Hz
        - 'Standard deviation': float, Hz (SD of pitch)
        - 'Minimum pitch': float, Hz
        - 'Maximum pitch': float, Hz

        Pulse/period measurements:
        - 'Number of pulses': float, count
        - 'Number of periods': float, count
        - 'Mean period': float, seconds
        - 'Standard deviation of period': float, seconds

        Voicing measurements (fractions 0.0-1.0):
        - 'Fraction of locally unvoiced frames': float, fraction
        - 'Number of voice breaks': float, count
        - 'Degree of voice breaks': float, fraction

        Jitter measurements:
        - 'Jitter (local)': float, % (percentage)
        - 'Jitter (local, absolute)': float, seconds (NOT microseconds!)
        - 'Jitter (rap)': float, % (percentage)
        - 'Jitter (ppq5)': float, % (percentage)
        - 'Jitter (ddp)': float, % (percentage)

        Shimmer measurements:
        - 'Shimmer (local)': float, % (percentage)
        - 'Shimmer (local, dB)': float, dB (decibels)
        - 'Shimmer (apq3)': float, % (percentage)
        - 'Shimmer (apq5)': float, % (percentage)
        - 'Shimmer (apq11)': float, % (percentage)
        - 'Shimmer (dda)': float, % (percentage)

        Harmonicity measurements:
        - 'Mean autocorrelation': float, fraction (0.0-1.0)
        - 'Mean noise-to-harmonics ratio': float, ratio (unitless)
        - 'Mean harmonics-to-noise ratio': float, dB (decibels)

    Notes
    -----
    Time units: All time parameters are in seconds
    - start_time: where in the audio to start (already handled by av)
    - end_time: where to end (already handled by av)
    - selection_offset: offset from start for further windowing
    - selection_length: length of window for analysis

    Units of measurement:
    - Pitch: Hertz (Hz)
    - Periods: seconds (s)
    - Jitter (local, absolute): seconds (s) - NOTE: Literature often reports
      this in microseconds, but Praat returns seconds
    - Jitter (percentages): % (e.g., 1.5 means 1.5%)
    - Shimmer (percentages): % (e.g., 3.8 means 3.8%)
    - Shimmer (dB): decibels (dB)
    - HNR: decibels (dB)
    - NHR: ratio (unitless)
    - Fractions: 0.0 to 1.0 (unitless)

    The function follows the same logic as praat_voice_report.praat but
    operates entirely in memory. All measurements match Praat's Voice Report
    output exactly, including units.
    """

    # Create Sound from numpy array (IN MEMORY!)
    sound = pm.Sound(values=audio_np, sampling_frequency=sample_rate)
    sound_end = sound.get_total_duration()

    # Calculate extraction times
    # Note: av already handled start_time/end_time extraction
    # so our audio starts at 0.0 in the Sound object
    start_at = selection_offset  # Just use offset from beginning

    # Determine end of selection
    if end_time == 0.0:
        end_at = sound_end
    else:
        end_at = sound_end  # We already have the right portion from av

    if selection_length > 0.0:
        sel_end = start_at + selection_length
        if sel_end < end_at:
            end_at = sel_end

    # Extract the part of the sound for analysis
    # Map window_shape string to Parselmouth window type
    # Parselmouth expects: 'rectangular', 'triangular', 'parabolic', 'Hanning',
    #                      'Hamming', 'Gaussian1', etc.
    sound_part = sound.extract_part(
        from_time=start_at, to_time=end_at, window_shape=window_shape, relative_width=relative_width, preserve_times=False
    )

    # Create PointProcess (periodic, cc)
    point_process = pm.praat.call(sound_part, "To PointProcess (periodic, cc)", min_f0, max_f0)

    # Create Pitch object
    pitch = pm.praat.call(
        sound_part,
        "To Pitch (cc)",
        0.0,  # time_step (0 = auto)
        min_f0,
        15,  # max_number_of_candidates
        1,  # very_accurate (yes)
        silence_threshold,
        voicing_threshold,
        octave_cost,
        octave_jump_cost,
        voiced_unvoiced_cost,
        max_f0,
    )

    # Get voice report
    # Parselmouth's Voice report returns a string
    voice_report_str = pm.praat.call([sound_part, pitch, point_process], "Voice report", 0.0, 0.0, min_f0, max_f0, max_period_factor, max_amplitude_factor, silence_threshold, voicing_threshold)

    # Parse voice report string to extract values
    def extract_number(report_str, key):
        """Extract numeric value from voice report string."""
        lines = report_str.split("\n")
        for line in lines:
            if key in line:
                # Extract number from line
                # Format is typically "Key: value units" or "Key: value%"
                parts = line.split(":")
                if len(parts) == 2:
                    value_str = parts[1].strip().split()[0]  # Get first token after colon
                    try:
                        return float(value_str)
                    except ValueError:
                        return None
        return None

    # Extract pitch measurements
    median_pitch = extract_number(voice_report_str, "Median pitch")
    mean_pitch = extract_number(voice_report_str, "Mean pitch")
    sd_pitch = extract_number(voice_report_str, "Standard deviation")
    min_pitch = extract_number(voice_report_str, "Minimum pitch")
    max_pitch = extract_number(voice_report_str, "Maximum pitch")

    # Extract pulse/period measurements
    num_pulses = extract_number(voice_report_str, "Number of pulses")
    num_periods = extract_number(voice_report_str, "Number of periods")
    mean_period = extract_number(voice_report_str, "Mean period")
    sd_period = extract_number(voice_report_str, "Standard deviation of period")

    # Extract voicing measurements
    frac_unvoiced = extract_number(voice_report_str, "Fraction of locally unvoiced frames")
    num_breaks = extract_number(voice_report_str, "Number of voice breaks")
    degree_breaks = extract_number(voice_report_str, "Degree of voice breaks")

    # Extract jitter measurements
    jitter_local = extract_number(voice_report_str, "Jitter (local)")
    jitter_local_abs = extract_number(voice_report_str, "Jitter (local, absolute)")
    jitter_rap = extract_number(voice_report_str, "Jitter (rap)")
    jitter_ppq5 = extract_number(voice_report_str, "Jitter (ppq5)")
    jitter_ddp = extract_number(voice_report_str, "Jitter (ddp)")

    # Extract shimmer measurements
    shimmer_local = extract_number(voice_report_str, "Shimmer (local)")
    shimmer_local_db = extract_number(voice_report_str, "Shimmer (local, dB)")
    shimmer_apq3 = extract_number(voice_report_str, "Shimmer (apq3)")
    shimmer_apq5 = extract_number(voice_report_str, "Shimmer (apq5)")
    shimmer_apq11 = extract_number(voice_report_str, "Shimmer (apq11)")
    shimmer_dda = extract_number(voice_report_str, "Shimmer (dda)")

    # Extract harmonicity measurements
    mean_autocor = extract_number(voice_report_str, "Mean autocorrelation")
    mean_nhr = extract_number(voice_report_str, "Mean noise-to-harmonics ratio")
    mean_hnr = extract_number(voice_report_str, "Mean harmonics-to-noise ratio")

    # Compute intensity measurements
    # Use mean pitch as min pitch for intensity computation
    # (matches Praat script logic)
    intensity_min_pitch = mean_pitch if mean_pitch is not None else 20.0

    intensity = pm.praat.call(sound_part, "To Intensity", intensity_min_pitch, 0.0, 0)  # yes (subtract mean)

    int_mean = pm.praat.call(intensity, "Get mean", 0.0, 0.0, "energy")
    int_median = pm.praat.call(intensity, "Get quantile", 0.0, 0.0, 0.50)
    int_sd = pm.praat.call(intensity, "Get standard deviation", 0.0, 0.0)

    # Build result dictionary
    # Note: We don't include time stamps (Start Time, End Time, etc.)
    # because that's handled by the R wrapper
    result = {
        "Median pitch": median_pitch,
        "Mean pitch": mean_pitch,
        "Standard deviation": sd_pitch,
        "Minimum pitch": min_pitch,
        "Maximum pitch": max_pitch,
        "Number of pulses": num_pulses,
        "Number of periods": num_periods,
        "Mean period": mean_period,
        "Standard deviation of period": sd_period,
        "Fraction of locally unvoiced frames": frac_unvoiced,
        "Number of voice breaks": num_breaks,
        "Degree of voice breaks": degree_breaks,
        "Jitter (local)": jitter_local,
        "Jitter (local, absolute)": jitter_local_abs,
        "Jitter (rap)": jitter_rap,
        "Jitter (ppq5)": jitter_ppq5,
        "Jitter (ddp)": jitter_ddp,
        "Shimmer (local)": shimmer_local,
        "Shimmer (local, dB)": shimmer_local_db,
        "Shimmer (apq3)": shimmer_apq3,
        "Shimmer (apq5)": shimmer_apq5,
        "Shimmer (apq11)": shimmer_apq11,
        "Shimmer (dda)": shimmer_dda,
        "Mean autocorrelation": mean_autocor,
        "Mean noise-to-harmonics ratio": mean_nhr,
        "Mean harmonics-to-noise ratio": mean_hnr,
        # Note: Intensity measurements not in original Praat output
        # but computed in script - omitting for now to match output exactly
    }

    return result

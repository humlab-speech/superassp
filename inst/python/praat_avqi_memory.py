"""
Parselmouth-based AVQI implementation (memory-based, av package integration)

This module provides a memory-based implementation of the Acoustic Voice Quality Index (AVQI)
using Parselmouth. It processes audio directly from numpy arrays loaded by the av package,
eliminating file I/O overhead.

Based on AVQI v3.01 by Youri Maryn and Paul Corthals.
"""

import parselmouth as pm
from parselmouth.praat import call
import numpy as np


def praat_avqi_memory(
    sv_audio_list,  # List of (audio_np, sample_rate) tuples for sustained vowels
    cs_audio_list,  # List of (audio_np, sample_rate) tuples for continuous speech
    speaker_name="",
    speaker_id="",
    speaker_dob="",
    assessment_date="",
):
    """
    Compute AVQI from numpy audio arrays (memory-based, no file I/O).

    This function replicates the AVQI 3.01 Praat script but operates entirely
    in memory using audio loaded by the av package.

    Parameters
    ----------
    sv_audio_list : list of dict
        List of sustained vowel audio segments, each dict contains:
        - 'audio_np': numpy array (float64, mono)
        - 'sample_rate': int, sampling frequency in Hz
    cs_audio_list : list of dict
        List of continuous speech audio segments, each dict contains:
        - 'audio_np': numpy array (float64, mono)
        - 'sample_rate': int, sampling frequency in Hz
    speaker_name : str
        Name of speaker (optional)
    speaker_id : str
        Speaker identifier (optional)
    speaker_dob : str
        Date of birth (optional)
    assessment_date : str
        Assessment date (optional)

    Returns
    -------
    dict
        Dictionary with AVQI measurements:
        - 'AVQI_VERSION': str, version identifier
        - 'Speaker': str, speaker name
        - 'ID': str, speaker ID
        - 'CPPS': float, Smoothed Cepstral Peak Prominence
        - 'HNR': float, Harmonics-to-Noise Ratio (dB)
        - 'Shim_local': float, Local Shimmer (%)
        - 'Shim_local_DB': float, Local Shimmer (dB)
        - 'LTAS_Slope': float, Long-Term Average Spectrum slope (dB)
        - 'LTAS_Tilt': float, LTAS tilt (dB)
        - 'AVQI': float, Acoustic Voice Quality Index score

    Notes
    -----
    The AVQI formula (Barsties & Maryn, 2015):
    AVQI = (4.152 - 0.177*CPPS - 0.006*HNR - 0.037*Shim + 0.941*ShimdB
            + 0.01*Slope + 0.093*Tilt) * 2.8902

    The function follows the AVQI 3.01 algorithm:
    1. Concatenate sustained vowel files
    2. Concatenate continuous speech files
    3. High-pass filter (remove < 34 Hz)
    4. Extract voiced segments from continuous speech
    5. Concatenate voiced CS + last 3s of SV
    6. Compute 6 acoustic measures
    7. Calculate AVQI score
    """

    # Validate inputs
    if not sv_audio_list or not cs_audio_list:
        raise ValueError("Both sustained vowel and continuous speech audio required")

    # Load and concatenate sustained vowels
    sv_sounds = []
    for audio_data in sv_audio_list:
        sound = pm.Sound(audio_data['audio_np'], sampling_frequency=audio_data['sample_rate'])
        sv_sounds.append(sound)

    # Concatenate all SV files
    if len(sv_sounds) == 1:
        sv = sv_sounds[0]
    else:
        sv = call(sv_sounds, "Concatenate")

    sv_total_duration = call(sv, "Get total duration")

    # Load and concatenate continuous speech files
    cs_sounds = []
    for audio_data in cs_audio_list:
        sound = pm.Sound(audio_data['audio_np'], sampling_frequency=audio_data['sample_rate'])
        cs_sounds.append(sound)

    # Concatenate all CS files
    if len(cs_sounds) == 1:
        cs = cs_sounds[0]
    else:
        cs = call(cs_sounds, "Concatenate")

    cs_total_duration = call(cs, "Get total duration")

    # === PART 0: High-pass filtering ===
    cs2 = call(cs, "Filter (stop Hann band)", 0, 34, 0.1)
    sv2 = call(sv, "Filter (stop Hann band)", 0, 34, 0.1)

    # === PART 1: Extract voiced segments from continuous speech ===
    original = cs2.copy()
    sampling_rate = call(original, "Get sampling frequency")
    intermediate_samples = call(original, "Get sampling period")

    # Create empty sound for voiced segments
    only_voice = call("Create Sound", "onlyVoice", 0, 0.001, sampling_rate, "0")

    # Detect silences
    textgrid = call(original, "To TextGrid (silences)", 50, 0.003, -25, 0.1, 0.1, "silence", "sounding")

    # Extract non-silent intervals
    intervals = call([original, textgrid], "Extract intervals where", 1, False, "does not contain", "silence")

    # Concatenate non-silent parts
    only_loud = call(intervals, "Concatenate")

    signal_end = call(only_loud, "Get end time")
    window_border_left = call(only_loud, "Get start time")

    window_width = 0.03
    extreme_right = signal_end - window_width

    global_power = call(only_loud, "Get power in air")
    voiceless_threshold = global_power * (30 / 100)

    # Slide window through signal, extract voiced parts
    while window_border_left < extreme_right:
        window_border_right = window_border_left + window_width

        part = call(only_loud, "Extract part", window_border_left, window_border_right, "rectangular", 1.0, False)
        partial_power = call(part, "Get power in air")

        if partial_power > voiceless_threshold:
            # Check zero crossing rate
            zero_crossings = call(part, "Get number of zero crossings")
            duration = window_border_right - window_border_left
            zero_crossing_rate = zero_crossings / duration

            if zero_crossing_rate < 3000:
                # This is voiced - append to only_voice
                only_voice = call([only_voice, part], "Concatenate")

        window_border_left += window_width

    # === PART 2: Compute acoustic measures ===
    duration_vowel = call(sv2, "Get total duration")

    # Extract last 3 seconds of sustained vowel (or all if shorter)
    if duration_vowel > 3:
        sv3 = call(sv2, "Extract part", duration_vowel - 3, duration_vowel, "rectangular", 1.0, False)
    else:
        sv3 = sv2.copy()

    # Concatenate voiced CS + SV
    avqi_sound = call([only_voice, sv3], "Concatenate")

    duration_only_voice = call(only_voice, "Get total duration")
    duration_all = call(avqi_sound, "Get total duration")

    # Compute CPPS
    power_cepstrogram = call(avqi_sound, "To PowerCepstrogram", 60, 0.002, 5000, 50)
    cpps = call(power_cepstrogram, "Get CPPS", False, 0.01, 0.001, 60, 330, 0.05, "Parabolic", 0.001, 0, "Straight", "Robust")

    # Compute LTAS slope
    ltas = call(avqi_sound, "To Ltas", 1)
    slope = call(ltas, "Get slope", 0, 1000, 1000, 10000, "energy")

    # Compute LTAS tilt
    ltas2 = call(avqi_sound, "To Ltas", 1)
    tilt_ltas = call(ltas2, "Compute trend line", 1, 10000)
    tilt = call(tilt_ltas, "Get slope", 0, 1000, 1000, 10000, "energy")

    # Compute shimmer
    point_process = call(avqi_sound, "To PointProcess (periodic, cc)", 50, 400)
    percent_shimmer = call([avqi_sound, point_process], "Get shimmer (local)", 0, 0, 0.0001, 0.02, 1.3, 1.6)
    shim = percent_shimmer * 100  # Convert to percentage
    shdB = call([avqi_sound, point_process], "Get shimmer (local, dB)", 0, 0, 0.0001, 0.02, 1.3, 1.6)

    # Compute HNR
    pitch = call(avqi_sound, "To Pitch (cc)", 0, 75, 15, False, 0.03, 0.45, 0.01, 0.35, 0.14, 600)
    point_process2 = call([avqi_sound, pitch], "To PointProcess (cc)")
    voice_report_str = call([avqi_sound, pitch, point_process2], "Voice report", 0, 0, 75, 600, 1.3, 1.6, 0.03, 0.45)

    # Extract HNR from voice report string
    lines = voice_report_str.split('\n')
    hnr = None
    for line in lines:
        if "Mean harmonics-to-noise ratio:" in line:
            parts = line.split(':')
            if len(parts) == 2:
                value_str = parts[1].strip().split()[0]
                try:
                    hnr = float(value_str)
                except ValueError:
                    hnr = 0.0
            break

    if hnr is None:
        hnr = 0.0

    # Calculate AVQI
    avqi_score = (4.152 - (0.177 * cpps) - (0.006 * hnr) - (0.037 * shim) +
                  (0.941 * shdB) + (0.01 * slope) + (0.093 * tilt)) * 2.8902

    # Build result dictionary
    result = {
        "AVQI_VERSION": "v03.01",
        "Speaker": speaker_name if speaker_name else "NA",
        "ID": speaker_id if speaker_id else "NA",
        "DOB": speaker_dob if speaker_dob else "NA",
        "Date": assessment_date if assessment_date else "NA",
        "Sustained_vowel_duration": round(sv_total_duration, 3),
        "Continuous_speech_duration": round(cs_total_duration, 3),
        "CPPS": round(cpps, 2),
        "HNR": round(hnr, 2),
        "Shim_local": round(shim, 2),
        "Shim_local_DB": round(shdB, 2),
        "LTAS_Slope": round(slope, 2),
        "LTAS_Tilt": round(tilt, 2),
        "AVQI": round(avqi_score, 2)
    }

    return result

"""
Parselmouth-based DSI implementation (memory-based, av package integration)

This module provides a memory-based implementation of the Dysphonia Severity Index (DSI)
using Parselmouth. It processes audio directly from numpy arrays loaded by the av package,
eliminating file I/O overhead.

Based on DSI v2.01 by Youri Maryn (implementation of Wuyts et al., 2000 algorithm).
"""

import parselmouth as pm
from parselmouth.praat import call
import numpy as np


def praat_dsi_memory(
    soft_audio_list,          # List of (audio_np, sample_rate) for softest voice samples
    highpitch_audio_list,     # List of (audio_np, sample_rate) for highest pitch samples
    maxprolonged_audio_list,  # List of (audio_np, sample_rate) for maximally prolonged vowels
    stable_audio_list,        # List of (audio_np, sample_rate) for stable vowel (jitter)
    apply_calibration=False,
    calibration_db=10,
    speaker_name="",
    speaker_id="",
    speaker_dob="",
    assessment_date="",
):
    """
    Compute DSI from numpy audio arrays (memory-based, no file I/O).

    This function replicates the DSI 2.01 Praat script but operates entirely
    in memory using audio loaded by the av package.

    Parameters
    ----------
    soft_audio_list : list of dict
        List of soft voice audio segments, each dict contains:
        - 'audio_np': numpy array (float64, mono)
        - 'sample_rate': int, sampling frequency in Hz
    highpitch_audio_list : list of dict
        List of high pitch audio segments
    maxprolonged_audio_list : list of dict
        List of maximally prolonged vowel segments
    stable_audio_list : list of dict
        List of stable vowel segments (for jitter computation)
    apply_calibration : bool
        Whether to apply calibration to intensity measurements
    calibration_db : float
        Calibration factor in dB to add to intensity measurements
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
        Dictionary with DSI measurements:
        - 'ID': str, speaker ID
        - 'Maximum_phonation_time': float, maximum phonation time (s)
        - 'Softest_intensity_of_voiced_speech': float, softest intensity (dB)
        - 'Maximum_fundamental_frequency': float, highest F0 (Hz)
        - 'Jitter_ppq5': float, five-point period perturbation quotient (%)
        - 'Dysphonia_Severity_Index': float, DSI score

    Notes
    -----
    The DSI formula (Wuyts et al., 2000):
    DSI = 1.127 + 0.164*MPT - 0.038*I_low + 0.0053*F0_high - 5.30*Jitter_ppq5

    The function follows the DSI 2.01 algorithm:
    1. Determine maximum phonation time from prolonged vowels
    2. Concatenate soft voice samples, extract voiced segments, find minimum intensity
    3. Concatenate high pitch samples, find maximum F0
    4. Concatenate stable vowel samples, extract last 3s, compute jitter ppq5
    5. Calculate DSI score
    """

    # Validate inputs
    if not soft_audio_list or not highpitch_audio_list or not maxprolonged_audio_list or not stable_audio_list:
        raise ValueError("All audio lists (soft, highpitch, maxprolonged, stable) must be non-empty")

    # === PART 1: Maximum Phonation Time (MPT) ===
    # Find the longest duration among maximally prolonged vowels
    mpt = 0.0
    for audio_data in maxprolonged_audio_list:
        sound = pm.Sound(audio_data['audio_np'], sampling_frequency=audio_data['sample_rate'])
        duration = call(sound, "Get total duration")
        if duration > mpt:
            mpt = duration

    # === PART 2: Softest Intensity (I-low) ===
    # Concatenate soft voice samples
    soft_sounds = []
    for audio_data in soft_audio_list:
        sound = pm.Sound(audio_data['audio_np'], sampling_frequency=audio_data['sample_rate'])
        soft_sounds.append(sound)

    if len(soft_sounds) == 1:
        im_sound = soft_sounds[0]
    else:
        im_sound = call(soft_sounds, "Concatenate")

    # Extract voiced segments (following DSI201.praat logic)
    # Create pitch object
    im_pitch = call(im_sound, "To Pitch (cc)", 0, 70, 15, False, 0.03, 0.8, 0.01, 0.35, 0.14, 600)

    # Create point process from pitch
    im_pp = call([im_sound, im_pitch], "To PointProcess (cc)")

    # Create TextGrid marking voiced/unvoiced
    im_tg = call(im_pp, "To TextGrid (vuv)", 0.02, 0.01)

    # Extract voiced intervals
    im_intervals = call([im_sound, im_tg], "Extract intervals where", 1, False, "is equal to", "V")

    # Concatenate voiced parts
    im_voiced = call(im_intervals, "Concatenate")

    # Compute intensity
    if apply_calibration:
        im_intensity = call(im_voiced, "To Intensity", 60, 0.0, True)
        # Apply calibration: Formula... 1*self+calibration
        call(im_intensity, "Formula", f"1*self+{calibration_db}")
    else:
        im_intensity = call(im_voiced, "To Intensity", 60, 0.0, True)

    minimum_intensity = call(im_intensity, "Get minimum", 0, 0, "None")

    # === PART 3: Highest F0 (F0-high) ===
    # Concatenate high pitch samples
    fh_sounds = []
    for audio_data in highpitch_audio_list:
        sound = pm.Sound(audio_data['audio_np'], sampling_frequency=audio_data['sample_rate'])
        fh_sounds.append(sound)

    if len(fh_sounds) == 1:
        fh_sound = fh_sounds[0]
    else:
        fh_sound = call(fh_sounds, "Concatenate")

    # Detect pitch with extended range for high F0
    fh_pitch = call(fh_sound, "To Pitch (cc)", 0, 70, 15, False, 0.03, 0.8, 0.01, 0.35, 0.14, 1300)
    maximum_f0 = call(fh_pitch, "Get maximum", 0, 0, "Hertz", "None")

    # === PART 4: Jitter ppq5 ===
    # Concatenate stable vowel samples
    ppq_sounds = []
    for audio_data in stable_audio_list:
        sound = pm.Sound(audio_data['audio_np'], sampling_frequency=audio_data['sample_rate'])
        ppq_sounds.append(sound)

    if len(ppq_sounds) == 1:
        ppq_sound = ppq_sounds[0]
    else:
        ppq_sound = call(ppq_sounds, "Concatenate")

    # Extract last 3 seconds (or use all if shorter)
    duration_vowel = call(ppq_sound, "Get total duration")
    duration_start = duration_vowel - 3.0

    if duration_vowel > 3.0:
        ppq2_sound = call(ppq_sound, "Extract part", duration_start, duration_vowel, "rectangular", 1.0, False)
    else:
        ppq2_sound = ppq_sound.copy()

    # Create pitch object
    ppq_pitch = call(ppq2_sound, "To Pitch", 0, 70, 600)

    # Create point process
    ppq_pp = call([ppq2_sound, ppq_pitch], "To PointProcess (cc)")

    # Get voice report to extract jitter ppq5
    voice_report_str = call(
        [ppq2_sound, ppq_pitch, ppq_pp],
        "Voice report",
        0, 0,  # time range (0 = all)
        70, 600,  # F0 range
        1.3, 1.6,  # period/amplitude factors
        0.03, 0.45  # silence/voicing thresholds
    )

    # Extract jitter ppq5 from voice report string
    lines = voice_report_str.split('\n')
    jitter_ppq5 = None
    for line in lines:
        if "Jitter (ppq5):" in line:
            parts = line.split(':')
            if len(parts) == 2:
                # Extract percentage value (e.g., "0.534%")
                value_str = parts[1].strip().split()[0]
                try:
                    jitter_ppq5_fraction = float(value_str)
                    # Convert to percentage (Praat returns as fraction, need to multiply by 100)
                    jitter_ppq5 = jitter_ppq5_fraction * 100
                except ValueError:
                    jitter_ppq5 = 0.0
            break

    if jitter_ppq5 is None:
        jitter_ppq5 = 0.0

    # === PART 5: Calculate DSI ===
    # Formula from DSI201.praat line 220:
    # dsi2 = 1.127+ 0.164*mpt - 0.038*minimumIntensity + 0.0053*maximumF0 - 5.30*jitterPpq5
    dsi_score = (1.127 +
                 0.164 * mpt -
                 0.038 * minimum_intensity +
                 0.0053 * maximum_f0 -
                 5.30 * jitter_ppq5)

    # Build result dictionary
    result = {
        "ID": speaker_id if speaker_id else "NA",
        "Maximum_phonation_time": round(mpt, 2),
        "Softest_intensity_of_voiced_speech": round(minimum_intensity, 2),
        "Maximum_fundamental_frequency": round(maximum_f0, 2),
        "Jitter_ppq5": round(jitter_ppq5, 2),
        "Dysphonia_Severity_Index": round(dsi_score, 2)
    }

    return result

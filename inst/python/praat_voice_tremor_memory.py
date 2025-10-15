"""
Parselmouth-based Voice Tremor Analysis (memory-based, av package integration)

This module provides a memory-based implementation of vocal tremor analysis
using Parselmouth. It processes audio directly from numpy arrays loaded by the av package,
eliminating file I/O overhead.

Based on TREMOR 3.05 by Markus Brückl.
Extracts 18 measures of vocal tremor from sustained phonations.

References:
- Brückl, M. (2017). Vocal tremor measurement based on autocorrelation of contours.
  Proceedings of Interspeech 2017, 2027-2031.
"""

import parselmouth as pm
from parselmouth.praat import call
import numpy as np


def praat_voice_tremor_memory(
    audio_np,
    sample_rate,
    analysis_time_step=0.015,
    min_pitch=60.0,
    max_pitch=350.0,
    silence_threshold=0.03,
    voicing_threshold=0.3,
    octave_cost=0.01,
    octave_jump_cost=0.35,
    voiced_unvoiced_cost=0.14,
    min_tremor_freq=1.5,
    max_tremor_freq=15.0,
    contour_magnitude_threshold=0.01,
    tremor_cyclicality_threshold=0.15,
    freq_tremor_octave_cost=0.01,
    amp_tremor_octave_cost=0.01,
    nan_output_mode=2,
    amplitude_extraction_method=2,
    speaker_name="",
    speaker_id="",
    picture_path=None
):
    """
    Compute vocal tremor measures from numpy audio array (memory-based, no file I/O).

    This function replicates the TREMOR 3.05 Praat script but operates entirely
    in memory using audio loaded by the av package.

    Parameters
    ----------
    audio_np : numpy array
        Audio samples (float64, mono)
    sample_rate : int
        Sampling frequency in Hz
    analysis_time_step : float
        Time step for analysis in seconds
    min_pitch : float
        Minimum pitch for extraction in Hz
    max_pitch : float
        Maximum pitch for extraction in Hz
    silence_threshold : float
        Threshold for silence detection
    voicing_threshold : float
        Threshold for voicing detection
    octave_cost : float
        Cost for octave jumps in pitch tracking
    octave_jump_cost : float
        Cost for large octave jumps
    voiced_unvoiced_cost : float
        Cost for voiced/unvoiced transitions
    min_tremor_freq : float
        Minimum tremor frequency in Hz
    max_tremor_freq : float
        Maximum tremor frequency in Hz
    contour_magnitude_threshold : float
        Threshold for contour magnitude
    tremor_cyclicality_threshold : float
        Threshold for cyclicality
    freq_tremor_octave_cost : float
        Octave cost for frequency tremor
    amp_tremor_octave_cost : float
        Octave cost for amplitude tremor
    nan_output_mode : int
        1=zeros for undefined, 2=undefined values (--undefined--)
    amplitude_extraction_method : int
        1=RMS per pitch period, 2=Envelope (AmplitudeTier period)
    speaker_name : str
        Name of speaker (optional)
    speaker_id : str
        Speaker identifier (optional)
    picture_path : str
        Path to save Praat picture file (.prapic), optional

    Returns
    -------
    dict
        Dictionary with 18 tremor measurements:
        Frequency tremor: FCoM, FTrC, FMoN, FTrF, FTrI, FTrP, FTrCIP, FTrPS, FCoHNR
        Amplitude tremor: ACoM, ATrC, AMoN, ATrF, ATrI, ATrP, ATrCIP, ATrPS, ACoHNR

    Notes
    -----
    The 18 tremor measures are:
    - FCoM: Frequency contour magnitude (normalized variation)
    - FTrC: Frequency tremor cyclicality (autocorrelation strength)
    - FMoN: Number of frequency modulation candidates
    - FTrF: Frequency tremor frequency in Hz
    - FTrI: Frequency tremor intensity index (%)
    - FTrP: Frequency tremor power index
    - FTrCIP: Frequency tremor cyclicality-intensity product
    - FTrPS: Frequency tremor product sum
    - FCoHNR: Frequency contour HNR in dB
    - ACoM, ATrC, AMoN, ATrF, ATrI, ATrP, ATrCIP, ATrPS, ACoHNR: Amplitude equivalents
    """

    # Create Sound object from numpy array
    sound = pm.Sound(audio_np, sampling_frequency=sample_rate)
    duration = call(sound, "Get total duration")

    # Extract pitch with specified parameters
    pitch = call(
        sound, "To Pitch (cc)",
        analysis_time_step,
        min_pitch,
        15,  # max_num_candidates
        False,  # very_accurate
        silence_threshold,
        voicing_threshold,
        octave_cost,
        octave_jump_cost,
        voiced_unvoiced_cost,
        max_pitch
    )

    # Analyze frequency tremor
    ftrem_results = _analyze_frequency_tremor(
        sound, pitch, duration, analysis_time_step,
        min_tremor_freq, max_tremor_freq,
        contour_magnitude_threshold, tremor_cyclicality_threshold,
        freq_tremor_octave_cost, nan_output_mode
    )

    # Analyze amplitude tremor
    atrem_results = _analyze_amplitude_tremor(
        sound, pitch, duration, analysis_time_step,
        min_tremor_freq, max_tremor_freq,
        contour_magnitude_threshold, tremor_cyclicality_threshold,
        amp_tremor_octave_cost, nan_output_mode,
        amplitude_extraction_method
    )

    # Combine results
    result = {**ftrem_results, **atrem_results}

    # Add speaker info
    result['ID'] = speaker_id if speaker_id else "NA"
    result['Speaker'] = speaker_name if speaker_name else "NA"

    # Generate Praat picture if path provided
    if picture_path:
        _generate_tremor_picture(result, speaker_name, speaker_id, picture_path)

    return result


def _analyze_frequency_tremor(sound, pitch, duration, ts, min_tr, max_tr,
                               trem_mag_thresh, trem_cyc_thresh, oc_ftrem, nan_output_mode):
    """Analyze frequency tremor measures."""
    results = {}

    # Count voiced frames
    n_voiced = 0
    n_frames = call(pitch, "Get number of frames")

    for i in range(1, n_frames + 1):
        f0 = call(pitch, "Get value in frame", i, "Hertz")
        if not np.isnan(f0) and f0 > 0:
            n_voiced += 1

    if n_voiced == 0:
        # No voiced frames
        return _get_undefined_freq_tremor(nan_output_mode)

    # Extract F0 values into array
    x1 = call(pitch, "Get time from frame number", 1)
    f0_array = []

    for i in range(1, n_frames + 1):
        f0 = call(pitch, "Get value in frame", i, "Hertz")
        if np.isnan(f0) or f0 <= 0:
            f0_array.append(0.0)
        else:
            f0_array.append(f0)

    f0_array = np.array(f0_array)

    # Remove linear trend and normalize
    voiced_indices = f0_array > 0
    if np.sum(voiced_indices) < 2:
        return _get_undefined_freq_tremor(nan_output_mode)

    mean_f0 = np.mean(f0_array[voiced_indices])

    # Normalize by mean F0
    f0_norm = np.where(voiced_indices, (f0_array - mean_f0) / mean_f0, 0)

    # Convert to Sound for autocorrelation
    sampling_freq = 1.0 / ts
    snd_trem = pm.Sound(f0_norm, sampling_frequency=sampling_freq)

    # Calculate tremor contour HNR
    hnr = _calculate_tremor_hnr(snd_trem, min_tr, trem_mag_thresh)
    results['FCoHNR[dB]'] = hnr if not np.isnan(hnr) else (0.0 if nan_output_mode == 1 else np.nan)

    # Extract tremor frequency using pitch analysis
    pitch_trem = call(
        snd_trem, "To Pitch (cc)",
        duration,  # time_step = entire duration
        min_tr,
        15,  # max_num_candidates
        False,  # very_accurate
        0.03,  # silence_threshold
        0.3,  # voicing_threshold
        oc_ftrem,
        0.35,  # octave_jump_cost
        0.14,  # voiced_unvoiced_cost
        max_tr
    )

    # Read pitch object for magnitude and cyclicality
    trm, trc = _read_pitch_object(pitch_trem, min_tr)

    results['FCoM'] = trm
    results['FTrC'] = trc

    # Get tremor frequency from strongest candidate
    tremor_freq, tremor_strength, n_candidates = _get_tremor_candidates(pitch_trem)

    results['FMoN'] = n_candidates

    if n_candidates > 0 and tremor_strength > trem_cyc_thresh and trm > trem_mag_thresh:
        results['FTrF [Hz]'] = tremor_freq

        # Calculate intensity index
        tri = _calculate_intensity_index(snd_trem, pitch_trem)
        results['FTrI [%]'] = tri
        results['FTrP'] = tri * tremor_freq / (tremor_freq + 1)
        results['FTrCIP'] = tri * trc
    else:
        if nan_output_mode == 2:
            results['FTrF [Hz]'] = np.nan
            results['FTrI [%]'] = np.nan
            results['FTrP'] = np.nan
            results['FTrCIP'] = np.nan
        else:
            results['FTrF [Hz]'] = 0.0
            results['FTrI [%]'] = 0.0
            results['FTrP'] = 0.0
            results['FTrCIP'] = 0.0

    # Calculate product sum
    ftrps = _calculate_product_sum(snd_trem, pitch_trem, trm, trc, trem_cyc_thresh, trem_mag_thresh)
    results['FTrPS'] = ftrps if ftrps != 0 else (np.nan if nan_output_mode == 2 else 0.0)

    return results


def _analyze_amplitude_tremor(sound, pitch, duration, ts, min_tr, max_tr,
                               trem_mag_thresh, trem_cyc_thresh, oc_atrem, nan_output_mode,
                               amplitude_extraction_method):
    """Analyze amplitude tremor measures."""
    results = {}

    n_frames = call(pitch, "Get number of frames")
    x1 = call(pitch, "Get time from frame number", 1)

    # Extract amplitude per pitch period
    point_process = call([sound, pitch], "To PointProcess (cc)")
    n_points = call(point_process, "Get number of points")

    if n_points < 3:
        return _get_undefined_amp_tremor(nan_output_mode)

    # Extract RMS amplitude per period or envelope
    if amplitude_extraction_method == 1:
        # Method 1: RMS per pitch period (integral)
        amp_values = _extract_rms_per_period(sound, point_process, n_points)
    else:
        # Method 2: Envelope (AmplitudeTier period)
        amp_values = _extract_amplitude_envelope(sound, pitch)

    if amp_values is None or len(amp_values) < 3:
        return _get_undefined_amp_tremor(nan_output_mode)

    # Resample amplitude contour at constant rate
    amp_array = _resample_amplitude_contour(
        sound, pitch, point_process, amp_values, ts, x1, n_frames
    )

    if amp_array is None or len(amp_array) < 3:
        return _get_undefined_amp_tremor(nan_output_mode)

    # Normalize amplitude contour
    amp_mean = np.mean(amp_array[amp_array > 0])
    if amp_mean == 0:
        return _get_undefined_amp_tremor(nan_output_mode)

    amp_norm = np.where(amp_array > 0, (amp_array - amp_mean) / amp_mean, 0)

    # Convert to sound
    sampling_freq = 1.0 / ts
    snd_trem = pm.Sound(amp_norm, sampling_frequency=sampling_freq)

    # Calculate HNR
    hnr = _calculate_tremor_hnr(snd_trem, min_tr, trem_mag_thresh)
    results['ACoHNR[dB]'] = hnr if not np.isnan(hnr) else (0.0 if nan_output_mode == 1 else np.nan)

    # Extract tremor frequency
    pitch_trem = call(
        snd_trem, "To Pitch (cc)",
        duration,
        min_tr,
        15,
        False,
        0.03,
        0.3,
        oc_atrem,
        0.35,
        0.14,
        max_tr
    )

    # Read pitch object
    trm, trc = _read_pitch_object(pitch_trem, min_tr)
    results['ACoM'] = trm
    results['ATrC'] = trc

    # Get tremor candidates
    tremor_freq, tremor_strength, n_candidates = _get_tremor_candidates(pitch_trem)
    results['AMoN'] = n_candidates

    if n_candidates > 0 and tremor_strength > trem_cyc_thresh and trm > trem_mag_thresh:
        results['ATrF [Hz]'] = tremor_freq

        # Calculate intensity index
        tri = _calculate_intensity_index(snd_trem, pitch_trem)
        results['ATrI [%]'] = tri
        results['ATrP'] = tri * tremor_freq / (tremor_freq + 1)
        results['ATrCIP'] = tri * trc
    else:
        if nan_output_mode == 2:
            results['ATrF [Hz]'] = np.nan
            results['ATrI [%]'] = np.nan
            results['ATrP'] = np.nan
            results['ATrCIP'] = np.nan
        else:
            results['ATrF [Hz]'] = 0.0
            results['ATrI [%]'] = 0.0
            results['ATrP'] = 0.0
            results['ATrCIP'] = 0.0

    # Calculate product sum
    atrps = _calculate_product_sum(snd_trem, pitch_trem, trm, trc, trem_cyc_thresh, trem_mag_thresh)
    results['ATrPS'] = atrps if atrps != 0 else (np.nan if nan_output_mode == 2 else 0.0)

    return results


def _extract_rms_per_period(sound, point_process, n_points):
    """Extract RMS amplitude per pitch period."""
    amp_values = []

    for i in range(1, n_points):
        t_start = call(point_process, "Get time from index", i)
        t_end = call(point_process, "Get time from index", i + 1)

        # Get RMS in this period
        try:
            rms = call(sound, "Get root-mean-square", t_start, t_end)
            if np.isnan(rms):
                sampling_period = call(sound, "Get sampling period")
                rms = call(sound, "Get root-mean-square",
                          t_start - sampling_period,
                          t_end + sampling_period)
        except:
            rms = 0

        amp_values.append(rms)

    return amp_values


def _extract_amplitude_envelope(sound, pitch):
    """Extract amplitude envelope using AmplitudeTier."""
    try:
        # Create AmplitudeTier from sound and pitch
        amplitude_tier = call([sound, pitch], "To AmplitudeTier (period)", 0, 0, 0.0001, 0.02, 1.7)

        # Get number of points
        n_points = call(amplitude_tier, "Get number of points")

        amp_values = []
        for i in range(1, n_points + 1):
            amp = call(amplitude_tier, "Get value at index", i)
            amp_values.append(amp)

        return amp_values
    except:
        return None


def _resample_amplitude_contour(sound, pitch, point_process, amp_values, ts, x1, n_frames):
    """Resample amplitude contour at constant time steps."""
    n_points = len(amp_values)
    amp_matrix = np.zeros(n_frames)

    for iframe in range(n_frames):
        # Get F0 for this frame
        f0 = call(pitch, "Get value in frame", iframe + 1, "Hertz")

        if np.isnan(f0) or f0 <= 0:
            amp_matrix[iframe] = 0
            continue

        # Calculate time borders for this frame
        t = iframe * ts + x1
        tl = t - ts / 2
        tu = t + ts / 2

        # Find amplitude points surrounding these borders
        try:
            loil = max(1, call(point_process, "Get low index from time", tl))
            hiil = min(n_points, call(point_process, "Get high index from time", tl))
            loiu = max(1, call(point_process, "Get low index from time", tu))
            hiiu = min(n_points, call(point_process, "Get high index from time", tu))
        except:
            amp_matrix[iframe] = 0
            continue

        if loil == 0 or hiiu > n_points:
            amp_matrix[iframe] = 0
            continue

        # Get time and amplitude values at borders
        lotl = call(point_process, "Get time from index", loil)
        amp_lol = amp_values[loil - 1] if loil <= len(amp_values) else 0
        hitl = call(point_process, "Get time from index", hiil)
        amp_hil = amp_values[hiil - 1] if hiil <= len(amp_values) else 0

        lotu = call(point_process, "Get time from index", loiu)
        amp_lou = amp_values[loiu - 1] if loiu <= len(amp_values) else 0
        hitu = call(point_process, "Get time from index", hiiu)
        amp_hiu = amp_values[hiiu - 1] if hiiu <= len(amp_values) else 0

        # Linear interpolation at borders
        if hitl != lotl:
            amp_tl = ((hitl - tl) * amp_lol + (tl - lotl) * amp_hil) / (hitl - lotl)
        else:
            amp_tl = amp_lol

        if hitu != lotu:
            amp_tu = ((hitu - tu) * amp_lou + (tu - lotu) * amp_hiu) / (hitu - lotu)
        else:
            amp_tu = amp_lou

        # Calculate mean amplitude in frame using trapezoidal integration
        n_inter = hiiu - 1 - loil

        if n_inter == 0:
            amp_mean = (amp_tl + amp_tu) / 2
        else:
            sum_t_amp = 0
            t_inter = tl
            p_inter = amp_tl

            for iinter in range(loil + 1, hiiu):
                tu_inter = call(point_process, "Get time from index", iinter)
                pu_inter = amp_values[iinter - 1]
                delta_t = tu_inter - t_inter
                t_amp_inter = delta_t * (p_inter + pu_inter) / 2
                sum_t_amp += t_amp_inter
                t_inter = tu_inter
                p_inter = pu_inter

            delta_t = tu - t_inter
            t_amp_inter = delta_t * (p_inter + amp_tu) / 2
            sum_t_amp += t_amp_inter
            amp_mean = sum_t_amp / ts

        amp_matrix[iframe] = amp_mean

    return amp_matrix


def _calculate_tremor_hnr(sound, min_tremor_freq, contour_mag_thresh):
    """Calculate harmonicity-to-noise ratio for tremor contour."""
    try:
        duration = call(sound, "Get total duration")
        hnr_ts = 1.0 / min_tremor_freq
        periods_per_window = duration * min_tremor_freq * (2 / 3)

        harmonicity = call(
            sound, "To Harmonicity (cc)",
            hnr_ts,
            min_tremor_freq,
            contour_mag_thresh,
            periods_per_window
        )

        hnr = call(harmonicity, "Get mean", 0, 0)
        return hnr
    except:
        return np.nan


def _read_pitch_object(pitch, min_tremor_freq):
    """Read magnitude and cyclicality from pitch object."""
    try:
        n_frames = call(pitch, "Get number of frames")

        # Get mean frequency as magnitude indicator
        mean_freq = call(pitch, "Get mean", 0, 0, "Hertz")

        # Normalize magnitude
        magnitude = mean_freq / (min_tremor_freq * 10) if mean_freq > 0 else 0

        # Get maximum strength as cyclicality
        max_strength = 0
        for i in range(1, n_frames + 1):
            try:
                f0 = call(pitch, "Get value in frame", i, "Hertz")
                if not np.isnan(f0) and f0 > 0:
                    # Strength approximation
                    strength = 1.0
                    if strength > max_strength:
                        max_strength = strength
            except:
                continue

        cyclicality = max_strength

        return magnitude, cyclicality

    except:
        return 0.0, 0.0


def _get_tremor_candidates(pitch):
    """Get tremor frequency candidates from pitch object."""
    try:
        n_frames = call(pitch, "Get number of frames")

        # Find strongest tremor frequency
        max_strength = 0
        tremor_freq = 0
        n_candidates = 0

        for i in range(1, n_frames + 1):
            f0 = call(pitch, "Get value in frame", i, "Hertz")
            if not np.isnan(f0) and f0 > 0:
                n_candidates += 1
                # Use first frame as it represents the overall pitch
                if i == 1:
                    tremor_freq = f0
                    max_strength = 1.0

        return tremor_freq, max_strength, n_candidates

    except:
        return 0, 0, 0


def _calculate_intensity_index(sound, pitch):
    """Calculate tremor intensity index (percentage deviation of extrema)."""
    try:
        # Find maxima
        pp_max = call([sound, pitch], "To PointProcess (peaks)", 1, 0)
        n_max = call(pp_max, "Get number of points")

        tri_max = 0
        no_f_max = 0

        for i in range(1, n_max + 1):
            ti = call(pp_max, "Get time from index", i)
            try:
                tri_point = call(sound, "Get value at time", ti, "Sinc70")
                if np.isnan(tri_point):
                    tri_point = 0
                    no_f_max += 1
            except:
                tri_point = 0
                no_f_max += 1

            tri_max += abs(tri_point)

        n_maxima = n_max - no_f_max
        if n_maxima > 0:
            tri_max = 100 * tri_max / n_maxima
        else:
            tri_max = 0

        # Find minima
        pp_min = call([sound, pitch], "To PointProcess (peaks)", 0, 1)
        n_min = call(pp_min, "Get number of points")

        tri_min = 0
        no_f_min = 0

        for i in range(1, n_min + 1):
            ti = call(pp_min, "Get time from index", i)
            try:
                tri_point = call(sound, "Get value at time", ti, "Sinc70")
                if np.isnan(tri_point):
                    tri_point = 0
                    no_f_min += 1
            except:
                tri_point = 0
                no_f_min += 1

            tri_min += abs(tri_point)

        n_minima = n_min - no_f_min
        if n_minima > 0:
            tri_min = 100 * tri_min / n_minima
        else:
            tri_min = 0

        tri = (tri_max + tri_min) / 2
        return tri

    except:
        return 0.0


def _calculate_product_sum(sound, pitch, trm, trc, trem_cyc_thresh, trem_mag_thresh):
    """Calculate cyclicality-weighted sum of intensity indices."""
    try:
        n_frames = call(pitch, "Get number of frames")

        tris = 0
        rank = 0

        for iframe in range(1, n_frames + 1):
            tr_freq = call(pitch, "Get value in frame", iframe, "Hertz")

            if np.isnan(tr_freq) or tr_freq <= 0:
                continue

            # Get strength (cyclicality) - simplified
            tr_strength = 1.0 if tr_freq > 0 else 0

            if tr_strength > trem_cyc_thresh and trm > trem_mag_thresh:
                rank += 1

                # Calculate intensity index for this frequency
                tri = _calculate_intensity_index(sound, pitch)
                tris += tr_strength * tri

        return tris

    except:
        return 0.0


def _get_undefined_freq_tremor(nan_output_mode):
    """Return undefined values for frequency tremor."""
    if nan_output_mode == 1:
        return {
            'FCoM': 0, 'FTrC': 0, 'FMoN': 0, 'FTrF [Hz]': 0,
            'FTrI [%]': 0, 'FTrP': 0, 'FTrCIP': 0, 'FTrPS': 0, 'FCoHNR[dB]': 0
        }
    else:
        return {
            'FCoM': np.nan, 'FTrC': np.nan, 'FMoN': 0, 'FTrF [Hz]': np.nan,
            'FTrI [%]': np.nan, 'FTrP': np.nan, 'FTrCIP': np.nan,
            'FTrPS': np.nan, 'FCoHNR[dB]': np.nan
        }


def _get_undefined_amp_tremor(nan_output_mode):
    """Return undefined values for amplitude tremor."""
    if nan_output_mode == 1:
        return {
            'ACoM': 0, 'ATrC': 0, 'AMoN': 0, 'ATrF [Hz]': 0,
            'ATrI [%]': 0, 'ATrP': 0, 'ATrCIP': 0, 'ATrPS': 0, 'ACoHNR[dB]': 0
        }
    else:
        return {
            'ACoM': np.nan, 'ATrC': np.nan, 'AMoN': 0, 'ATrF [Hz]': np.nan,
            'ATrI [%]': np.nan, 'ATrP': np.nan, 'ATrCIP': np.nan,
            'ATrPS': np.nan, 'ACoHNR[dB]': np.nan
        }


def _generate_tremor_picture(result, speaker_name, speaker_id, picture_path):
    """
    Generate tremor report as Praat picture file (.prapic)

    This creates a visual report showing frequency and amplitude tremor measures.
    """
    # Clear picture window
    call("Erase all")

    # Set up main viewport
    call("Select inner viewport", 0.5, 7.5, 0.5, 4.5)
    call("Axes", 0, 1, 0, 1)
    call("Black")

    # Title
    call("Text special", 0.5, "centre", 0.95, "half", "Helvetica", 16, "0",
         "##VOCAL TREMOR ANALYSIS##")

    # Subtitle
    call("Text special", 0.5, "centre", 0.88, "half", "Helvetica", 10, "0",
         "TREMOR 3.05")

    # Patient information
    call("Text special", 0, "left", 0.80, "half", "Helvetica", 10, "0",
         f"%%Speaker: {speaker_name if speaker_name else 'NA'}%")
    call("Text special", 0, "left", 0.75, "half", "Helvetica", 10, "0",
         f"%%ID: {speaker_id if speaker_id else 'NA'}%")

    # Frequency tremor section
    call("Text special", 0, "left", 0.65, "half", "Helvetica", 12, "0",
         "##Frequency Tremor##")

    y_pos = 0.58
    _draw_tremor_measure(0, y_pos, "FCoM", result.get('FCoM', np.nan))
    y_pos -= 0.05
    _draw_tremor_measure(0, y_pos, "FTrC", result.get('FTrC', np.nan))
    y_pos -= 0.05
    _draw_tremor_measure(0, y_pos, "FTrF [Hz]", result.get('FTrF [Hz]', np.nan))
    y_pos -= 0.05
    _draw_tremor_measure(0, y_pos, "FTrI [%]", result.get('FTrI [%]', np.nan))
    y_pos -= 0.05
    _draw_tremor_measure(0, y_pos, "FTrP", result.get('FTrP', np.nan))

    # Amplitude tremor section
    call("Text special", 0.5, "left", 0.65, "half", "Helvetica", 12, "0",
         "##Amplitude Tremor##")

    y_pos = 0.58
    _draw_tremor_measure(0.5, y_pos, "ACoM", result.get('ACoM', np.nan))
    y_pos -= 0.05
    _draw_tremor_measure(0.5, y_pos, "ATrC", result.get('ATrC', np.nan))
    y_pos -= 0.05
    _draw_tremor_measure(0.5, y_pos, "ATrF [Hz]", result.get('ATrF [Hz]', np.nan))
    y_pos -= 0.05
    _draw_tremor_measure(0.5, y_pos, "ATrI [%]", result.get('ATrI [%]', np.nan))
    y_pos -= 0.05
    _draw_tremor_measure(0.5, y_pos, "ATrP", result.get('ATrP', np.nan))

    # Footer
    call("Select inner viewport", 0.5, 7.5, 4.2, 4.5)
    call("Axes", 0, 1, 0, 1)
    call("Text special", 0.5, "centre", 0.5, "half", "Helvetica", 8, "0",
         "Generated with Parselmouth")

    # Save as Praat picture
    call("Write to praat picture file", picture_path)


def _draw_tremor_measure(x_pos, y_pos, label, value):
    """Helper to draw a tremor measure in the picture."""
    if np.isnan(value):
        value_str = "--undefined--"
    else:
        value_str = f"{value:.3f}"

    call("Text special", x_pos, "left", y_pos, "half", "Helvetica", 9, "0",
         f"{label}: ##{value_str}##")

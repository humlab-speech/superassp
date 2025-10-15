"""
Parselmouth-based PraatSauce implementation (memory-based, av package integration)

This module provides a memory-based implementation of PraatSauce voice quality analysis
using Parselmouth. It processes audio directly from numpy arrays loaded by the av package,
eliminating file I/O overhead.

PraatSauce extracts spectral voice quality measures including:
- Pitch (f0)
- Formants (F1-F3) and bandwidths (B1-B3)
- Spectral measures (H1, H2, H4, A1-A3, CPP, HNR, etc.)

Based on PraatSauce by James Kirby and Fredrik Karlsson
"""

import parselmouth as pm
from parselmouth.praat import call
import numpy as np


def praat_sauce_memory(
    audio_np,
    sample_rate,
    window_shift=5.0,  # milliseconds between measurements
    window_size=25.0,  # window length in milliseconds
    f0_min=50.0,
    f0_max=300.0,
    max_formant_hz=5000.0,
    num_formants=5,
    pre_emph_from=50.0,
    formant_tracking=True,
    f1_ref=500.0,
    f2_ref=1500.0,
    f3_ref=2500.0,
    use_bandwidth_formula=False,
    resample_to_16k=True,
    time_step=0.0,
    spectral_measures=True,
    formant_measures=True,
    pitch_tracking=True
):
    """
    Extract voice quality measures from audio using Parselmouth (memory-based).

    This function replicates the PraatSauce workflow but operates entirely
    in memory using audio loaded by the av package.

    Parameters
    ----------
    audio_np : numpy array
        Audio samples (float64, mono)
    sample_rate : int
        Sampling frequency in Hz
    window_shift : float
        Time step between measurements in milliseconds (default 5ms)
    window_size : float
        Analysis window length in milliseconds (default 25ms)
    f0_min : float
        Minimum F0 in Hz (default 50)
    f0_max : float
        Maximum F0 in Hz (default 300)
    max_formant_hz : float
        Maximum formant frequency in Hz (default 5000)
    num_formants : int
        Number of formants to track (default 5)
    pre_emph_from : float
        Pre-emphasis frequency in Hz (default 50)
    formant_tracking : bool
        Use formant tracking for cleaner tracks (default True)
    f1_ref : float
        Reference F1 frequency for tracking (default 500)
    f2_ref : float
        Reference F2 frequency for tracking (default 1500)
    f3_ref : float
        Reference F3 frequency for tracking (default 2500)
    use_bandwidth_formula : bool
        Use Hawks & Miller bandwidth formula (default False)
    resample_to_16k : bool
        Resample audio to 16kHz (default True)
    time_step : float
        Time step for formant analysis (0 = auto)
    spectral_measures : bool
        Compute spectral measures (default True)
    formant_measures : bool
        Compute formant measures (default True)
    pitch_tracking : bool
        Compute pitch (default True)

    Returns
    -------
    dict
        Dictionary with measurement arrays:
        - 't': time points in seconds
        - 'f0': fundamental frequency in Hz (if pitch_tracking=True)
        - 'F1', 'F2', 'F3': formant frequencies (if formant_measures=True)
        - 'B1', 'B2', 'B3': formant bandwidths (if formant_measures=True)
        - Spectral measures (if spectral_measures=True):
          H1u, H2u, H4u, H2Ku, H5Ku, A1u, A2u, A3u,
          H1H2u, H2H4u, H1A1u, H1A2u, H1A3u, H2KH5Ku,
          H1c, H2c, H4c, A1c, A2c, A3c,
          H1H2c, H2H4c, H1A1c, H1A2c, H1A3c,
          CPP, HNR05, HNR15, HNR25, HNR35

    Notes
    -----
    Spectral measures require both pitch and formant analysis.
    """

    # Create Sound object from numpy array
    sound = pm.Sound(audio_np, sampling_frequency=sample_rate)

    # Optionally resample to 16kHz (speeds up processing)
    if resample_to_16k and sample_rate != 16000:
        sound = call(sound, "Resample", 16000, 50)

    duration = call(sound, "Get total duration")

    # Calculate number of timepoints
    window_shift_sec = window_shift / 1000.0
    timepoints = int(duration / window_shift_sec)

    if timepoints < 1:
        timepoints = 1

    # Initialize results dictionary
    results = {
        't': []
    }

    # Create Pitch object if needed
    pitch = None
    if pitch_tracking or spectral_measures:
        pitch = call(
            sound, "To Pitch (ac)",
            time_step if time_step > 0 else 0,
            f0_min,
            15,  # max_num_candidates
            False,  # very_accurate
            0.03,  # silence_threshold
            0.45,  # voicing_threshold
            0.01,  # octave_cost
            0.35,  # octave_jump_cost
            0.14,  # voiced_unvoiced_cost
            f0_max
        )

        if pitch_tracking:
            results['f0'] = []

    # Create Formant object if needed
    formant = None
    if formant_measures or spectral_measures:
        window_length = window_size / 1000.0  # Convert to seconds
        formant = call(
            sound, "To Formant (burg)",
            time_step if time_step > 0 else 0,
            num_formants,
            max_formant_hz,
            window_length,
            pre_emph_from
        )

        # Apply formant tracking if requested
        if formant_tracking:
            min_formants = call(formant, "Get minimum number of formants")
            if min_formants >= 3:
                formant = call(formant, "Track", 3, f1_ref, f2_ref, f3_ref, 3850, 4950, 1, 1, 1)
            elif min_formants == 2:
                formant = call(formant, "Track", 2, f1_ref, f2_ref, f3_ref, 3850, 4950, 1, 1, 1)

        if formant_measures:
            results['F1'] = []
            results['F2'] = []
            results['F3'] = []
            results['B1'] = []
            results['B2'] = []
            results['B3'] = []

    # Create Harmonicity objects for HNR if needed
    hnr_objects = {}
    if spectral_measures:
        # Initialize spectral measure arrays
        spectral_keys = [
            'H1u', 'H2u', 'H4u', 'H2Ku', 'H5Ku', 'A1u', 'A2u', 'A3u',
            'H1H2u', 'H2H4u', 'H1A1u', 'H1A2u', 'H1A3u', 'H2KH5Ku',
            'H1c', 'H2c', 'H4c', 'A1c', 'A2c', 'A3c',
            'H1H2c', 'H2H4c', 'H1A1c', 'H1A2c', 'H1A3c',
            'CPP', 'HNR05', 'HNR15', 'HNR25', 'HNR35'
        ]
        for key in spectral_keys:
            results[key] = []

        # Create HNR objects for different frequency bands
        for band, max_freq in [(500, '05'), (1500, '15'), (2500, '25'), (3500, '35')]:
            filtered = call(sound, "Filter (pass Hann band)", 0, band, 100)
            hnr = call(filtered, "To Harmonicity (cc)", 0.01, f0_min, 0.1, 1.0)
            hnr_objects[f'HNR{max_freq}'] = hnr

    # Extract measurements at each timepoint
    for i in range(timepoints):
        t = i * window_shift_sec
        results['t'].append(t)

        # Extract F0
        if pitch_tracking and pitch is not None:
            f0 = call(pitch, "Get value at time", t, "Hertz", "Linear")
            if f0 is None or np.isnan(f0):
                f0 = np.nan
            results['f0'].append(f0)

        # Extract formants and bandwidths
        if formant_measures and formant is not None:
            for formant_num in [1, 2, 3]:
                freq = call(formant, "Get value at time", formant_num, t, "Hertz", "Linear")
                if freq is None or np.isnan(freq):
                    freq = np.nan
                results[f'F{formant_num}'].append(freq)

                bw = call(formant, "Get bandwidth at time", formant_num, t, "Hertz", "Linear")
                if bw is None or np.isnan(bw):
                    bw = np.nan
                results[f'B{formant_num}'].append(bw)

        # Extract spectral measures
        if spectral_measures and pitch is not None and formant is not None:
            spectral_vals = _extract_spectral_measures(
                sound, pitch, formant, t,
                window_size / 1000.0,
                hnr_objects,
                use_bandwidth_formula
            )
            for key in spectral_keys:
                results[key].append(spectral_vals.get(key, np.nan))

    return results


def _extract_spectral_measures(sound, pitch, formant, t, window_length, hnr_objects, use_bw_formula):
    """
    Extract spectral measures at a given timepoint.

    This includes harmonics (H1, H2, H4, etc.), formant amplitudes (A1-A3),
    CPP, and HNR measures.
    """
    measures = {}

    # Get F0 at this timepoint
    f0 = call(pitch, "Get value at time", t, "Hertz", "Linear")

    if f0 is None or np.isnan(f0) or f0 <= 0:
        # No voicing - return NaN for all measures
        for key in ['H1u', 'H2u', 'H4u', 'H2Ku', 'H5Ku', 'A1u', 'A2u', 'A3u',
                    'H1H2u', 'H2H4u', 'H1A1u', 'H1A2u', 'H1A3u', 'H2KH5Ku',
                    'H1c', 'H2c', 'H4c', 'A1c', 'A2c', 'A3c',
                    'H1H2c', 'H2H4c', 'H1A1c', 'H1A2c', 'H1A3c',
                    'CPP', 'HNR05', 'HNR15', 'HNR25', 'HNR35']:
            measures[key] = np.nan
        return measures

    # Get formant frequencies
    F1 = call(formant, "Get value at time", 1, t, "Hertz", "Linear")
    F2 = call(formant, "Get value at time", 2, t, "Hertz", "Linear")
    F3 = call(formant, "Get value at time", 3, t, "Hertz", "Linear")

    # Get formant bandwidths
    B1 = call(formant, "Get bandwidth at time", 1, t, "Hertz", "Linear")
    B2 = call(formant, "Get bandwidth at time", 2, t, "Hertz", "Linear")
    B3 = call(formant, "Get bandwidth at time", 3, t, "Hertz", "Linear")

    # Create LPC spectrum at this timepoint for harmonic amplitude extraction
    # Extract a portion of sound around timepoint
    t_start = max(0, t - window_length / 2)
    t_end = min(call(sound, "Get total duration"), t + window_length / 2)

    try:
        sound_part = call(sound, "Extract part", t_start, t_end, "rectangular", 1.0, False)
        spectrum = call(sound_part, "To Spectrum", True)  # FFT
        ltas = call(sound_part, "To Ltas", 1)
    except:
        # If extraction fails, return NaN
        for key in ['H1u', 'H2u', 'H4u', 'H2Ku', 'H5Ku', 'A1u', 'A2u', 'A3u',
                    'H1H2u', 'H2H4u', 'H1A1u', 'H1A2u', 'H1A3u', 'H2KH5Ku',
                    'H1c', 'H2c', 'H4c', 'A1c', 'A2c', 'A3c',
                    'H1H2c', 'H2H4c', 'H1A1c', 'H1A2c', 'H1A3c',
                    'CPP', 'HNR05', 'HNR15', 'HNR25', 'HNR35']:
            measures[key] = np.nan
        return measures

    # Extract uncorrected harmonic amplitudes
    H1u = _get_harmonic_amplitude(spectrum, f0 * 1)
    H2u = _get_harmonic_amplitude(spectrum, f0 * 2)
    H4u = _get_harmonic_amplitude(spectrum, f0 * 4)
    H2Ku = _get_harmonic_amplitude(spectrum, 2000)
    H5Ku = _get_harmonic_amplitude(spectrum, 5000)

    # Extract uncorrected formant amplitudes
    A1u = _get_harmonic_amplitude(spectrum, F1) if not np.isnan(F1) else np.nan
    A2u = _get_harmonic_amplitude(spectrum, F2) if not np.isnan(F2) else np.nan
    A3u = _get_harmonic_amplitude(spectrum, F3) if not np.isnan(F3) else np.nan

    # Store uncorrected values
    measures['H1u'] = H1u
    measures['H2u'] = H2u
    measures['H4u'] = H4u
    measures['H2Ku'] = H2Ku
    measures['H5Ku'] = H5Ku
    measures['A1u'] = A1u
    measures['A2u'] = A2u
    measures['A3u'] = A3u

    # Calculate uncorrected differences
    measures['H1H2u'] = H1u - H2u
    measures['H2H4u'] = H2u - H4u
    measures['H1A1u'] = H1u - A1u if not np.isnan(A1u) else np.nan
    measures['H1A2u'] = H1u - A2u if not np.isnan(A2u) else np.nan
    measures['H1A3u'] = H1u - A3u if not np.isnan(A3u) else np.nan
    measures['H2KH5Ku'] = H2Ku - H5Ku

    # Apply Iseli corrections for formant influence
    # Simplified version - full correction requires complex formant bandwidth modeling
    # For now, we'll compute basic corrections
    if not np.isnan(F1) and not np.isnan(B1):
        H1c = H1u - _formant_correction(f0 * 1, F1, B1)
        A1c = A1u - _formant_correction(F1, F1, B1)
    else:
        H1c = H1u
        A1c = A1u

    if not np.isnan(F2) and not np.isnan(B2):
        H2c = H2u - _formant_correction(f0 * 2, F2, B2)
        A2c = A2u - _formant_correction(F2, F2, B2)
    else:
        H2c = H2u
        A2c = A2u

    if not np.isnan(F3) and not np.isnan(B3):
        H4c = H4u - _formant_correction(f0 * 4, F3, B3)
        A3c = A3u - _formant_correction(F3, F3, B3)
    else:
        H4c = H4u
        A3c = A3u

    # Store corrected values
    measures['H1c'] = H1c
    measures['H2c'] = H2c
    measures['H4c'] = H4c
    measures['A1c'] = A1c
    measures['A2c'] = A2c
    measures['A3c'] = A3c

    # Calculate corrected differences
    measures['H1H2c'] = H1c - H2c
    measures['H2H4c'] = H2c - H4c
    measures['H1A1c'] = H1c - A1c if not np.isnan(A1c) else np.nan
    measures['H1A2c'] = H1c - A2c if not np.isnan(A2c) else np.nan
    measures['H1A3c'] = H1c - A3c if not np.isnan(A3c) else np.nan

    # Calculate CPP (Cepstral Peak Prominence)
    try:
        power_cepstrogram = call(sound_part, "To PowerCepstrogram", 60, 0.002, 5000, 50)
        cpp = call(power_cepstrogram, "Get CPPS", False, 0.01, 0.001, 60, 330, 0.05,
                   "Parabolic", 0.001, 0, "Straight", "Robust")
        measures['CPP'] = cpp
    except:
        measures['CPP'] = np.nan

    # Extract HNR measures
    for band in ['05', '15', '25', '35']:
        key = f'HNR{band}'
        if key in hnr_objects:
            hnr_val = call(hnr_objects[key], "Get value at time", t, "Linear")
            measures[key] = hnr_val if not np.isnan(hnr_val) else np.nan
        else:
            measures[key] = np.nan

    return measures


def _get_harmonic_amplitude(spectrum, frequency):
    """
    Get amplitude at a given frequency from spectrum (in dB).
    """
    try:
        # Get amplitude in dB
        amplitude = call(spectrum, "Get real value at frequency", frequency, "Sinc70")
        if amplitude is None or np.isnan(amplitude):
            return np.nan

        # Convert to dB (spectrum amplitudes are linear)
        if amplitude > 0:
            db = 20 * np.log10(amplitude)
        else:
            db = -100  # Very small value

        return db
    except:
        return np.nan


def _formant_correction(harmonic_freq, formant_freq, formant_bw):
    """
    Calculate formant influence correction (simplified Iseli method).

    This is a simplified version. Full correction requires modeling
    the formant as a resonance filter.
    """
    if np.isnan(formant_freq) or np.isnan(formant_bw):
        return 0.0

    # Calculate distance from harmonic to formant
    distance = abs(harmonic_freq - formant_freq)

    # Simplified correction based on formant bandwidth
    # Full Iseli correction uses complex formant response calculation
    if distance < formant_bw:
        # Close to formant - significant boost
        correction = 10 * np.log10(1 + (formant_bw / max(distance, 1)))
    elif distance < formant_bw * 2:
        # Moderate distance
        correction = 5 * np.log10(1 + (formant_bw / distance))
    else:
        # Far from formant - minimal influence
        correction = 0.0

    return correction

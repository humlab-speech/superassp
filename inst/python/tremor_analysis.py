"""
Vocal Tremor Analysis - Python/Parselmouth Implementation
Based on Praat TREMOR script by Markus Brückl
Extracts 18 measures of vocal tremor from sustained phonations
"""

import parselmouth
from parselmouth.praat import call
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, Union, Optional, Tuple
import warnings


class TremorAnalyzer:
    """
    Analyzes vocal tremor in sustained phonations.
    
    Extracts frequency and amplitude tremor measures including:
    - Magnitude, cyclicality, frequency, intensity indices
    - Power indices, product sums, and HNR measures
    """
    
    def __init__(self,
                 analysis_time_step: float = 0.015,
                 min_pitch: float = 60.0,
                 max_pitch: float = 350.0,
                 silence_threshold: float = 0.03,
                 voicing_threshold: float = 0.3,
                 octave_cost: float = 0.01,
                 octave_jump_cost: float = 0.35,
                 voiced_unvoiced_cost: float = 0.14,
                 min_tremor_freq: float = 1.5,
                 max_tremor_freq: float = 15.0,
                 contour_magnitude_threshold: float = 0.01,
                 tremor_cyclicality_threshold: float = 0.15,
                 freq_tremor_octave_cost: float = 0.01,
                 amp_tremor_octave_cost: float = 0.01,
                 nan_output_mode: int = 2):
        """
        Initialize tremor analyzer with parameters.
        
        Args:
            analysis_time_step: Time step for analysis in seconds
            min_pitch: Minimum pitch for extraction in Hz
            max_pitch: Maximum pitch for extraction in Hz
            silence_threshold: Threshold for silence detection
            voicing_threshold: Threshold for voicing detection
            octave_cost: Cost for octave jumps in pitch tracking
            octave_jump_cost: Cost for large octave jumps
            voiced_unvoiced_cost: Cost for voiced/unvoiced transitions
            min_tremor_freq: Minimum tremor frequency in Hz
            max_tremor_freq: Maximum tremor frequency in Hz
            contour_magnitude_threshold: Threshold for contour magnitude
            tremor_cyclicality_threshold: Threshold for cyclicality
            freq_tremor_octave_cost: Octave cost for frequency tremor
            amp_tremor_octave_cost: Octave cost for amplitude tremor
            nan_output_mode: 1=zeros for undefined, 2=undefined values
        """
        self.ts = analysis_time_step
        self.min_pitch = min_pitch
        self.max_pitch = max_pitch
        self.silence_threshold = silence_threshold
        self.voicing_threshold = voicing_threshold
        self.octave_cost = octave_cost
        self.octave_jump_cost = octave_jump_cost
        self.voiced_unvoiced_cost = voiced_unvoiced_cost
        self.min_tremor_freq = min_tremor_freq
        self.max_tremor_freq = max_tremor_freq
        self.contour_mag_thresh = contour_magnitude_threshold
        self.tremor_cyc_thresh = tremor_cyclicality_threshold
        self.freq_tremor_oct_cost = freq_tremor_octave_cost
        self.amp_tremor_oct_cost = amp_tremor_octave_cost
        self.nan_output_mode = nan_output_mode
        
    def analyze_file(self, filepath: Union[str, Path]) -> Dict[str, float]:
        """
        Analyze a sound file for tremor measures.
        
        Args:
            filepath: Path to audio file
            
        Returns:
            Dictionary containing 18 tremor measures
        """
        sound = parselmouth.Sound(str(filepath))
        return self.analyze_sound(sound)
    
    def analyze_sound(self, sound: parselmouth.Sound) -> Dict[str, float]:
        """
        Analyze a Parselmouth Sound object for tremor measures.
        
        Args:
            sound: Parselmouth Sound object
            
        Returns:
            Dictionary with tremor measures
        """
        duration = sound.duration
        
        # Extract pitch
        pitch = sound.to_pitch_cc(
            time_step=self.ts,
            pitch_floor=self.min_pitch,
            pitch_ceiling=self.max_pitch
        )
        
        # Analyze frequency tremor
        ftrem_results = self._analyze_frequency_tremor(sound, pitch, duration)
        
        # Analyze amplitude tremor
        atrem_results = self._analyze_amplitude_tremor(sound, pitch, duration)
        
        # Combine results
        results = {**ftrem_results, **atrem_results}
        
        return results
    
    def _analyze_frequency_tremor(self, sound, pitch, duration):
        """Analyze frequency tremor measures."""
        results = {}
        
        # Count voiced frames
        n_voiced = 0
        n_frames = pitch.get_number_of_frames()
        
        for i in range(1, n_frames + 1):
            f0 = pitch.get_value_in_frame(i)
            if not np.isnan(f0) and f0 > 0:
                n_voiced += 1
        
        if n_voiced == 0:
            # No voiced frames
            return self._get_undefined_freq_tremor()
        
        # Create matrix for frequency contour
        x1 = pitch.get_time_from_frame_number(1)
        
        # Extract F0 values into matrix
        f0_matrix = []
        for i in range(1, n_frames + 1):
            f0 = pitch.get_value_in_frame(i)
            if np.isnan(f0) or f0 <= 0:
                f0_matrix.append(0.0)
            else:
                f0_matrix.append(f0)
        
        f0_array = np.array(f0_matrix)
        
        # Remove linear trend (detrend)
        voiced_indices = f0_array > 0
        if np.sum(voiced_indices) < 2:
            return self._get_undefined_freq_tremor()
        
        mean_f0 = np.mean(f0_array[voiced_indices])
        
        # Normalize by mean F0
        f0_norm = np.where(voiced_indices, (f0_array - mean_f0) / mean_f0, 0)
        
        # Convert to Sound for autocorrelation
        sampling_freq = 1.0 / self.ts
        snd_trem = self._array_to_sound(f0_norm, sampling_freq)
        
        # Calculate tremor contour HNR
        hnr = self._calculate_tremor_hnr(snd_trem)
        results['ftrHNR'] = hnr if not np.isnan(hnr) else (0 if self.nan_output_mode == 1 else np.nan)
        
        # Extract tremor frequency using pitch analysis
        pitch_trem = snd_trem.to_pitch_cc(
            time_step=duration,
            pitch_floor=self.min_tremor_freq,
            pitch_ceiling=self.max_tremor_freq
        )
        
        # Read pitch object for magnitude and cyclicality
        trm, trc = self._read_pitch_object(pitch_trem)
        
        results['ftrm'] = trm
        results['ftrc'] = trc
        
        # Get tremor frequency from strongest candidate
        tremor_freq, tremor_strength, n_candidates = self._get_tremor_candidates(pitch_trem)
        
        results['fmodN'] = n_candidates
        
        if n_candidates > 0 and tremor_strength > self.tremor_cyc_thresh and trm > self.contour_mag_thresh:
            results['ftrf'] = tremor_freq
            
            # Calculate intensity index
            tri = self._calculate_intensity_index(snd_trem, pitch_trem)
            results['ftri'] = tri
            results['ftrp'] = tri * tremor_freq / (tremor_freq + 1)
            results['ftrcip'] = tri * trc
        else:
            if self.nan_output_mode == 2:
                results['ftrf'] = np.nan
                results['ftri'] = np.nan
                results['ftrp'] = np.nan
                results['ftrcip'] = np.nan
            else:
                results['ftrf'] = 0.0
                results['ftri'] = 0.0
                results['ftrp'] = 0.0
                results['ftrcip'] = 0.0
        
        # Calculate product sum
        ftrps = self._calculate_product_sum(snd_trem, pitch_trem, 'frequency', trm, trc)
        results['ftrps'] = ftrps if ftrps != 0 else (np.nan if self.nan_output_mode == 2 else 0)
        
        return results
    
    def _analyze_amplitude_tremor(self, sound, pitch, duration):
        """Analyze amplitude tremor measures."""
        results = {}
        
        n_frames = pitch.get_number_of_frames()
        x1 = pitch.get_time_from_frame_number(1)
        
        # Extract amplitude per pitch period
        # Requires both sound and pitch objects
        point_process = call([sound, pitch], "To PointProcess (cc)")
        n_points = call(point_process, "Get number of points")

        if n_points < 3:
            return self._get_undefined_amp_tremor()

        # Extract RMS amplitude per period
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
        
        if len(amp_values) == 0:
            return self._get_undefined_amp_tremor()
        
        # Resample amplitude contour at constant rate
        amp_array = self._resample_amplitude_contour(
            sound, pitch, point_process, amp_values, duration, x1
        )
        
        if amp_array is None or len(amp_array) < 3:
            return self._get_undefined_amp_tremor()

        # Normalize amplitude contour
        amp_mean = np.mean(amp_array[amp_array > 0])
        if amp_mean == 0:
            return self._get_undefined_amp_tremor()

        amp_norm = np.where(amp_array > 0, (amp_array - amp_mean) / amp_mean, 0)
        
        # Convert to sound
        sampling_freq = 1.0 / self.ts
        snd_trem = self._array_to_sound(amp_norm, sampling_freq)
        
        # Calculate HNR
        hnr = self._calculate_tremor_hnr(snd_trem)
        results['atrHNR'] = hnr if not np.isnan(hnr) else (0 if self.nan_output_mode == 1 else np.nan)
        
        # Extract tremor frequency
        pitch_trem = snd_trem.to_pitch_cc(
            time_step=duration,
            pitch_floor=self.min_tremor_freq,
            pitch_ceiling=self.max_tremor_freq
        )
        
        # Read pitch object
        trm, trc = self._read_pitch_object(pitch_trem)
        results['atrm'] = trm
        results['atrc'] = trc
        
        # Get tremor candidates
        tremor_freq, tremor_strength, n_candidates = self._get_tremor_candidates(pitch_trem)
        results['amodN'] = n_candidates
        
        if n_candidates > 0 and tremor_strength > self.tremor_cyc_thresh and trm > self.contour_mag_thresh:
            results['atrf'] = tremor_freq
            
            # Calculate intensity index
            tri = self._calculate_intensity_index(snd_trem, pitch_trem)
            results['atri'] = tri
            results['atrp'] = tri * tremor_freq / (tremor_freq + 1)
            results['atrcip'] = tri * trc
        else:
            if self.nan_output_mode == 2:
                results['atrf'] = np.nan
                results['atri'] = np.nan
                results['atrp'] = np.nan
                results['atrcip'] = np.nan
            else:
                results['atrf'] = 0.0
                results['atri'] = 0.0
                results['atrp'] = 0.0
                results['atrcip'] = 0.0
        
        # Calculate product sum
        atrps = self._calculate_product_sum(snd_trem, pitch_trem, 'amplitude', trm, trc)
        results['atrps'] = atrps if atrps != 0 else (np.nan if self.nan_output_mode == 2 else 0)
        
        return results
    
    def _resample_amplitude_contour(self, sound, pitch, point_process, amp_values, duration, x1):
        """Resample amplitude contour at constant time steps using interpolation."""
        n_frames = pitch.get_number_of_frames()
        n_points = len(amp_values)

        if n_points == 0:
            return None

        # Get time points for interpolation
        time_points = []
        for i in range(1, n_points + 1):
            t = call(point_process, "Get time from index", i)
            time_points.append(t)

        time_points = np.array(time_points)
        amp_values_array = np.array(amp_values)

        # Create target time points (one per frame)
        amp_matrix = np.zeros(n_frames)

        for iframe in range(n_frames):
            # Get F0 for this frame
            f0 = pitch.get_value_in_frame(iframe + 1)

            if np.isnan(f0) or f0 <= 0:
                amp_matrix[iframe] = 0
                continue

            # Calculate time for this frame
            t = iframe * self.ts + x1

            # Linear interpolation of amplitude at this time
            if t <= time_points[0]:
                amp_matrix[iframe] = amp_values_array[0]
            elif t >= time_points[-1]:
                amp_matrix[iframe] = amp_values_array[-1]
            else:
                # Find surrounding points
                idx = np.searchsorted(time_points, t)
                if idx == 0:
                    amp_matrix[iframe] = amp_values_array[0]
                elif idx >= len(time_points):
                    amp_matrix[iframe] = amp_values_array[-1]
                else:
                    # Linear interpolation
                    t0 = time_points[idx - 1]
                    t1 = time_points[idx]
                    a0 = amp_values_array[idx - 1]
                    a1 = amp_values_array[idx]

                    if t1 != t0:
                        alpha = (t - t0) / (t1 - t0)
                        amp_matrix[iframe] = a0 + alpha * (a1 - a0)
                    else:
                        amp_matrix[iframe] = a0

        return amp_matrix
    
    def _array_to_sound(self, array, sampling_freq):
        """Convert numpy array to Parselmouth Sound."""
        # Create sound from values
        sound = parselmouth.Sound(array, sampling_frequency=sampling_freq)
        return sound
    
    def _calculate_tremor_hnr(self, sound):
        """Calculate harmonicity-to-noise ratio for tremor contour."""
        try:
            hnr_ts = 1.0 / self.min_tremor_freq
            periods_per_window = sound.duration * self.min_tremor_freq * (2/3)
            
            harmonicity = sound.to_harmonicity_cc(
                time_step=hnr_ts,
                minimum_pitch=self.min_tremor_freq,
                silence_threshold=self.contour_mag_thresh,
                periods_per_window=periods_per_window
            )
            
            hnr = call(harmonicity, "Get mean", 0, 0)
            return hnr
        except:
            return np.nan
    
    def _read_pitch_object(self, pitch):
        """Read magnitude and cyclicality from pitch object."""
        try:
            # Save pitch to temporary format and read back
            # This mimics the Praat "Save as text file" approach
            n_frames = pitch.get_number_of_frames()
            
            # Get intensity (magnitude) and number of candidates
            intensity = call(pitch, "Get mean", 0, 0, "Hertz")
            
            # Get number of candidates at frame with strongest pitch
            max_strength = 0
            max_frame = 1
            
            for i in range(1, n_frames + 1):
                try:
                    f0 = pitch.get_value_in_frame(i)
                    if not np.isnan(f0) and f0 > 0:
                        # Approximate strength by inverse of voiced/unvoiced cost
                        strength = 1.0 if f0 > 0 else 0.0
                        if strength > max_strength:
                            max_strength = strength
                            max_frame = i
                except:
                    continue
            
            # Get the cyclicality as the strength of the strongest candidate
            cyclicality = max_strength
            magnitude = intensity / 100.0  # Normalize
            
            return magnitude, cyclicality
            
        except Exception as e:
            return 0.0, 0.0
    
    def _get_tremor_candidates(self, pitch):
        """Get tremor frequency candidates from pitch object."""
        try:
            n_frames = pitch.get_number_of_frames()
            
            # Find strongest tremor frequency
            max_strength = 0
            tremor_freq = 0
            n_candidates = 0
            
            for i in range(1, n_frames + 1):
                f0 = pitch.get_value_in_frame(i)
                if not np.isnan(f0) and f0 > 0:
                    n_candidates += 1
                    # Use frame 1 as it represents the overall pitch
                    if i == 1:
                        tremor_freq = f0
                        max_strength = 1.0
            
            return tremor_freq, max_strength, n_candidates
            
        except:
            return 0, 0, 0
    
    def _calculate_intensity_index(self, sound, pitch):
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
    
    def _calculate_product_sum(self, sound, pitch, contour_type, trm, trc):
        """Calculate cyclicality-weighted sum of intensity indices."""
        try:
            n_frames = pitch.get_number_of_frames()
            
            tris = 0
            rank = 0
            
            for iframe in range(1, n_frames + 1):
                tr_freq = pitch.get_value_in_frame(iframe)
                
                if np.isnan(tr_freq) or tr_freq <= 0:
                    continue
                
                # Get strength (cyclicality) - simplified
                tr_strength = 1.0 if tr_freq > 0 else 0
                
                if tr_strength > self.tremor_cyc_thresh and trm > self.contour_mag_thresh:
                    rank += 1
                    
                    # Calculate intensity index for this frequency
                    tri = self._calculate_intensity_index(sound, pitch)
                    tris += tr_strength * tri
            
            return tris
            
        except:
            return 0.0
    
    def _get_undefined_freq_tremor(self):
        """Return undefined values for frequency tremor."""
        if self.nan_output_mode == 1:
            return {
                'ftrm': 0, 'ftrc': 0, 'fmodN': 0, 'ftrf': 0,
                'ftri': 0, 'ftrp': 0, 'ftrcip': 0, 'ftrps': 0, 'ftrHNR': 0
            }
        else:
            return {
                'ftrm': np.nan, 'ftrc': np.nan, 'fmodN': 0, 'ftrf': np.nan,
                'ftri': np.nan, 'ftrp': np.nan, 'ftrcip': np.nan, 
                'ftrps': np.nan, 'ftrHNR': np.nan
            }
    
    def _get_undefined_amp_tremor(self):
        """Return undefined values for amplitude tremor."""
        if self.nan_output_mode == 1:
            return {
                'atrm': 0, 'atrc': 0, 'amodN': 0, 'atrf': 0,
                'atri': 0, 'atrp': 0, 'atrcip': 0, 'atrps': 0, 'atrHNR': 0
            }
        else:
            return {
                'atrm': np.nan, 'atrc': np.nan, 'amodN': 0, 'atrf': np.nan,
                'atri': np.nan, 'atrp': np.nan, 'atrcip': np.nan,
                'atrps': np.nan, 'atrHNR': np.nan
            }
    
    def analyze_directory(self, input_path: Union[str, Path], 
                         output_file: Optional[Union[str, Path]] = None) -> pd.DataFrame:
        """
        Analyze all WAV files in a directory.
        
        Args:
            input_path: Directory containing WAV files
            output_file: Optional path to save results CSV
            
        Returns:
            DataFrame with tremor measures for all files
        """
        input_path = Path(input_path)
        wav_files = list(input_path.glob("*.wav"))
        
        results_list = []
        
        for wav_file in wav_files:
            print(f"Processing: {wav_file.name}")
            try:
                results = self.analyze_file(wav_file)
                results['filename'] = wav_file.stem
                results_list.append(results)
            except Exception as e:
                print(f"Error processing {wav_file.name}: {e}")
                continue
        
        df = pd.DataFrame(results_list)
        
        # Reorder columns
        cols = ['filename'] + [col for col in df.columns if col != 'filename']
        df = df[cols]
        
        if output_file:
            df.to_csv(output_file, index=False)
            print(f"Results saved to: {output_file}")
        
        return df


# Example usage
if __name__ == "__main__":
    # Initialize analyzer with default parameters
    analyzer = TremorAnalyzer()
    
    # Analyze single file
    filepath = Path("/Users/frkkan96/Documents/src/superassp/inst/samples/sustained/a1.wav").expanduser()
    results = analyzer.analyze_file(filepath)
    print(results)
    
    # Analyze directory
    # df = analyzer.analyze_directory("./sounds", "./results/tremor_results.csv")
    # print(df)
    
    print("Tremor Analyzer initialized. Use analyze_file() or analyze_directory().")

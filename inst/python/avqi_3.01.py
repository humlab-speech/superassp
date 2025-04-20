import parselmouth
from parselmouth.praat import call
import os

def parselmouth_avqi(
    Simple_version=True,
    name_patient="Fredrik Karlsson",
    Date_of_birth="1975-12-31",
    Assessment_date="2021-12-31",
    Input_directory="/Users/frkkan96/Documents/src/superassp/tests/signalfiles/AVQI/input/",
    Generate_PDF_files=True,
    Speaker_ID="1",
    Output_directory="/Users/frkkan96/Documents/src/superassp/tests/signalfiles/AVQI/output",
    Output_file="/Users/frkkan96/Documents/src/superassp/tests/signalfiles/AVQI/output/parselmouth_avqi_3_01.csv"
):
    # Helper function to get file list
    def get_file_list(directory, pattern):
        return [os.path.join(directory, f) for f in os.listdir(directory) if f.startswith(pattern)]

    # Load and concatenate sustained vowels
    sv_files = get_file_list(Input_directory, 'sv')
    sv_sounds = [parselmouth.Sound(f) for f in sv_files]
    sv = call(sv_sounds, "Concatenate")

    # Load and concatenate continuous speech files
    cs_files = get_file_list(Input_directory, 'cs')
    cs_sounds = [parselmouth.Sound(f) for f in cs_files]
    sv = call(cs_sounds, "Concatenate")

    # High-pass filtering of the sound files
    cs2 = call(cs, "Filter (stop Hann band)", 0, 34, 0.1)
    sv2 = call(sv, "Filter (stop Hann band)", 0, 34, 0.1)

    # Detection, extraction and concatenation of the voiced segments in the recording of continuous speech
    original = cs2.copy()
    sampling_rate = call(original, "Get sampling frequency")
    intermediate_samples = call(original, "Get sampling period")
    
    only_voice = call("Create Sound", "onlyVoice", 0, 0.001, sampling_rate,"0")
    
    textgrid = call(original, "To TextGrid (silences)", 50, 0.003, -25, 0.1, 0.1, "silence", "sounding")
    
    intervals = call([original, textgrid], "Extract intervals where", 1, False, "does not contain", "silence")
    
    only_loud= call(intervals, "Concatenate")
    
    signal_end = call(only_loud,"Get end time")
    window_border_left = call(only_loud,"Get start time")
    window_width = 0.03
    extreme_right = signal_end - window_width
      
    global_power = call(only_loud, "Get power in air")
    
    voiceless_threshold = global_power * (30 / 100)
    

    while window_border_left < extreme_right:
   
        part = call(only_loud, "Extract part", window_border_left, window_border_right, "rectangular", 1.0, False)
        partial_power = call(part, "Get power in air")
        
        if partial_power > voiceless_threshold:
            zero_crossing_rate = call(part, "Get number of zero crossings") / (window_border_right - window_border_left)
            if zero_crossing_rate < 3000:
                only_voice.append(part)
        
        window_border_left += window_width
    
    # Determination of the six acoustic measures and calculation of the Acoustic Voice Quality Index
    duration_vowel = call(sv2, "Get total duration")
    
    if duration_vowel > 3:
        sv3 = call(sv2, "Extract part", duration_vowel - 3, duration_vowel, "rectangular", 1.0, False)
    else:
        sv3 = sv2.copy()
    
    avqi_sound = only_voice.concatenate(sv3)
    
    cpps = call(avqi_sound.to_power_cepstrogram(), "Get CPPS", False, 0.01, 0.001, 60, 330, 0.05)
    
    ltas_slope = call(avqi_sound.to_ltas(), "Get slope", 0, 1000, 1000, 10000)
    
    ltas_tilt = call(avqi_sound.to_ltas(), "Compute trend line", 1, 10000).get_slope(0, 1000)
    
    shimmer_local = call(avqi_sound.to_point_process_periodic_cc(), "Get shimmer (local)", 0.0001, 0.02)
    
    shimmer_local_db = call(avqi_sound.to_point_process_periodic_cc(), "Get shimmer (local_dB)", 0.0001, 0.02)
    
    hnr = call(avqi_sound.to_pitch_cc(), "Voice report", 75).extract_number("Mean harmonics-to-noise ratio:")
    
    avqi_score = (4.152 - (0.177 * cpps) - (0.006 * hnr) - (0.037 * shimmer_local) + 
                  (0.941 * shimmer_local_db) + (0.01 * ltas_slope) + (0.093 * ltas_tilt)) * 2.8902
    
    # Create output dictionary
    outTab = {
        "AVQI_VERSION": "v03.01",
        "Speaker": name_patient,
        "ID": Speaker_ID,
        "CPPS": round(cpps, 2),
        "HNR": round(hnr, 2),
        "Shim_local": round(shimmer_local * 100, 2),
        "Shim_local_DB": round(shimmer_local_db, 2),
        "LTAS_Slope": round(ltas_slope, 2),
        "LTAS_Tilt": round(ltas_tilt, 2),
        "AVQI": round(avqi_score, 2)
    }
    
    return outTab

# Example usage
result = parselmouth_avqi()
print(result)

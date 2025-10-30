import parselmouth as pm
import math
import pandas as pd
import numpy as np
import io as io



def praat_formant_burg(
    sound,  # Changed from soundFile - now accepts Sound object
    beginTime=0.0,
    endTime=0.0,
    timeStep = 0.005,
    number_of_formants=5.0,
    maxHzFormant=5500.0,
    windowLength=0.025,
    pre_emphasis=50.0,
    track_formants=False,
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
    spectrogram_window_shape="Gaussian",  # Accept string, will convert to enum
    spectrogram_resolution=40.0):

    # Accept Sound object directly (in-memory)
    snd = sound
    dur = snd.get_total_duration()
    if  beginTime > 0.0 and endTime > 0.0 and beginTime >= 0.0 and endTime <= dur :
        snd = snd.extract_part(beginTime, endTime, windowShape, relativeWidth, True)

    form = snd.to_formant_burg(time_step=timeStep,
        max_number_of_formants=number_of_formants,
        maximum_formant=maxHzFormant,
        window_length=windowLength,
        pre_emphasis_from=pre_emphasis)

    if track_formants:
        form = pm.praat.call(form,"Track",number_of_tracks,reference_F1, reference_F2, reference_F3, reference_F4, reference_F5, frequency_cost, bandwidth_cost, transition_cost)
        number_of_formants = min(number_of_formants, number_of_tracks)

    formantTable = pm.praat.call(form,"Down to Table", True, True, 10, True, 3, True, 3, True)
    nFormantRows = pm.praat.call(formantTable,"Get number of rows")
    for formant in range(1,math.ceil(number_of_formants) + 1):
        lab = f"L{formant}(dB)"
        pm.praat.call(formantTable,"Append column",lab)

    maxHz = maxHzFormant + 2000.0

    # Convert string to enum if needed (Parselmouth requires enum)
    if isinstance(spectrogram_window_shape, str):
        if spectrogram_window_shape.lower() == "gaussian":
            spec_window = pm.SpectralAnalysisWindowShape.GAUSSIAN
        else:
            spec_window = pm.SpectralAnalysisWindowShape.GAUSSIAN  # Default
    else:
        spec_window = spectrogram_window_shape

    spec = snd.to_spectrogram(windowLength, maxHz, timeStep, spectrogram_resolution, spec_window)
    
    for r in range(1,nFormantRows + 1):
        for f in range(1,math.ceil(number_of_formants) + 1):
            outLab = f"L{f}(dB)"
            inLab = f"F{f}(Hz)"
            currFreq = pm.praat.call(formantTable,"Get value",r,inLab)
            if currFreq != "--undefined--":
                currTime = pm.praat.call(formantTable,"Get value",r,"time(s)")
                currIntensity = pm.praat.call(spec, "Get power at", float(currTime), float(currFreq))
                pm.praat.call(formantTable,"Set numeric value",r, outLab,currIntensity)

    return pd.read_table(io.StringIO(pm.praat.call(formantTable, "List", True)))

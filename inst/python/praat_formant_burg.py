soundFile = "/Users/frkkan96/Desktop/a1.wav"

import parselmouth as pm
import math
import pandas as pd
import numpy as np
import io as io

#pt = pm.praat.call("Create formant table (Peterson & Barney 1952)")
def PraatTableToPandas(pt):
    ncols = pm.praat.call(pt,"Get number of columns") #Note that +1 needs to be used on the lext line
    columnNames = [pm.praat.call(pt,"Get column label", i) for i in range(1,ncols + 1)] 
    outDataFrame = pd.DataFrame(columns=columnNames)
    nrows = pm.praat.call(pt,"Get number of rows")
    for col in columnNames:
        rowdata = [ pm.praat.call(pt,"Get value",i, col) for i in range(1,nrows)]
        outDataFrame[col] = rowdata 
    return outDataFrame.replace('--undefined--',np.nan).replace('',np.nan)

def PraatTableToPandas2(pt):
    return pd.read_table(io.StringIO(pm.praat.call(pt, "List", True)))



#pb = PraatTableToPandas(pt)

def praat_formant_burg(
    soundFile,
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
    spectrogram_window_shape= pm.SpectralAnalysisWindowShape.GAUSSIAN,
    spectrogram_resolution=40.0):

    snd = pm.Sound(soundFile)
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
    for formant in range(1,math.ceil(number_of_formants)):
        lab = f"L{formant}(dB)"
        pm.praat.call(formantTable,"Append column",lab)

    maxHz = maxHzFormant + 2000.0

        
    spec = snd.to_spectrogram(windowLength, maxHz, timeStep, spectrogram_resolution,  spectrogram_window_shape )
    
    for r in range(1,nFormantRows):
        for f in range(1,math.ceil(number_of_formants)):
            outLab = f"L{f}(dB)"
            inLab = f"F{f}(Hz)"
            currFreq = pm.praat.call(formantTable,"Get value",r,inLab)
            if currFreq != "--undefined--":
                currTime = pm.praat.call(formantTable,"Get value",r,"time(s)")
                currIntensity = pm.praat.call(spec, "Get power at", float(currTime), float(currFreq))
                pm.praat.call(formantTable,"Set numeric value",r, outLab,currIntensity)

    return pd.read_table(io.StringIO(pm.praat.call(formantTable, "List", True)))

res = praat_formant_burg(soundFile)


import parselmouth as pm
import math
import pandas as pd
import numpy as np
import io as io

# #pt = pm.praat.call("Create formant table (Peterson & Barney 1952)")
# def PraatTableToPandas(pt):
#     ncols = pm.praat.call(pt,"Get number of columns") #Note that +1 needs to be used on the lext line
#     columnNames = [pm.praat.call(pt,"Get column label", i) for i in range(1,ncols + 1)] 
#     outDataFrame = pd.DataFrame(columns=columnNames)
#     nrows = pm.praat.call(pt,"Get number of rows")
#     for col in columnNames:
#         rowdata = [ pm.praat.call(pt,"Get value",i, col) for i in range(1,nrows)]
#         outDataFrame[col] = rowdata 
#     return outDataFrame.replace('--undefined--',np.nan).replace('',np.nan)

# def PraatTableToPandas2(pt):
#     return pd.read_table(io.StringIO(pm.praat.call(pt, "List", True)))



#pb = PraatTableToPandas(pt)

def praat_voice_report(
    soundFile = "/Users/frkkan96/Desktop/a1.wav",
    beginTime=0.0,
    endTime=0.0,
    selectionOffset=0.0,
    selectionLength=2.0,
    timeStep = None,
    minF= 75.0,
    maxF=600.0,
    max_number_of_candidates=15,
    maximum_period_factor=1.3,
    maximum_amplitude_factor=1.6,
    silence_threshold=0.03,
    voicing_threshold=0.45,
    octave_cost=0.01,
    octave_jump_cost=0.35,
    voiced_unvoiced_cost=0.14,
    windowShape=pm.WindowShape.GAUSSIAN1,
    relativeWidth=1.0):

    snd = pm.Sound(soundFile)
    dur = snd.get_total_duration()
    if  beginTime > 0.0 and endTime > 0.0 and beginTime >= 0.0 and endTime <= dur :
        snd = snd.extract_part(beginTime, endTime, windowShape, relativeWidth, True)

    currPP=  pm.praat.call(snd,"To PointProcess (periodic, cc)",minF,maxF)
    currPitch = snd.to_pitch_cc(time_step=timeStep,
        pitch_floor=minF,
        max_number_of_candidates=int(max_number_of_candidates),
        very_accurate=True,
        silence_threshold=silence_threshold,
        voicing_threshold=voicing_threshold,
        octave_cost=octave_cost,
        octave_jump_cost=octave_jump_cost,
        voiced_unvoiced_cost=voiced_unvoiced_cost,
        pitch_ceiling=maxF)

    vr = pm.praat.call([snd,currPP,currPitch], "Voice report",0.0,0.0, minF, maxF, maximum_period_factor, maximum_amplitude_factor,silence_threshold,voicing_threshold)
    out = pd.read_table(io.StringIO(vr),delimiter=":",header=None,skip_blank_lines=True,on_bad_lines='skip').dropna()
    return(out)

#out = praat_voice_report.praat_voice_report("/Users/frkkan96/Desktop/a1.wav")

#pd.read_table(io.StringIO(out),delimiter=":\s+",header=None,skip_blank_lines=True,on_bad_lines='skip')

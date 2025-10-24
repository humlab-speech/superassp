# momelintsint.py

from pathlib import Path
import parselmouth
import math
import numpy as np
import pandas as pd
import platform
import subprocess
import io
import re
from scipy import stats
import glob
import statistics
import time
import os


# About the code
# 
# This code makes heavy use of the parselmouth library.
# For computing the prosodic measures used in the paper (doi: 10.3389/fnhum.2025.1566274)
# just call 
# 
# prosody_measures(<FULL PATH TO SOUND FILE>)
# 
# in your python code, wich the modules imported above installed in your system.
#
# To apply the function across all files in a directory and get a pandas dataframe out
# use something like the following (assumes .wav files and recursive search for sound files):
#  
# li = glob.glob("<PATH TO BASE DIRECTORY>/**/*.wav",recursive=True)
# df = pd.DataFrame(li,columns=["soundPath"],index=[os.path.basename(l).removesuffix(".wav") for l in li])
# result_df = df.apply(lambda x: prosody_measures(x['soundPath']),result_type="expand",axis=1)
# 
# Like the paper, the code is distributed under CC BY 4.0 license (https://creativecommons.org/licenses/by/4.0/).

# To perform Momel and INTSINT labeling steps, the intsint.pl (and requires a Perl installation) and a momel_* binary needs to be present.
# which is not under copyright. Please make sure that the the authors are acknowledged appropriately, e.g,:
#
#   Campione, E., Hirst, D. & Véronis, J. (n.d.). Automatic stylisation and symbolic coding of F0: implementations of the INTSINT model. Intonation. Research and Applications.
#   Hirst, D. (n.d.). A Praat plugin for Momel and INTSINT with improved algorithms for modelling and coding intonation. Proceedings of the XVIth International Conference of Phonetic Sciences, 12331236, 1223–1236.
#   Hirst, D. (2013). Melody metrics for prosodic typology: comparing English, French and Chinese. Interspeech 2013, 572–576. https://doi.org/10.21437/interspeech.2013-158
#   Hirst, D. (2019). INTSINT: a new algorithm using the OMe scale. ExLing 2018: Proceedings of 9th Tutorial and Research Workshop on Experimental Linguistics, 53–56. https://doi.org/10.36505/exling-2018/09/0012/000345
#   Hirst, D. J. (2015). Speech Prosody in Speech Synthesis: Modeling and generation of prosody for high quality and flexible speech synthesis. Prosody, Phonology and Phonetics, 3–17. https://doi.org/10.1007/978-3-662-45258-5_1
#   Hirst, D. J. & Espesser, R. (1993). Automatic Modelling Of Fundamental Frequency Using A Quadratic Spline Function. 75–85.
#
# when acknowledging the use of the code below.



# About the initial part of the code : 
# (A Python Parselmouth reimplementation of Praat script automatic_min_max_fo.praat)
#
#version:   2016-11-21
#author     Daniel Hirst
#email: daniel.hirst@lpl-aix.fr
#
#purpose: calculate fo using two passes with automatic estimation of optimal max and min
#       min fo is 0.75 * 1st quartile of fo distribution
#       max fo is min fo * max pitch span
#       values rounded down/up to nearest 10.
# Version History
# [2016-11-21]      added  variable pitch-span to control range
# [2011-02-16]      max_pitch raised to 2.5*q75 to allow for expressive speech (was 1.5*q75)
# [2008-07-11]
# 
# Reimplementation by 
#   Fredrik Nylén 

def intsint_linear (x):
    return pow(2,x)

def intsint_octave(x):
    return math.log(x)/math.log(2)

def intsint_round(x):
    return int(x + 0.5)

def getSound(file,beginTime=0.0 ,endTime=0.0 ):
    name = Path(file).stem
    if beginTime <= 0.0 and endTime <= 0.0:
        sound = parselmouth.Sound(file)
    else:
        beginTime = max(beginTime,0.0) 
        endTime = max(endTime,0.0)
        ls = parselmouth.praat.call("Open long sound file","~/Desktop/input/test.wav")
        sound = parselmouth.praat.call(ls,"Extract part",beginTime,endTime)
    parselmouth.praat.call(sound,"Rename",name)
    return sound


def automatic_min_max_fo(sound,maximum_pitch_span=1.5, minimum_fo=60, maximum_fo=750): 
    firstPass = sound.to_pitch(0.01,minimum_fo,maximum_fo)
    q25 = parselmouth.praat.call(firstPass, "Get quantile", 0, 0 , 0.25, "Hertz")
    min_fo = math.floor(q25 * 0.75/10)*10
    max_fo = math.ceil((min_fo * pow(2,maximum_pitch_span))/10)*10
    secondPass = sound.to_pitch(0.01,min_fo,max_fo)
    parselmouth.praat.call(secondPass,"Rename",str(sound.name))
    return (secondPass, min_fo, max_fo)



def momel(pitchobj,window_length=30, minimum_fo=60, maximum_fo=750, maximum_error=1.04,reduced_window_length=20, minimal_distance=5, minimal_frequency_ratio=0.05):
    momel_parameters=[window_length, minimum_fo, maximum_fo, maximum_error,reduced_window_length, minimal_distance, minimal_frequency_ratio]
    momel_parameters_str = ["%.2f" % x for x in momel_parameters]
    if platform.system() == "Darwin":
        momel_exec = "./momel_osx_intel"
    else:
        momel_exec = "./momel_linux_intel"
    pitchvalues = parselmouth.praat.call(pitchobj,"List values in all frames","Hertz")
    np.nan_to_num(pitchvalues,nan=0.0,copy=False)
    pitchvalues_str = "\n".join(pitchvalues.astype(str))
    process = subprocess.Popen([momel_exec] + momel_parameters_str, stdout=subprocess.PIPE, stdin=subprocess.PIPE,text=True)
    stdout, stderr = process.communicate(input=pitchvalues_str)
    momel_pitch_pairs = stdout.split("\n")
    momel_pitch_values = [tuple(x.split(" ")) for x in momel_pitch_pairs if x != '']
    momelPitchTier = parselmouth.praat.call("Create PitchTier",pitchobj.name,0.0,float(pitchobj.duration))
    for x in momel_pitch_values:
        # In the Praat plugin code, this is done in the perl command, but he logic is moved here since it makes more sense
        # and is more convenient to check at insert time than later.
        pitchValue = min(max(float(x[1]),float(minimum_fo)), float(maximum_fo)) 
        timeValue = float(x[0])/1000
        parselmouth.praat.call(momelPitchTier,"Add point", timeValue, pitchValue)
    
    parselmouth.praat.call(momelPitchTier,"Rename",str(pitchobj.name))
    return (momelPitchTier,momel_pitch_values)

def code_with_intsint(sound, pitchTierObj):
    epsilon_time = 0.0001
    nPoints = parselmouth.praat.call(pitchTierObj,"Get number of points")
    tout = [str(parselmouth.praat.call(pitchTierObj,"Get time from index",i) *1000) + " " + str(parselmouth.praat.call(pitchTierObj,"Get value at index",i)) for i in range(1,nPoints)]
    pitchvalues_str = "\n".join(tout)
    intsint_exec= "intsint.pl"
    perl_exec = "perl"
    process = subprocess.Popen([perl_exec, intsint_exec] , stdout=subprocess.PIPE, stdin=subprocess.PIPE,text=True)
    stdout, stderr = process.communicate(input=pitchvalues_str)

    intsint_lines = stdout.split("\n")
    intsint_header = intsint_lines[:6]

    nvalues, pmean = list(map(float, re.findall(r'\d+\.\d+|\d+', intsint_header[1]) ))
    orlow, orhigh, orstep = list(map(float, re.findall(r'\d+\.\d+|\d+', intsint_header[2]) ))
    omlow, omhigh, omstep = list(map(float, re.findall(r'\d+\.\d+|\d+', intsint_header[3]) ))
    prange = float(re.findall(r'\d+\.\d+|\d+', intsint_header[4])[0])
    key = float(re.findall(r'\d+\.\d+|\d+', intsint_header[5])[0])
    intsint_values = [ tuple( x.split(" "))  for x in intsint_lines[6:] if len(x.split(" ")) == 4]
    prev_time = 0.0
    tgMomelTier = 1
    tgIntsintTier = 2
    tgMomelIntsintTier = 3
    tg = parselmouth.praat.call(sound, "To TextGrid","Momel Intsint IntsintMomel","Momel Intsint IntsintMomel")
    parselmouth.praat.call(tg,"Rename",str(sound.name)) #May not be needed
    #nTiers = parselmouth.praat.call(tg,"Get number of tiers")
    for time, intsint, momel, intsintmomel in intsint_values:
        time = float(time)
        if abs( time- prev_time) <= epsilon_time:
            time += epsilon_time
        prev_time = time
        parselmouth.praat.call(tg,"Insert point", tgMomelTier, time, str(momel))
        parselmouth.praat.call(tg,"Insert point", tgIntsintTier, time, str(intsint))
        parselmouth.praat.call(tg,"Insert point", tgMomelIntsintTier, time, str(intsintmomel))    
    return (tg,nvalues, pmean, orlow, orhigh, orstep,omlow, omhigh, omstep ,prange, key) 

# Code below this point is not an reimplementation of the Praat plugin code

    
#Port from praatsauce to python
def correct_iseli_z (f, fx, bx, fs):
   r = math.exp(-math.pi*bx/fs)
   omega_x = 2*math.pi*fx/fs
   omega  = 2*math.pi*f/fs
   a = r ** 2 + 1 - 2*r*math.cos(omega_x + omega)
   b = r ** 2 + 1 - 2*r*math.cos(omega_x - omega)
   corr = -10*(math.log10(a)+math.log10(b))  # not normalized: H(z=0)~=0
   numerator = r ** 2 + 1 - 2 * r * math.cos(omega_x)   # sqrt of numerator, omega=0
   corr = -10*(math.log10(a)+math.log10(b)) + 20*math.log10(numerator);  # normalized: H(z=0)=1
   return corr

## bandwidth scaling factor as a function of f0, 
## to accommodate the wider bandwidths of female speech
def getbw_HawksMiller(fo, fmt):
    s = 1 + 0.25*(fo-132)/88
    if fmt < 500:
        ## coefficients for when f0<500 Hz 
        k = 165.327516
        coef =  [-6.73636734e-1, 1.80874446e-3, -4.52201682e-6, 7.49514000e-9, -4.70219241e-12 ]
    else:
        ## coefficients for when f0>=500 Hz
        k = 15.8146139
        coef = [ 8.10159009e-2, -9.79728215e-5, 5.28725064e-8, -1.07099364e-11, 7.91528509e-16 ]
    fbw = s * (k + (coef[0] * fmt) + (coef[1] * fmt**2) + (coef[2] * fmt**3) + (coef[3] * fmt**4) + (coef[4] * fmt**5) )
    return fbw

def correction_iseli_i(f, F_i, B_i, fs):
    """Return the i-th correction (dB) to the harmonic amplitude using the
       algorithm developed by Iseli and Alwan. Note this correction should be
       *subtracted* from the amplitude. The total correction is computed by
       subtracting all of the i-th corrections from the amplitude.

       Reference -- M. Iseli and A. Alwan, An improved correction formula for
       the estimation of harmonic magnitudes and its application to open
       quotient estimation.

    Args:
        f      - frequency/harmonic to be corrected (Hz) [NumPy vector]
        F_i    - i-th formant frequency (Hz) [NumPy vector]
        B_i    - i-th formant bandwidth (Hz) [NumPy vector]
        fs     - sampling frequency (Hz)
    Returns:
        corr_i - i-th correction to harmonic amplitude in dB [NumPy vector]
     
     Implementation from OpenSauce https://github.com/voicesauce/opensauce-python/blob/e633608137d4e733304321843b324db3c0b17983/opensauce/harmonics.py#L11
     
    """
    # These variable names are from the Iseli-Alwan paper
    # Normalize frequencies to sampling frequency
    rooInd = min([F_i.size, B_i.size,f.size]) #Safeguard against incompatible lengths
    corrPadLength = max([F_i.size, B_i.size,f.size]) - rooInd
    r_i = np.exp(- np.pi * B_i[:rooInd] / fs)
    omega_i = 2 * np.pi * F_i[:rooInd] / fs
    omega  = 2 * np.pi * f[:rooInd] / fs

    # Factors needed to compute correction
    numerator_sqrt = r_i**2 + 1 - 2 * r_i * np.cos(omega_i)
    denom_factor1 = r_i**2 + 1 - 2 * r_i * np.cos(omega_i + omega)
    denom_factor2 = r_i**2 + 1 - 2 * r_i * np.cos(omega_i - omega)

    # Correction in the z-domain
    # corr = 10 * log10(numerator_sqrt**2 / (denom_factor1 * denom_factor2))
    # Formula simplifies due to logarithm arithmetic
    corr_i = 20 * np.log10(numerator_sqrt) - 10 * np.log10(denom_factor1) - 10 * np.log10(denom_factor2)
    #Assume correction to be zero if too few formants have been identified 
    corr_i = np.pad(corr_i, (0, corrPadLength), 'constant', constant_values=(0,0))
    
    return corr_i

def bandwidth_hawks_miller(F_i, F0):
    """Return formant bandwidth estimated from the formant frequency and the
       fundamental frequency

       For each formant frequency, estimate the bandwidth from a 5th order
       power series with coefficients C1 or C2 depending on whether the
       frequency is less or greater than 500 Hz, then scale by a factor that
       depends on the fundamental frequency.

       Reference -- J.W. Hawks and J.D. Miller, A formant bandwidth estimation
       procedure for vowel synthesis, JASA, Vol. 97, No. 2, 1995

    Args:
        F_i - i-th formant frequency (Hz) [NumPy vector]
        F0  - Fundamental frequency (Hz) [NumPy vector]
    Returns:
        B_i - Bandwidth corresponding to i-th formant (Hz) [NumPy vector]
    Implementation from OpenSauce https://github.com/voicesauce/opensauce-python/blob/e633608137d4e733304321843b324db3c0b17983/opensauce/harmonics.py#L11
    """
    # Bandwidth scaling factor as a function of F0,
    # to accommodate the wider bandwidths of female speech
    S = 1 + 0.25 * (F0 - 132) / 88

    # Coefficients C1 (for F_i < 500 Hz) and C2 (F_i >= 500 Hz)
    C1 = np.array([165.327516, -6.73636734e-1, 1.80874446e-3, -4.52201682e-6, 7.49514000e-9, -4.70219241e-12])
    C2 = np.array([15.8146139, 8.10159009e-2, -9.79728215e-5, 5.28725064e-8, -1.07099364e-11, 7.91528509e-16])

    # Construct matrix that is a 5th order power series
    # of the formant frequency
    F_i_mat = np.vstack((F_i**0, F_i**1, F_i**2, F_i**3, F_i**4, F_i**5))

    # Construct mask for formant frequency < 500 Hz
    #
    # Set NaN values in F_i to 0, so that when we do the boolean operation
    # F_i < 500, it doesn't throw a runtime error about trying to do boolean
    # operations on NaN, which is an invalid value.
    # It doesn't matter what value we replace NaN with, because regardless of
    # the values in the Boolean mask corresponding to F_i = NaN, these will
    # get multiplied by NaN again in the formant bandwidth estimation below
    # and NaN * False = NaN * True = NaN. Any arithmetic operation on NaN
    # results in another NaN value.
    F_i_dummy = F_i.copy()
    F_i_dummy[np.isnan(F_i_dummy)] = 0
    # Tile/repeat the mask for each of the 6 terms in the power series
    mask_less_500 = np.tile(F_i_dummy < 500, (len(C1), 1))

    # Formant bandwidth estimation
    #
    # For each formant frequency, estimate the bandwidth from a 5th order power
    # series with coefficients C1 or C2 depending on whether the frequency is
    # less or greater than 500 Hz, then scale by a factor that depends on the
    # fundamental frequency
    B_i = S * (np.dot(C1, F_i_mat * mask_less_500) + np.dot(C2, F_i_mat * np.logical_not(mask_less_500)))

    return B_i




def spectral_tilt(sound,momel_pitch, formantObj, time, windowSize=60,minimum_fo=60,maximum_fo=750):
    pitch = momel_pitch
    soundDur = parselmouth.praat.call(sound,"Get total duration")
    
    windowSize = float(windowSize/1000)  # in seconds
    beginTime = max(time - (windowSize /2),0.0) # 
    endTime = min(time + (windowSize / 2),soundDur)
    soundpart = parselmouth.praat.call(sound,"Extract part",beginTime,endTime,"Gaussian1",1.0,False)
    #for C1
    #Tsiakoulis, P., Potamianos, A. & Dimitriadis, D. (2010). Spectral Moment Features Augmented by Low Order Cepstral Coefficients for Robust ASR. IEEE Signal Processing Letters, 17(6), 551–554. https://doi.org/10.1109/lsp.2010.2046349
    mfccSound = parselmouth.praat.call(soundpart,"Resample",4000,50)
    stepSize = 0.005

    #In the publication, a bandwidth of 236 Mels is reported, so the distance is between bands is 118 mels.
    mfcc = parselmouth.praat.call(mfccSound,"To MFCC",12,  0.015, stepSize, 118, 118, 0.0) #2840 
    mfccFrame = math.ceil(parselmouth.praat.call(mfcc,"Get number of frames")/2)
    C1 = parselmouth.praat.call(mfcc,"Get value in frame",mfccFrame,1)

    spct = parselmouth.praat.call(soundpart,"To Spectrum","yes")
    ltas = parselmouth.praat.call(spct,"To Ltas (1-to-1)")

    Fi = [parselmouth.praat.call(formantObj,"Get value at time",i+1,time, "Hertz","Linear") for i in range(parselmouth.praat.call(formantObj,"Get minimum number of formants"))]

    Bi = bandwidth_hawks_miller(np.array(Fi),pitch)

    #Spectral balance
    #Sluijter, A. M. C. & Heuven, V. J. van. (1996). Spectral balance as an acoustic correlate of linguistic stress. The Journal of the Acoustical Society of America, 100(4), 2471–2485. https://doi.org/10.1121/1.417955
    spectralbalance = parselmouth.praat.call(ltas,"Get slope",0,500,500,1000,"energy")
    #SLF / Spectral tilt
    #Spectral tilt : Schweitzer, A. (2019). Exemplar-theoretic integration of phonetics and phonology: Detecting prominence categories in phonetic space. Journal of Phonetics, 77, 100915. https://doi.org/10.1016/j.wocn.2019.100915
    slR = parselmouth.praat.call(ltas,"Report spectral trend",100,5000, "logarithmic","least squares")

     # Kakouros, S., Räsänen, O. & Alku, P. (2018). Comparison of spectral tilt measures for sentence prominence in speech—Effects of dimensionality and adverse noise conditions. Speech Communication, 103, 11–26. https://doi.org/10.1016/j.specom.2018.08.002
    SLF = float(slR.split()[11])
    soundpart.to_spectrogram()
    spectrogram= sound.to_spectrogram()
    log_magnitude_spectrum = np.log1p(spectrogram.values)
    mean_log_magnitude = np.mean(log_magnitude_spectrum, axis=1)
    SLF6D_coefficients = np.polyfit(np.arange(len(mean_log_magnitude)), mean_log_magnitude, 6)
    SLF6D_dict = zip(["SLF6D." + str(i) for i in range(1,6)], SLF6D_coefficients)
    
    lowerbh = [ (pitch * (n+1)) - (pitch /10) for n in range(4)]
    upperbh = [ (pitch * (n+1)) + (pitch /10) for n in range(4)]

    Ln = [parselmouth.praat.call(ltas,"Get maximum", lowerbh[i],upperbh[i],"Parabolic") for i in range(4)]
    L_Fn = [parselmouth.praat.call(ltas,"Get value at frequency", Fi[i],"Linear") for i in range(len(Fi))]
    nfo = [parselmouth.praat.call(ltas,"Get frequency of maximum", lowerbh[i],upperbh[i],"Parabolic") for i in range(4)]
    corr = correction_iseli_i(np.array(nfo[0:2] + nfo[3:]), np.array(Fi[0:3]), Bi[0:3], sound.sampling_frequency)
    
    Ln_c = Ln
    for i in range(len(corr)):
        Ln_c = Ln_c - corr[i] #Apply correction
    #Formant aplitude correction
    corr = correction_iseli_i(np.array(L_Fn), np.array(Fi), Bi, sound.sampling_frequency)
    L_Fn_c = L_Fn

    for i in range(len(corr)):
        L_Fn_c = L_Fn_c - corr[i]
    #H1-H2 (dB)
    #Campbell, N. & Beckman, M. (1997). Stress, prominence, and spectral tilt. In A. Botinis, G. Kouroupetroglou & G. Carayiannis (Eds.), Intonation: Theory, models and applications (pp. 67–70). ESCA and University of Athens Department of Informatics.
    L2L1 = Ln[0] - Ln[1]
    L2cL1c = Ln_c[0] - Ln_c[1]
    #H1*-A3* (dB), which is Ln[0] - 
    #Okobi, A. O. (2006). Acoustic correlates of word stress in American English. Massachusetts Institute of Technology.
    L1cLF3c = float(Ln_c[0] - L_Fn_c[2]) if len(L_Fn_c) > 2 else None
    L1LF3 = float(Ln[0] - L_Fn[2]) if len(L_Fn) > 2 else None

    out = {'L2L1': float(L2L1),'L2cL1c':float(L2cL1c),'L1cLF3c': L1cLF3c,'L1LF3': L1LF3,'SLF': float(SLF),'C1':float(C1),'Spectral Balance':float(spectralbalance)}
    out.update(SLF6D_dict)
    return out


def prosody_measures(soundPath,minF=60,maxF=750,windowShift=10):
    windowShift = float(windowShift/1000)  # in seconds
    soundObj = getSound(soundPath)
    print("Processing: " + soundPath)
    duration = float(parselmouth.praat.call(soundObj,"Get total duration"))
    if duration < 1.0:
        print("Skopping utterance with less than 1 second duration: " + soundPath)
        return None # Skipp sentences that have less than 1 second duration
    pitchObj, min_fo, max_fo = automatic_min_max_fo(soundObj,maximum_pitch_span=1.5, minimum_fo=minF, maximum_fo=maxF)
    momelPitchTier,momel_pitch_values = momel(pitchObj,minimum_fo=min_fo, maximum_fo=max_fo)
    textGrid,nintsint, pitchmean, optim_r_low, optim_r_high, optim_r_step, optim_m_low, optim_m_high, optim_m_step, pitchrange, key = code_with_intsint(soundObj,momelPitchTier)
    formantObj  = parselmouth.praat.call(soundObj,"To Formant (burg)",0.001, 5, 5000, 0.025, 50)
    intensityObj = parselmouth.praat.call(soundObj,"To Intensity",min_fo, windowShift,"yes")
    tgTable = pd.read_table(io.StringIO(parselmouth.praat.call(textGrid, "List", False,20,True,False))) #Pandas version
    tgTabWide = pd.pivot(tgTable[['tmin','tier','text']],index="tmin",columns="tier",values="text")
    
    tgTabWide['Intensity'] = tgTabWide.apply(lambda x: parselmouth.praat.call(intensityObj,"Get value at time",x.name,"cubic"),axis=1)
     #Merge in computed momel_pitch values
    momelPitch = pd.DataFrame([(float(a)/1000.0, float(b)) for a, b in momel_pitch_values],columns=["time","momel_pitch"],index=pd.Index([float(a)/1000 for a, b in momel_pitch_values],name="tmin"))
    tgTabWide = pd.merge_asof(tgTabWide,momelPitch,left_index=True, right_index=True)
    tgTabWide.dropna(axis=0,inplace=True) #Drop initial or end markers which may not have times associated with them

    #Merge in spectral balance measures
    specTilt = tgTabWide.apply(lambda x: spectral_tilt(soundObj, float(x['momel_pitch']),formantObj,x.name), axis=1,result_type="expand")
    tgTabWide = tgTabWide.join(specTilt)
    
    #Compute Inter-intsintlabel differences for all associated spectral properties across the utterance
    tgTabWide_toDiff = tgTabWide[tgTabWide.columns[3:18]]

    # Apply differentiation
    tgTabDiff = pd.DataFrame(np.diff(tgTabWide_toDiff, axis=0), columns=tgTabWide_toDiff.columns)
    tgTabDiff = pd.concat([pd.DataFrame(np.zeros((1, tgTabWide_toDiff.shape[1])), columns=tgTabWide_toDiff.columns), tgTabDiff], ignore_index=True)

    # Display the differentiated DataFrame
    tgTabDiff.columns = tgTabDiff.columns + "_diff"
    tgTabDiff.index = tgTabWide_toDiff.index
    
    tgWideTab =  tgTabDiff.join(tgTabWide_toDiff)

    # Define a function to compute statistics safely
    def safe_statistics(series):
        results = {}
        try:
            results['tstd'] = stats.tstd(series)
            results['tmean'] = stats.tmean(series)
            results['variation'] = stats.variation(series)
            results['iqr'] = stats.iqr(series)
            results['tmax'] = stats.tmax(series)
            results['tmin'] = stats.tmin(series)
        except ValueError:
            # If a ValueError occurs, return NaN for all statistics
            results = {'tstd': None, 'tmean': None, 'variation': None, 'iqr': None, 'tmax': None, 'tmin': None}
        return pd.Series(results)

    # Safe application of stats functions across all columns 
    tgTabSummary = tgWideTab.apply(lambda x: safe_statistics(x), axis=0, result_type='expand')
    tgTabSummary.columns = ['_'.join(col).strip() for col in tgTabSummary.columns.values]

    tgTabSummary.reset_index(names="statistic",inplace=True)
    tgTabSummary_long =tgTabSummary.melt(id_vars="statistic")
    tgTabSummary_long["variable_name"] =  tgTabSummary_long.variable + "_" + tgTabSummary_long.statistic
 
    #Fix momel_pitch problem and rename incidies
    tgTabSummary.rename(index={"tstd":"std","tmean":"mean","hdmedian":"median", "variation":"var","iqr":"iqr","tmax":"max","tmin":"min"},columns={"momel_pitch":"momelpitch"},inplace=True)
    # flat indicies works better in multi-language contexts and when merging
    tgTabSummary.columns = tgTabSummary.columns.to_flat_index() + "_diff"

    #Here we start preparing the results dictionary
    tgSummary = dict(zip(tgTabSummary_long.variable_name, tgTabSummary_long.value))
    nUniqueIntsint = int(tgTabWide[["Intsint"]].nunique().iloc[0])
    duration = float(parselmouth.praat.call(soundObj,"Get total duration"))
    additionalInfo = {"Duration": duration ,"PitchKey":key,"PitchRange":pitchrange,"PitchMean":pitchmean,"IntsIntLabels":nintsint,"UniqueIntsInt":nUniqueIntsint,"IntsIntConcentration":float(nintsint/duration),"OptimizationRangeLow":optim_r_low,"OptimizationRangeHigh":optim_r_high,"OptimizationStep":optim_r_step, "OptimizationMidLow":optim_m_low, "OptimizationMidHigh": optim_m_high, "OptimizationMidStep":optim_m_step}

    return pd.Series(tgSummary)


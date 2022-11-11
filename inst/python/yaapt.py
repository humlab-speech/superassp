soundFile = "/Users/frkkan96/Desktop/a1.wav"
minF = 70.0
maxF =400.0
windowShift = 10.0
windowSize = 35
beginTime = 1.2
endTime = 1.6
tda_frame_length = 35
fft_length = 8192
bp_forder=  150
bp_low= 50 
bp_high=  1500 
nlfer_thresh1= 0.75
nlfer_thresh2=  0.1
shc_numharms= 3
shc_window= 40 
shc_maxpeaks=  4
shc_pwidth=  50 
shc_thresh1=  5
shc_thresh2=  1.25
f0_double= 150 
f0_half= 150 
dp5_k1= 11
dec_factor=  1
nccf_thresh1=  0.3
nccf_thresh2= 0.9
nccf_maxcands=  3
nccf_pwidth=  5
merit_boost=  0.20
merit_pivot=  0.99
merit_extra=  0.4
median_value=  7
dp_w1=  0.15
dp_w2=  0.5
dp_w3= 0.1
dp_w4= 0.9

import amfm_decompy.pYAAPT as pYAAPT
import amfm_decompy.basic_tools as basic
import math as m

signal = basic.SignalObj(soundFile)
if endTime > 0.0 or beginTime > 0.0:
	startSample = m.floor(beginTime * signal.fs)
	endSample = m.ceil(endTime * signal.fs)
	subsignal = basic.SignalObj(signal.data[startSample:endSample],signal.fs)
else:
	subsignal = signal

pitch = pYAAPT.yaapt(subsignal, **{'f0_min' : minF, 
	'tda_frame_length' : tda_frame_length,
	'f0_max' : maxF, 
	'frame_length' : windowSize, 
	'frame_space' : windowShift,
	'tda_frame_length' : tda_frame_length,
	'fft_length' : fft_length,
	'bp_forder': bp_forder,
	'bp_low' : bp_low,
	'bp_high' : bp_high,
	'nlfer_thresh1' : nlfer_thresh1,
	'nlfer_thresh2' : nlfer_thresh2,
	'shc_numharms' : shc_numharms,
	'shc_window' : shc_window,
	'shc_maxpeaks' : shc_maxpeaks,
	'shc_pwidth' : shc_pwidth,
	'shc_thresh1' : shc_thresh1,
	'shc_thresh2' : shc_thresh2,
	'f0_double' : f0_double,
	'f0_half' : f0_half,
	'dp5_k1' : dp5_k1,
	'dec_factor' : dec_factor,
	'nccf_thresh1' : nccf_thresh1,
	'nccf_thresh2' : nccf_thresh2,
	'nccf_maxcands' : nccf_maxcands,
	'nccf_pwidth' : nccf_pwidth,
	'merit_boost' : merit_boost,
	'merit_pivot' : merit_pivot,
	'merit_extra'  : merit_extra,
	'median_value' : median_value,
	'dp_w1' : dp_w1,
	'dp_w2' : dp_w2,
	'dp_w3' : dp_w3,
	'dp_w4' : dp_w4 })

f0 =  pitch.samp_values
vuv =   pitch.vuv



py$endTime <- reticulate::r_to_py(endTime)

py$tda_frame_length <- reticulate::r_to_py(tda_frame_length)
py$fft_length <- reticulate::r_to_py(fft_length )
py$bp_forder <- reticulate::r_to_py(bp_forder)
py$bp_low <- reticulate::r_to_py(bp_low)
py$bp_high <- reticulate::r_to_py(bp_high)
py$nlfer_thresh1 <- reticulate::r_to_py(nlfer_thresh1)
py$nlfer_thresh2 <- reticulate::r_to_py(nlfer_thresh2)
py$shc_numharms <- reticulate::r_to_py(shc_numharms)
py$shc_window <- reticulate::r_to_py(shc_window)
py$shc_maxpeaks <- reticulate::r_to_py(shc_maxpeaks)
py$shc_pwidth <- reticulate::r_to_py(shc_pwidth)
py$shc_thresh1 <- reticulate::r_to_py(shc_thresh1)
py$shc_thresh2 <- reticulate::r_to_py(shc_thresh2)
py$f0_double <- reticulate::r_to_py(f0_double)
py$f0_half <- reticulate::r_to_py(f0_half)
py$dp5_k1 <- reticulate::r_to_py(dp5_k1)
py$dec_factor <- reticulate::r_to_py(dec_factor)
py$nccf_thresh1 <- reticulate::r_to_py(nccf_thresh1)
py$nccf_thresh2 <- reticulate::r_to_py(nccf_thresh2)
py$nccf_maxcands <- reticulate::r_to_py(nccf_maxcands)
py$nccf_pwidth <- reticulate::r_to_py(nccf_pwidth)
py$merit_boost <- reticulate::r_to_py(merit_boost)
py$merit_pivot <- reticulate::r_to_py(merit_pivot)
py$merit_extra <- reticulate::r_to_py(merit_extra)
py$median_value <- reticulate::r_to_py(median_value)
py$dp_w1 <- reticulate::r_to_py(dp_w1)
py$dp_w2 <- reticulate::r_to_py(dp_w2)
py$dp_w3 <- reticulate::r_to_py(dp_w3)
py$dp_w4 <- reticulate::r_to_py(dp_w4)
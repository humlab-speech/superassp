soundFile = "/Users/frkkan96/Desktop/a1.wav"
soundFile = "Data/VISP_emuDB/Svenska_ses/svenska_bndl/svenska.wav"
minF = 70.0
maxF =400.0
voiced_voiceless_threshold = 0.1
windowShift = 10.0
beginTime = 1.2
endTime = 1.6
targetSampleRate=16000

duration = endTime - beginTime

if duration < (windowShift / 1000) :
	duration = None


import tensorflow as tf
import tensorflow_hub as hub

import numpy as np
import librosa as lr


audio_samples, fs = lr.load(soundFile,
	dtype= tf.float32.as_numpy_dtype,
	sr=targetSampleRate,
	offset= beginTime,
	duration= duration
	)



# We now feed the audio to the SPICE tf.hub model to obtain pitch and uncertainty outputs as tensors.
model_output = model.signatures["serving_default"](tf.constant(audio_samples, tf.float32))

pitch_outputs = model_output["pitch"]
uncertainty_outputs = model_output["uncertainty"]


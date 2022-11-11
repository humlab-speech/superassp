soundFile = "/Users/frkkan96/Desktop/a1.wav"
beginTime = 1.2
endTime = 1.6


if endTime == 0.0
	endTime = None


import opensmile as o
import numpy as np


smile = opensmile.Smile(
    feature_set=opensmile.FeatureSet.GeMAPS,
    feature_level=opensmile.FeatureLevel.Functionals,
)

smile_results = smile.process_file(file=soundFile,
	start=beginTime,
	end=endTime)
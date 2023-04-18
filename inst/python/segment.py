#!/Users/frkkan96/opt/anaconda3/envs/pyannote/bin/python

import sys
import pandas as pd
import os

wavfile = sys.argv[1]
token = sys.argv[2]

# instantiate pretrained speaker diarization pipeline
#from pyannote.audio import Pipeline
#pipeline = Pipeline.from_pretrained("pyannote/segmentation")

from pyannote.audio import Model
model = Model.from_pretrained("pyannote/segmentation", 
                              use_auth_token=token)


from pyannote.audio.pipelines import VoiceActivityDetection
pipeline = VoiceActivityDetection(segmentation=model)
HYPER_PARAMETERS = {
  # onset/offset activation thresholds
  "onset": 0.5, "offset": 0.5,
  # remove speech regions shorter than that many seconds.
  "min_duration_on": 0.0,
  # fill non-speech regions shorter than that many seconds.
  "min_duration_off": 0.0
}
pipeline.instantiate(HYPER_PARAMETERS)

# apply pretrained pipeline
vad = pipeline(wavfile)

arr = []

# print the result
for turn, _, speaker in vad.itertracks(yield_label=True):
    # print(f"start={turn.start:.1f}s stop={turn.end:.1f}s speaker_{speaker}")
    arr.append( {"start":turn.start, "end": turn.end, "speaker":speaker, "file":wavfile} )

outpd = pd.DataFrame(arr)


csvFile = os.path.splitext(wavfile)[0]+'_segment.csv'


outpd.to_csv(csvFile,sep=",", columns=['start','end','speaker','file'],encoding='utf-8',index=False)
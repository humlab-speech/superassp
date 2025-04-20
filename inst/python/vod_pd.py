#!/Users/frkkan96/opt/anaconda3/envs/pyannote/bin/python


wavfile = sys.argv[1]
YOUR_AUTH_TOKEN = sys.argv[2]

import sys
import pandas as pd

from pyannote.audio import Model
from pyannote.audio.pipelines import VoiceActivityDetection

def voice_activity_detection(wavfile, session, level, attribute, YOUR_AUTH_TOKEN):
    # instantiate the model
    model = Model.from_pretrained(
      "pyannote/segmentation-3.0", 
      use_auth_token=YOUR_AUTH_TOKEN)
    
    pipeline = VoiceActivityDetection(segmentation=model)
    HYPER_PARAMETERS = {
      # remove speech regions shorter than that many seconds.
      "min_duration_on": 0.0,
      # fill non-speech regions shorter than that many seconds.
      "min_duration_off": 0.0
    }
    pipeline.instantiate(HYPER_PARAMETERS)
    vad = pipeline(wavfile)
    
    vaddf = pd.DataFrame.from_records(vad.itersegments(),columns=["start","end"])
    vaddf["bundle"] = wavfile.removesuffix(".wav")
    vaddf["session"] = session
    vaddf["level"] = level
    vaddf["attribute"] = attribute
    vaddf["labels"] = "SPEECH"
    
    vaddf['start_item_seq_idx'] = vaddf.index
    return vaddf


conda create -n pysuperassp python=3.8
conda activate pysuperassp

conda install pytorch torchaudio -c pytorch

pip install https://github.com/pyannote/pyannote-audio/archive/develop.zip
pip install speechbrain
pip install opensmile
pip install pyreaper
pip install pysptk
pip install librosa
# 
# 
# # R code
# library(reticulate)
reticulate::conda_create("pysuperassp")
reticulate::py_install("pytorch",channel="pytorch")
reticulate::py_install("torchaudio",channel="pytorch")
# 
reticulate::py_install("https://github.com/pyannote/pyannote-audio/archive/develop.zip",pip=TRUE)
reticulate::py_install("speechbrain",pip=TRUE)
reticulate::py_install("opensmile",pip=TRUE)
reticulate::py_install("pyreaper",pip=TRUE)
reticulate::py_install("pysptk",pip=TRUE)
reticulate::py_install("librosa",pip=FALSE)
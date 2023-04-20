conda create -n pysuperassp

conda install pytorch torchaudio -c pytorch
conda install torchcrepe

pip install https://github.com/pyannote/pyannote-audio/archive/develop.zip
pip install speechbrain
pip install opensmile
pip install pyreaper
pip install pysptk
pip install librosa
pip install pyworld
pip install amfm_decompy
pip install tensorflow
pip install tensorflow_hub

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
reticulate::py_install("pyworld",pip=TRUE)
reticulate::py_install("amfm_decompy",pip=TRUE)
reticulate::py_install("tensorflow",pip=TRUE)



conda create conda create --prefix -n pysuperassp python=3.8 
conda activate pysuperassp

#conda create --prefix ./envs python=3.8 numpy pandas
#conda activate ./envs
conda install pytorch torchaudio -c pytorch

pip install https://github.com/pyannote/pyannote-audio/archive/develop.zip

pip install speechbrain

pip install opensmile

pip install pyreaper
#git clone https://github.com/MLSpeech/Dr.VOT.git
#cd Dr.VOT 
#./check_installations.sh 

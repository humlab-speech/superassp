Flytta tar-fil till där det ska köras:

 

sudo docker load < voice_analysis_docker.tar

 

 

Utskriften gjordes till /mnt/temp för att där mountar jag med kör-kommandot:

docker run --rm -e "DISPLAY=:0" -v /tmp/.11-unix:/tmp/.X11-unix -v /mnt/c/Users/tosabn02/Documents/MATLAB/VoiceAnalysisToolbox/VoiceAnalysisToolbox/:/mnt/temp/ voice_analysis "/mnt/temp/aaa.wav"

 

Eller minimalt:

 

 docker run --rm -v /mnt/c/Users/tosabn02/Documents/MATLAB/:/mnt/matlab_directory/ voice_analysis "/mnt/matlab_directory/aaa.wav"

 

 


	
/mnt/temp - volymen innehåller wav-filen som ska köras
	
Argumentet I slutet blir då alltid /mnt/temp/filnamn
	
Output kommer skrivas I den mountade katalogen som filnamn_voice_analysis_docker.mat
	

		
Notera att scriptet skriver alltid I /mnt/matlab_directory/ även om annan input path ges; då kommer ej resultatet ut från containern!
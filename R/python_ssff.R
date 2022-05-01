

swipe <- function(listOfFiles,
                  beginTime=0,
                  endTime=0,
                  windowShift=5, 
                  minF=70, 
                  maxF=200, 
                  voicing.threshold=0.3,
                  explicitExt="swi",
                  outputDirectory=NULL,
                  toFile=TRUE, 
                  conda.env="pysuperassp"){
  
  if(!conda.env %in%  reticulate::conda_list()$name){
    stop("The conda environment ",conda.env, " does not exist.\n Please ")
  }
  
  if(length(listOfFiles) > 1 & ! toFile){
    stop("length(listOfFiles) is > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.")
  }
  
  tryCatch({
    fileBeginEnd <- data.frame(
      listOfFiles = listOfFiles, 
      beginTime = beginTime,
      endTime=endTime
    )
  },error=function(e){stop("The beginTime and endTime must either be a single value or the same length as listOfFiles")})
  #Check that all files exists before we begin
  filesEx <- file.exists(listOfFiles)
  if(!all(filesEx)){
    filedNotExists <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ",paste(filedNotExists, collapse = ", "))
  }
  #The empty vector of file names that should be returned
  outListOfFiles <- c()
  
  for(i in 1:nrow(fileBeginEnd)){ 
    origSoundFile <- fileBeginEnd[i, "listOfFiles"]
    
    beginTime <- fileBeginEnd[i, "beginTime"]
    endTime <- fileBeginEnd[i, "endTime"]
        

    reticulate::use_condaenv(conda.env)
    
    py$soundFile <- r_to_py(origSoundFile)
    py$ws <- r_to_py(windowShift)
    py$fMax <- r_to_py(maxF)
    py$fMin <- r_to_py(minF)
    py$bt <- r_to_py(beginTime)
    py$et <- r_to_py(beginTime)
    
    reticulate::py_run_string("import numpy as np\
import pysptk\
from scipy.io import wavfile\
fs, x = wavfile.read(soundFile)\
begin = int(bt * fs)\
end = int(x.size-1) if et == 0 else min(int(et * fs), x.size-1 ) \

w = x[range(begin,end)]\
f0_swipe = pysptk.swipe(w.astype(np.float64), fs=fs, hopsize=ws / 1000 * fs, min=fMin, max=fMax, otype=\"f0\")\
pitch_swipe = pysptk.swipe(w.astype(np.float64), fs=fs, hopsize=ws / 1000 * fs, min=fMin, max=fMax, otype=\"pitch\")")
    
    inTable <- data.frame( "f0" = py$f0_swipe,
                           "pitch"=py$pitch_swipe)
  
    startTime = windowShift
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16", "INT16")
    #Use the time separation between second and pitch measurement time stamps to compute a sample frequency.
    
    sampleRate <-  1/ windowShift * 1000
    attr(outDataObj, "sampleRate") <- sampleRate
    
    attr(outDataObj, "origFreq") <-  as.numeric(py$fs) 
    startTime <- 1/sampleRate
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
    class(outDataObj) = "AsspDataObj"
    
    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    # Cross-correlation track
    f0Table <- inTable %>%
      dplyr::select(f0) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    nof0Values <- nrow(f0Table)
    names(f0Table) <- NULL
    outDataObj = wrassp::addTrack(outDataObj, "f0", as.matrix(f0Table[,1]), "INT16")
    
    # Auto-correlation track
    pitchTable <- inTable %>%
      dplyr::select(pitch) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    noPitchValues <- nrow(pitchTable)
    names(pitchTable) <- NULL
    outDataObj = wrassp::addTrack(outDataObj, "pitch", as.matrix(pitchTable[,1]), "INT16")
    
    
    
    #return(outDataObj)
    ## Apply fix from Emu-SDMS manual
    ##https://raw.githubusercontent.com/IPS-LMU/The-EMU-SDMS-Manual/master/R/praatToFormants2AsspDataObj.R
    
    # add missing values at the start as Praat sometimes
    # has very late start values which causes issues
    # in the SSFF file format as this sets the startRecord
    # depending on the start time of the first sample
    if( startTime > (1/sampleRate) ){
      
      nr_of_missing_samples = as.integer(floor(startTime / (1/sampleRate)))
      
      missing_f0_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$f0))
      missing_pitch_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$pitch))
      
      # prepend values
      outDataObj$f0 = rbind(missing_f0_vals, outDataObj$f0)
      outDataObj$pitch = rbind(missing_pitch_vals, outDataObj$pitch)
      
      
      # fix start time
      attr(outDataObj, "startTime") = startTime - nr_of_missing_samples * (1/sampleRate)
    }
    
    assertthat::assert_that(wrassp::is.AsspDataObj(outDataObj),
                            msg = "The AsspDataObj created by the swipe function is invalid.")
    
    ssff_file <- sub("wav$",explicitExt,origSoundFile)
    if(!is.null(outputDirectory)){
      ssff_file <- file.path(outputDirectory,basename(ssff_file))
    }
    
    attr(outDataObj,"filePath") <- as.character(ssff_file)
    if(toFile){
      wrassp::write.AsspDataObj(dobj=outDataObj,file=ssff_file)
      #Here we can be sure that the list is a valid SSFF object, so the
      # so we add TRUE to the out vector
      outListOfFiles <- c(listOfFiles,TRUE)
    }
  }
  if(toFile){
    return(length(outListOfFiles))
  }else{
    return(outDataObj)
  }
  
}



attr(swipe,"ext") <-  c("swi") 
attr(swipe,"tracks") <-  c("f0","pitch")
attr(swipe,"outputType") <-  c("SSFF")


## FOR INTERACTIVE TESTING
library(superassp)
library(reticulate)
library(dplyr)
f <- "/Users/frkkan96/Desktop/a1.wav"
swipe(f,beginTime=0,endTime=1,toFile=FALSE) -> outportion
swipe(f,beginTime=0,endTime=0,toFile=FALSE) -> out
praat_pitch(f,toFile=FALSE,corr.only=TRUE) -> out2 
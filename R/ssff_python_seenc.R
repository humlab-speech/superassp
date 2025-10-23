#' (En)coded spectral envelope estimation
#'
#' The Spectral envelope is estimated using function for assessing band-aperiodicities
#' using the CheapTrick algorithm
#' \insertCite{Morise:2015ia}{superassp} implemented
#' in the WORLD vocoder
#' \insertCite{MORISE.2016.10.1587/transinf.2015edp7457}{superassp}. The the `harvest()`
#' pitch algorithm \insertCite{Morise.2017.10.21437/interspeech.2017-68}{superassp} is used to calculate
#' the periodic component.
#'
#'
#' @inheritParams trk_harvest
#' @param dimensions Number of dimensions of coded spectral envelope
#'
#' @return An SSFF track object containing two tracks (f0 and corr) that are
#'   either returned (toFile == FALSE) or stored on disk.
#' @references \insertAllCited{}
#' 
#' @export
#'
#'   
trk_seenc <- function(listOfFiles,
                 beginTime=0,
                 endTime=0,
                 windowShift=5,
                 minF=70, 
                 maxF=200, 
                 dimensions = 1,
                 explicitExt="sec",
                 outputDirectory=NULL,
                 toFile=TRUE){
  
  
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
    origSoundFile <- normalizePath(fileBeginEnd[i, "listOfFiles"],mustWork = TRUE)
    
    beginTime <- fileBeginEnd[i, "beginTime"]
    endTime <- fileBeginEnd[i, "endTime"]
    
    py$soundFile <- reticulate::r_to_py(origSoundFile)
    py$windowShift <- reticulate::r_to_py(windowShift)
    py$maxF <- reticulate::r_to_py(as.double(maxF))
    py$minF <- reticulate::r_to_py(as.double(minF))
    py$beginTime <- reticulate::r_to_py(as.double(beginTime))
    py$endTime <- reticulate::r_to_py(as.double(endTime))
    py$dimensions <- reticulate::r_to_py(as.integer(dimensions))
    
    reticulate::py_run_string("duration = endTime - beginTime \
import gc\
import pyworld as pw\
import librosa as lr\
import numpy as np\
\
if duration < (windowShift / 1000) :\
	duration = None\
\
x, fs = lr.load(soundFile,\
	dtype=np.float64,\
	offset= beginTime,\
	duration= duration\
	)\

f0, t = pw.trk_harvest(x,\
	fs,\
	f0_floor=minF,\
	f0_ceil=maxF,\
	frame_period=windowShift )\
\
sp = pw.cheaptrick(x,f0, t,fs, f0_floor=minF)\
\
sl = pw.code_spectral_envelope(sp,fs,dimensions)\
del x\
del f0\
del sp\
gc.collect()")
    
    seencTable <- as.data.frame( py$sl) %>% 
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.double))
    
    
    startTime = as.numeric(py$t)[[1]]
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("REAL32")
    #Use the time separation between second and pitch measurement time stamps to compute a sample frequency.
    
    sampleRate <-  1/ windowShift * 1000
    attr(outDataObj, "sampleRate") <- sampleRate
    
    attr(outDataObj, "origFreq") <-  as.numeric(py$fs) 
    #startTime <- 1/sampleRate
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(seencTable))
    class(outDataObj) = "AsspDataObj"
    
    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    
    
    noSeencValues <- nrow(seencTable)
    names(seencTable) <- NULL
    
    outDataObj = addTrack(outDataObj, "trk_seenc", as.matrix(seencTable), "REAL32")
    
    
    ## Apply fix from Emu-SDMS manual
    ##https://raw.githubusercontent.com/IPS-LMU/The-EMU-SDMS-Manual/master/R/praatToFormants2AsspDataObj.R
    
    # add missing values at the start as Praat sometimes
    # has very late start values which causes issues
    # in the SSFF file format as this sets the startRecord
    # depending on the start time of the first sample
    if( startTime > (1/sampleRate) ){
      
      nr_of_missing_samples = as.integer(floor(startTime / (1/sampleRate)))
      
      missing_seenc_vals = matrix(0,
                                         nrow = nr_of_missing_samples,
                                         ncol = ncol(outDataObj$seenc))
      
      
      # prepend values
      outDataObj$seenc = rbind(missing_seenc_vals, outDataObj$seenc)   
      
      # fix start time
      attr(outDataObj, "startTime") = startTime - nr_of_missing_samples * (1/sampleRate)
    }
    
    assertthat::assert_that(is.AsspDataObj(outDataObj),
                            msg = "The AsspDataObj created by the swipe function is invalid.")
    
    ssff_file <- sub("wav$",explicitExt,origSoundFile)
    if(!is.null(outputDirectory)){
      ssff_file <- file.path(outputDirectory,basename(ssff_file))
    }
    
    attr(outDataObj,"filePath") <- as.character(ssff_file)
    if(toFile){
      write.AsspDataObj(dobj=outDataObj,file=ssff_file)
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



attr(trk_seenc,"ext") <-  c("sec") 
attr(trk_seenc,"tracks") <-  c("trk_seenc")
attr(trk_seenc,"outputType") <-  c("SSFF")

## FOR INTERACTIVE TESTING
#library(superassp)
#library(reticulate)
#library(dplyr)
#f <- "/Users/frkkan96/Desktop/a1.wav"
#trk_swipe(f,beginTime=0,endTime=1,toFile=FALSE) -> outportion
#trk_swipe(f,beginTime=0,endTime=0,toFile=FALSE) -> outsa
#trk_reaper(f,beginTime=0,endTime=1,toFile=FALSE) -> outportion
#trk_excite(f,beginTime=0,endTime=0,toFile=FALSE) -> outex
#reaper_pm(f,beginTime=0,endTime=0,toFile=FALSE) -> outpm
#trk_kaldi_pitch(f,beginTime=0,endTime=0,toFile=FALSE) -> outkaldi
#trk_kaldi_pitch(f,beginTime=0,endTime=0,toFile=TRUE)
#trk_dio(f,beginTime=0,endTime=0,toFile=FALSE) -> out
#trk_dio(f,beginTime=0,endTime=0,toFile=TRUE)
#trk_kaldi_pitch(f,beginTime=0,endTime=0,toFile=FALSE) -> out
#trk_harvest(f,beginTime=0,endTime=0,toFile=FALSE) -> out1
#trk_seenc(f,beginTime=0,endTime=0,toFile=FALSE) -> out2

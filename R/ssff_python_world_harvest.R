#' 
#' The Harvest algorithm \insertCite{Morise.2017.10.21437/interspeech.2017-68}{superassp} was developed for the WORLD vocoder \insertCite{MORISE.2016.10.1587/transinf.2015edp7457}{superassp} and aims to obtain a reliable F0 contour and reduce erroneously identified voice frames. 
#' 
#' The algorithm consists of two steps. In the first step, the algorithm uses fundamental component extraction by many band-pass filters with different center frequencies and obtains the basic f0 candidates from filtered signals. In the second step, basic f0 candidates are refined and scored by using the instantaneous frequency, and then several f0 candidates in each frame are estimated.
#'  
#' @inheritParams trk_swipe
#'
#' @return An SSFF track object containing two tracks (f0 and corr) that are
#'   either returned (toFile == FALSE) or stored on disk.
#' @references 
#'   \insertAllCited{}
#' @keywords internal
#' @noRd
#'
harvest_python<- function(listOfFiles,
               beginTime=0,
               endTime=0,
               windowShift=5,
               minF=70, 
               maxF=200, 
               explicitExt="hf0",
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
    py$maxF <- reticulate::r_to_py(maxF)
    py$minF <- reticulate::r_to_py(minF)
    py$beginTime <- reticulate::r_to_py(beginTime)
    py$endTime <- reticulate::r_to_py(endTime)
    
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
\
f0, t = pw.trk_harvest(x,\
	fs,\
	f0_floor=minF,\
	f0_ceil=maxF,\
	frame_period=windowShift)\
del x\
gc.collect()")
    
    inTable <- data.frame( "f0" = py$f0)
    
    
    startTime = as.numeric(py$t)[[1]]
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16")
    #Use the time separation between second and pitch measurement time stamps to compute a sample frequency.
    
    sampleRate <-  1/ windowShift * 1000
    attr(outDataObj, "sampleRate") <- sampleRate
    
    attr(outDataObj, "origFreq") <-  as.numeric(py$fs) 
    #startTime <- 1/sampleRate
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
    class(outDataObj) = "AsspDataObj"
    
    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    # Pitch track
    f0Table <- inTable %>%
      dplyr::select(f0) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    nof0Values <- nrow(f0Table)
    names(f0Table) <- NULL
    
    outDataObj = addTrack(outDataObj, "f0", as.matrix(f0Table[,1]), "INT16")
    
    
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
      
      
      # prepend values
      outDataObj$f0 = rbind(missing_f0_vals, outDataObj$f0)   
      
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

attr(trk_harvest,"ext") <-  c("hf0") 
attr(trk_harvest,"tracks") <-  c("f0")
attr(trk_harvest,"outputType") <-  c("SSFF")



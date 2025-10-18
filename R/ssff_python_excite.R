excite <- function(listOfFiles,
                   beginTime=0,
                   endTime=0,
                   windowShift=5, 
                   minF=40, 
                   maxF=500, 
                   voicing.threshold=0.3,
                   use.gaussian=FALSE,
                   interpolation.period=1,
                   gaussian.seed=1,
                   explicitExt="xte",
                   outputDirectory=NULL,
                   toFile=TRUE, 
                   conda.env=NULL){
  
  if(!is.null(conda.env) && !conda.env %in%  reticulate::conda_list()$name){
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
    origSoundFile <- normalizePath(fileBeginEnd[i, "listOfFiles"],mustWork = TRUE)
    
    beginTime <- fileBeginEnd[i, "beginTime"]
    endTime <- fileBeginEnd[i, "endTime"]
    
 
    py$soundFile <- reticulate::r_to_py(origSoundFile)
    py$ws <- reticulate::r_to_py(windowShift/1000) # reaper takes seconds
    py$fMax <- reticulate::r_to_py(maxF)
    py$fMin <- reticulate::r_to_py(minF)
    py$bt <- reticulate::r_to_py(beginTime)
    py$et <- reticulate::r_to_py(endTime)
    py$vt <- reticulate::r_to_py(voicing.threshold)
    py$gaussian <- reticulate::r_to_py(use.gaussian)
    py$interpp <- reticulate::r_to_py(interpolation.period)
    py$gs <- reticulate::r_to_py(gaussian.seed)
        
    "import numpy as np\
import gc\
import pysptk as sp\
import librosa as lr\
import pyreaper\
if et > 0:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt, duration = et - bt)\
else:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt)\

hs = ws * fs\
pitch_swipe = sp.swipe(x.astype(np.float64), fs=fs, hopsize=hs, min=fMin, max=fMax, otype=\"pitch\",threshold=vt)\
exct = sp.excite(pitch_swipe, hopsize=hs)\
del x\
gc.collect()" -> script
    Residual_symetry_string <- reticulate::py_suppress_warnings( reticulate::py_run_string(script))
    
    inTable <- data.frame( "excitation" = py$exct)
    
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16")
    #Use the time separation between second and pitch measurement time stamps to compute a sample frequency.
    
    sampleRate <-  1/ windowShift * 1000
    attr(outDataObj, "sampleRate") <- sampleRate
    
    attr(outDataObj, "origFreq") <-  as.numeric(py$fs) 
    startTime <- 1/sampleRate
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
    class(outDataObj) = "AsspDataObj"
    
    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    # excitation track
    extTable <- inTable %>%
      dplyr::select(excitation) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    noExtValues <- nrow(extTable)
    names(extTable) <- NULL
    outDataObj = addTrack(outDataObj, "excitation", as.matrix(extTable[,1]), "INT16")
    
    
    ## Apply fix from Emu-SDMS manual
    ##https://raw.githubusercontent.com/IPS-LMU/The-EMU-SDMS-Manual/master/R/praatToFormants2AsspDataObj.R
    
    # add missing values at the start as Praat sometimes
    # has very late start values which causes issues
    # in the SSFF file format as this sets the startRecord
    # depending on the start time of the first sample
    if( startTime > (1/sampleRate) ){
      
      nr_of_missing_samples = as.integer(floor(startTime / (1/sampleRate)))
      
      missing_excitation_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$excitation))

      
      
      # prepend values
      outDataObj$excitation = rbind(missing_excitation_vals, outDataObj$excitation)

      
      
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

attr(excite,"ext") <-  c("xte") 
attr(excite,"tracks") <-  c("excitation")
attr(excite,"outputType") <-  c("SSFF")



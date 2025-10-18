#' Estimates aperiodicity of a speech signal (DEPRECATED)
#'
#' @description **DEPRECATED**: This function is deprecated. Please use [d4c()] instead,
#'   which provides a faster C++ implementation with the same functionality.
#'
#' Aperiodicity is estimated using function for assessing band-aperiodicities
#' using the Definitive Decomposition Derived Dirt-Cheap (D4C) algorithm
#' \insertCite{Morise.2016.10.1016/j.specom.2016.09.001}{superassp} implemented
#' in the WORLD vocoder
#' \insertCite{MORISE.2016.10.1587/transinf.2015edp7457}{superassp}. The the [dio][DIO]
#' \insertCite{morise2010rapid}{superassp} pitch algorithm is used to calculate
#' the periodic component.
#'
#'
#' @inheritParams dio
#'
#' @return An SSFF track object containing two tracks (f0 and corr) that are
#'   either returned (toFile == FALSE) or stored on disk.
#' @references \insertAllCited{}
#'
#'
#'   
aperiodicities<- function(listOfFiles,
               beginTime=0,
               endTime=0,
               windowShift=5,
               minF=70, 
               maxF=200, 
               voiced_voiceless_threshold = 0.01,
               explicitExt="wap",
               outputDirectory=NULL,
               toFile=TRUE){
  
  .Deprecated("d4c", 
              msg = "aperiodicities() is deprecated. Please use d4c() for faster C++ implementation.")
  
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
    py$voiced_voiceless_threshold <- reticulate::r_to_py(voiced_voiceless_threshold)
    
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
_f0, t = pw.dio(x,\
	fs,\
	f0_floor=minF,\
	f0_ceil=maxF,\
	frame_period=windowShift,\
	allowed_range=voiced_voiceless_threshold)\
f0 = pw.stonemask(x, _f0, t, fs)\
ap = pw.d4c(x, f0, t, fs)\
apc = pw.code_aperiodicity(ap,fs)\
del x\
del _f0\
del ap\
del f0\
gc.collect()")
    
    aperiodicityTable <- as.data.frame( py$apc) %>% 
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
    attr(outDataObj, "endRecord") <- as.integer(nrow(aperiodicityTable))
    class(outDataObj) = "AsspDataObj"
    
    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    
    
    noAperiodicitiesValues <- nrow(aperiodicityTable)
    names(aperiodicityTable) <- NULL
    
    outDataObj = addTrack(outDataObj, "aperiod", as.matrix(aperiodicityTable), "REAL32")
    
    
    ## Apply fix from Emu-SDMS manual
    ##https://raw.githubusercontent.com/IPS-LMU/The-EMU-SDMS-Manual/master/R/praatToFormants2AsspDataObj.R
    
    # add missing values at the start as Praat sometimes
    # has very late start values which causes issues
    # in the SSFF file format as this sets the startRecord
    # depending on the start time of the first sample
    if( startTime > (1/sampleRate) ){
      
      nr_of_missing_samples = as.integer(floor(startTime / (1/sampleRate)))
      
      missing_aperiodicity_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$aperiodicities))
      
      
      # prepend values
      outDataObj$aperiodicities = rbind(missing_aperiodicity_vals, outDataObj$aperiodicities)   
      
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



attr(aperiodicities,"ext") <-  c("wap") 
attr(aperiodicities,"tracks") <-  c("aperiod")
attr(aperiodicities,"outputType") <-  c("SSFF")



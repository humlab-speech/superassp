
#' Extract pitch marks using the REAPER algoritm
#'
#' Robust Epoch And Pitch EstimatoR (REAPER) algorithm
#' \insertCite{talkin2019reaper}{superassp} uses an EpochTracker class to
#' simultaneously estimate the location of voiced-speech "epochs" or glottal
#' closure instants (GCI), voicing state (voiced or unvoiced) and fundamental
#' frequency (F0 or "pitch"). This function returns the voicing state of each
#' windowed `windowShift` (ms) portion of the signal.
#'
#' DC bias and low-frequency noise are removed by high-pass filtering, and the
#' signal is converted to floating point. If the input is known to have phase
#' distortion that is impacting tracker performance, a Hilbert transform,
#' optionally done at this point, may improve performance.
#'
#'
#' @details The function uses the python library `pyreaper` combined with the R
#'   package [reticulate] to compute the tracks, and the user therefore has to
#'   make sure that `pyreaper` and python is available on the machine. It is
#'   recommended to set up an anaconda ("conda") environment for the superassp
#'   library, like this:
#'
#'   ``` 
#'   conda create conda create --prefix -n pysuperassp python=3.8 
#'   conda activate pysuperassp 
#'   pip install librosa
#'   pip install pyreaper 
#'   #Not used by
#'   this function but by other functions in this package 
#'   pip install pysptk 
#'   ```
#' @references \insertAllCited{}
#'
#' @inheritParams reaper
#'
#'
#' @return An SSFF track object containing two tracks (f0 and corr) that are
#'   either returned (toFile == FALSE) or stored on disk.
#'
#' @export
#' 
reaper_pm <- function(listOfFiles,
                   beginTime=0,
                   endTime=0,
                   windowShift=10, 
                   minF=40, 
                   maxF=500, 
                   unvoiced_cost=0.9,
                   high.pass=TRUE,
                   hilbert.transform=FALSE,
                   explicitExt="rpm",
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
    py$uc <- reticulate::r_to_py(unvoiced_cost)
    py$hp <- reticulate::r_to_py(high.pass)
    py$ht <- reticulate::r_to_py(hilbert.transform)
    
    "import numpy as np\
import gc\
import pysptk as sp\
import librosa as lr\
import pyreaper\
if et > 0:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt, duration = et - bt)\
else:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt)\

raw_x = x * 2**15\
int_x = raw_x.astype(np.int16)\

pm_times, pm, f0_times, f0, corr = pyreaper.trk_reaper(x=int_x, fs=fs, minf0 = fMin, maxf0 = fMax, do_high_pass=hp, do_hilbert_transform= ht,  frame_period=ws, inter_pulse=ws, unvoiced_cost =uc)\
del x\
del raw_x\
del int_x\
gc.collect()" -> script
    Residual_symetry_string <- reticulate::py_capture_output( reticulate::py_run_string(script),type = "stderr")
    
    inTable <- data.frame( "time" = py$pm_times,
                           "pm" = py$pm)
    
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16")
    #Use the time separation between second and pitch measurement time stamps to compute a sample frequency.
    
    sampleRate <-  1/ windowShift * 1000
    attr(outDataObj, "sampleRate") <- sampleRate
    
    attr(outDataObj, "origFreq") <-  as.numeric(py$fs) 
    startTime <- 1/sampleRate
    attr(outDataObj, "startTime") <- as.numeric(py$f0_times[[1]])
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
    class(outDataObj) = "AsspDataObj"
    
    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    # f0 track
    pmTable <- inTable %>%
      dplyr::select(pm) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    noPMValues <- nrow(pmTable)
    names(pmTable) <- NULL
    outDataObj = addTrack(outDataObj, "pm", as.matrix(pmTable[,1]), "INT16")
    
    
    ## Apply fix from Emu-SDMS manual
    ##https://raw.githubusercontent.com/IPS-LMU/The-EMU-SDMS-Manual/master/R/praatToFormants2AsspDataObj.R
    
    # add missing values at the start as Praat sometimes
    # has very late start values which causes issues
    # in the SSFF file format as this sets the startRecord
    # depending on the start time of the first sample
    if( startTime > (1/sampleRate) ){
      
      nr_of_missing_samples = as.integer(floor(startTime / (1/sampleRate)))
      
      missing_pm_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$pm))

      
      
      # prepend values
      outDataObj$pm = rbind(missing_pm_vals, outDataObj$pm)
      
      
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



attr(reaper_pm,"ext") <-  c("rpm") 
attr(reaper_pm,"tracks") <-  c("pm")
attr(reaper_pm,"outputType") <-  c("SSFF")



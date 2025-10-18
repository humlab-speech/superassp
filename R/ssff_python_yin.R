#' Pitch detection using the YIN metod
#'
#' This function applies the YIN
#' \insertCite{Cheveigné.2002.10.1121/1.1458024}{superassp} method to estimate
#' the fundamental frequency.
#' 
#' This function calls the librosa \insertCite{brian_mcfee_2022_6097378}{superassp} Python library to load the audio data an
#' make pitch related estimates.
#'
#' @inheritParams rapt
#' @param trough_threshold The absolute threshold for peak estimation.
#' @param center Should analysis windows be centered around the time of the
#'   window (`TRUE`, the default) or should the window be considered to have
#'   started at the indicated time point (`FALSE`).
#' @param pad_mode The mode in which padding occurs. Ignored if `center` is not
#'   `TRUE`. Padding occurs in the python library librosa, and the user should
#'   therefore consult the manual of the NumPy library function 
#'   [numpy.pad](https://numpy.org/doc/stable/reference/generated/numpy.pad.html#numpy-pad) for other options.
#'
#' @return
#'  An SSFF track object containing two tracks (f0 and pitch) that are either returned (toFile == FALSE) or stored on disk.
#'  
#' @export
#' 
#' @seealso pyin
#'
#' @references 
#' \insertAllCited{}
#' 
#' 
trk_yin <- function(listOfFiles,
                  beginTime=0,
                  endTime=0,
                  windowShift=5,
                windowSize=30,
                  minF=70, 
                  maxF=200, 
                trough_threshold=0.1,
                center=TRUE,
                pad_mode="constant",
                  explicitExt="yip",
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

    # Initialize Python environment
    py <- reticulate::import_main()

    py$soundFile <- reticulate::r_to_py(origSoundFile)
    py$windowShift <- reticulate::r_to_py(windowShift)
    py$windowSize <- reticulate::r_to_py(windowSize)
    py$fMax <- reticulate::r_to_py(maxF)
    py$fMin <- reticulate::r_to_py(minF)
    py$beginTime <- reticulate::r_to_py(beginTime)
    py$endTime <- reticulate::r_to_py(endTime)
    py$center <- reticulate::r_to_py(center)
    py$trough_threshold <- reticulate::r_to_py(trough_threshold)
    py$pad_mode <- reticulate::r_to_py(pad_mode)

    reticulate::py_run_string("import librosa\
import gc\
import numpy as np\
duration = None\
if endTime > (windowSize/1000) and (endTime-beginTime) >= (windowSize/1000) :\
	duration =  (endTime - beginTime)\
\
waveform, fs = librosa.load(soundFile,\
	offset=beginTime,\
	duration=duration,\
	mono=True)\
\
hop_length = librosa.time_to_samples(windowShift/1000,sr=fs)\
win_length = librosa.time_to_samples(windowSize/1000,sr=fs)\
frame_length= win_length * 2\
\
pitch = librosa.trk_yin(waveform,\
	fmin=fMin,\
	fmax=fMax,\
	hop_length = hop_length,\
	frame_length = frame_length,\
	win_length = win_length,\
	sr=fs,\
	trough_threshold=trough_threshold,\
	center=center,\
	pad_mode=pad_mode)\
del waveform\
gc.collect()")
    
    inTable <- data.frame( "f0" = py$pitch)
    

    startTime = windowSize/2
    
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
    
    # Cross-correlation track
    f0Table <- inTable %>%
      dplyr::select(f0) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    nof0Values <- nrow(f0Table)
    names(f0Table) <- NULL
    outDataObj = addTrack(outDataObj, "f0", as.matrix(f0Table[,1]), "INT16")

    
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



attr(trk_yin,"ext") <-  c("yip") 
attr(trk_yin,"tracks") <-  c("f0")
attr(trk_yin,"outputType") <-  c("SSFF")



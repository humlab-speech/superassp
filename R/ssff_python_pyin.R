

#' Estimate pitch using the probabilistic YIN algorithm
#'
#' The probabilistic YIN algorithm
#' \insertCite{Mauch.2014.10.1109/icassp.2014.6853678}{superassp} is an
#' extension of [YIN][superassp::yin] \insertCite{Cheveigné.2002.10.1121/1.1458024}{superassp} that considers multiple pitch
#' candidates in a hidden Markov model that is Viterbi-decoded to deduce the
#' final pitch estimate. The function also returns a track encoding whether the
#' track was considered voiced or not, and a track containing the probability of
#' voicing in the analysis frame.
#'
#' This function calls the librosa \insertCite{brian_mcfee_2022_6097378}{superassp} Python library to load the audio data an
#' make pitch related estimates.
#'
#' @inheritParams yin
#' @param max_transition_rate The maximum pitch transition rate in octaves per second.
#' @param beta_parameters The shape parameters for the beta distribution prior over thresholds.
#' @param boltzmann_parameter The shape parameter for the Boltzmann distribution prior over troughs. Larger values will assign more mass to smaller periods.
#' @param resolution The resolution of the pitch bins. 0.01 corresponds to cents.
#' @param thresholds The number of thresholds for peak estimation.
#' @param switch_probability The probability of switching from voiced to unvoiced or vice versa.
#' @param no_trough_probability The maximum probability to add to global minimum if no trough is below threshold.
#'
#' @return
#'  An SSFF track object containing two tracks (f0 and pitch) that are either returned (toFile == FALSE) or stored on disk.
#'  
#' @export
#'
#'
#' @references 
#' \insertAllCited{}
#' 
pyin <- function(listOfFiles,
                beginTime=0,
                endTime=0,
                windowShift=5,
                windowSize=30,
                minF=70, 
                maxF=200, 
                #trough_threshold=0.1,
                max_transition_rate=35.92,
                beta_parameters=c(2,18),
                center=TRUE,
                boltzmann_parameter = 2,
                resolution = 0.1,
                thresholds=100,
                switch_probability = 0.01,
                no_trough_probability = 0.01,
                pad_mode="constant",
                explicitExt="pyp",
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
    py$thresholds <- reticulate::r_to_py(as.integer(thresholds))
    py$pad_mode <- reticulate::r_to_py(pad_mode)
    py$max_transition_rate <- reticulate::r_to_py(max_transition_rate)
    py$switch_probability <- reticulate::r_to_py(switch_probability)
    py$no_trough_probability <- reticulate::r_to_py(no_trough_probability)
    py$beta_parameters <- reticulate::r_to_py(beta_parameters)
    py$resolution <- reticulate::r_to_py(resolution)

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
pitch, voiced_flag, voiced_prob = librosa.pyin(waveform, 
	fmin=fMin, 
	fmax=fMax,
	hop_length = hop_length,
	frame_length = frame_length,
	win_length = win_length,
	sr=fs,
	n_thresholds=thresholds, 
	beta_parameters=beta_parameters, 
	boltzmann_parameter=2, 
	resolution=resolution, 
	max_transition_rate=max_transition_rate, 
	switch_prob=switch_probability, 
	no_trough_prob=no_trough_probability, 
	fill_na=0, 
	center=center, 
	pad_mode=pad_mode)\
del waveform\
gc.collect()")
    
    inTable <- data.frame( "f0" = py$pitch,
                           "voiced"=ifelse(py$voiced_flag,1,0),
                           "vprob"=py$voiced_prob)
    
    startTime = windowSize/2
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16","INT16","REAL32")
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
    
    # Voiced / unvoiced (1,0) track
    voicedTable <- inTable %>%
      dplyr::select(voiced) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    noVoicedValues <- nrow(voicedTable)
    names(voicedTable) <- NULL
    
    outDataObj = addTrack(outDataObj, "voiced", as.matrix(voicedTable[,1]), "INT16")
    
    # Voicing probability track
    
    vprobTable <- inTable %>%
      dplyr::select(vprob) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    noVprobValues <- nrow(vprobTable)
    names(vprobTable) <- NULL
    
    outDataObj = addTrack(outDataObj, "vprob", as.matrix(vprobTable[,1]), "REAL32") 
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
      missing_voiced_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$voiced))
      missing_vprob_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$vprob))
      
      # prepend values
      outDataObj$f0 = rbind(missing_f0_vals, outDataObj$f0)
      outDataObj$voiced = rbind(missing_voiced_vals, outDataObj$voiced)     
      outDataObj$vprob = rbind(missing_vprob_vals, outDataObj$vprob)     
      
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



attr(pyin,"ext") <-  c("pyp") 
attr(pyin,"tracks") <-  c("f0","voiced","vprob")
attr(pyin,"outputType") <-  c("SSFF")

#' Compute f0 using the Harvest algorithm



#' Compute f0 using the SWIPE algorithm
#' 
#' This function takes a sound file and computes f$_0$ and an estimate of pitch 
#' using the Sawtooth Waveform Inspired Pitch Estimator (SWIPE) algorithm \insertCite{Camacho.2008.10.1121/1.2951592}{superassp}. 
#' 
#' @details 
#' The implementation of SWIPE in the Speech Signal Processing Toolkit (SPTK) \insertCite{sptkspeech}{superassp} is used, and called via its Python interface and the [retiulate] R package to compute the signal track.
#' Therefore, the user will have to make sure that a python environment is present and can be attached by the [reticulate]. An anaconda environment is recommended, and can set up by the user by a setup procedure that involve at least these commands:
#' 
#' ```
#' conda create conda create --prefix -n pysuperassp python=3.8 
#' conda activate pysuperassp
#' pip install librosa
#' pip install pysptk
#' ```
#' to make the functionality that this function requires available. 
#'  
#'
#' @param listOfFiles A vector of file paths to wav files.
#' @param beginTime The start time of the section of the sound file that should be processed.
#' @param endTime The end time of the section of the sound file that should be processed.
#' @param windowShift  The measurement interval (frame duration), in seconds.
#' @param minF Candidate f0 frequencies below this frequency will not be considered. 
#' @param maxF Candidates above this frequency will be ignored.
#' @param voicing.threshold Voice/unvoiced threshold. Default is 0.3.
#' @param conda.env The name of the conda environment in which Python and its
#'   required packages are stored. Please make sure that you know what you are
#'   doing if you change this. Defaults to `NULL`, which means that the default enviroment or the environment set in the 
#'   `RETICULATE_PYTHON` environment variable will be used.
#' @inheritParams praat_formant_burg
#'
#' @return
#'  An SSFF track object containing two tracks (f0 and pitch) that are either returned (toFile == FALSE) or stored on disk.
#' 
#' @export
#' 
#' 
#' @references 
#'   \insertAllCited{}
#'   
swipe <- function(listOfFiles,
                  beginTime=0,
                  endTime=0,
                  windowShift=5, 
                  minF=70, 
                  maxF=200, 
                  voicing.threshold=0.3,
                  explicitExt="swi",
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
    py$ws <- reticulate::r_to_py(windowShift)
    py$fMax <- reticulate::r_to_py(maxF)
    py$fMin <- reticulate::r_to_py(minF)
    py$bt <- reticulate::r_to_py(beginTime)
    py$et <- reticulate::r_to_py(endTime)
    
    
    reticulate::py_run_string("import numpy as np\
import gc\
import pysptk as sp\
import librosa as lr\
if et > 0:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt, duration = et - bt)\
else:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt)\

f0_swipe = sp.swipe(x.astype(np.float64), fs=fs, hopsize=ws / 1000 * fs, min=fMin, max=fMax, otype=\"f0\")\
pitch_swipe = sp.swipe(x.astype(np.float64), fs=fs, hopsize=ws / 1000 * fs, min=fMin, max=fMax, otype=\"pitch\")\
del x\
gc.collect()")
    
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


#' Compute f0 using the RAPT algorithm
#' 
#' This function takes a sound file and computes f$_0$ and an estimate of pitch 
#' using the "A robust algorithm for pitch tracking" (RAPT) algorithm \insertCite{talkin1995robust}{superassp}. 
#' 
#' @details 
#' The implementation of RAPT in the Speech Signal Processing Toolkit (SPTK) \insertCite{sptkspeech}{superassp} is used, and called via its Python interface and the [retiulate] R package to compute the signal track.
#' Therefore, the user will have to make sure that a python environment is present and can be attached by the [reticulate]. An anaconda environment is recommended, and can set up by the user by a setup procedure that involve at least these commands:
#' 
#' ```
#' conda create conda create --prefix -n pysuperassp python=3.8 
#' conda activate pysuperassp
#' pip install librosa
#' pip install pysptk
#' #Not used by this function but by other functions in this package
#' pip install pyreaper 
#' ```
#' to make the functionality that this function requires available. 
#'  
#'
#' @param listOfFiles A vector of file paths to wav files.
#' @param beginTime The start time of the section of the sound file that should be processed.
#' @param endTime The end time of the section of the sound file that should be processed.
#' @param windowShift  The measurement interval (frame duration), in seconds.
#' @param minF Candidate f0 frequencies below this frequency will not be considered. 
#' @param maxF Candidates above this frequency will be ignored.
#' @param voicing.threshold Voice/unvoiced threshold. Default is 0.3.
#' @param conda.env The name of the conda environment in which Python and its required packages are stored. Please make sure that you know what you are doing if you change this.
#' @inheritParams praat_formant_burg
#'
#' @return
#'  An SSFF track object containing two tracks (f0 and pitch) that are either returned (toFile == FALSE) or stored on disk.
#' 
#' @export
#' 
#' 
#' @references 
#'   \insertAllCited{}
#'   
rapt <- function(listOfFiles,
                  beginTime=0,
                  endTime=0,
                  windowShift=5, 
                  minF=70, 
                  maxF=200, 
                  voicing.threshold=0.3,
                  explicitExt="swi",
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
    py$ws <- reticulate::r_to_py(windowShift)
    py$fMax <- reticulate::r_to_py(maxF)
    py$fMin <- reticulate::r_to_py(minF)
    py$bt <- reticulate::r_to_py(beginTime)
    py$et <- reticulate::r_to_py(endTime)
    
    
    reticulate::py_run_string("import numpy as np\
import gc\
import pysptk as sp\
import librosa as lr\
if et > 0:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt, duration = et - bt)\
else:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt)\

f0_rapt = sp.rapt(x.astype(np.float64), fs=fs, hopsize=ws / 1000 * fs, min=fMin, max=fMax, otype=\"f0\")\
pitch_rapt = sp.rapt(x.astype(np.float64), fs=fs, hopsize=ws / 1000 * fs, min=fMin, max=fMax, otype=\"pitch\")\
del x\
gc.collect()")
    
    inTable <- data.frame( "f0" = py$f0_rapt,
                           "pitch"=py$pitch_rapt)
    
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



attr(rapt,"ext") <-  c("rpt") 
attr(rapt,"tracks") <-  c("f0","pitch")
attr(rapt,"outputType") <-  c("SSFF")



#' Extract f0 tracks using the REAPER algoritm
#'
#' Robust Epoch And Pitch EstimatoR (REAPER) algorithm
#' \insertCite{talkin2019reaper}{superassp} uses an EpochTracker class to
#' simultaneously estimate the location of voiced-speech "epochs" or glottal
#' closure instants (GCI), voicing state (voiced or unvoiced) and fundamental
#' frequency (F0 or "pitch"). The local (instantaneous) f0 is defined as the
#' inverse of the time between successive GCI. This function returns the f0 and
#' normalized CGI cross-correlation in each windowed `windowShift` (ms) portion
#' of the signal.
#'
#' DC bias and low-frequency noise are removed by high-pass filtering, and the
#' signal is converted to floating point. If the input is known to have phase
#' distortion that is impacting tracker performance, a Hilbert transform,
#' optionally done at this point, may improve performance.
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
#'   
#' 
#'
#' @param unvoiced_cost Set the cost for unvoiced segments. Default is 0.9, the higher the value the more f0 estimates in noise.
#' @param high.pass Perform high-pass filtering to remove DC and low-frequency noise?
#' @param hilbert.transform Remove phase distortion using Hilbert transform?
#' @inheritParams swipe
#'
#'
#' @return An SSFF track object containing two tracks (f0 and corr) that are
#'   either returned (toFile == FALSE) or stored on disk.
#' @references \insertAllCited{}
#' 
#' @export
#' 
reaper <- function(listOfFiles,
                 beginTime=0,
                 endTime=0,
                 windowShift=5, 
                 minF=40, 
                 maxF=500, 
                 unvoiced_cost=0.9,
                 high.pass=TRUE,
                 hilbert.transform=FALSE,
                 explicitExt="rp0",
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

pm_times, pm, f0_times, f0, corr = pyreaper.reaper(x=int_x, fs=fs, minf0 = fMin, maxf0 = fMax, do_high_pass=hp, do_hilbert_transform= ht,  frame_period=ws, inter_pulse=ws, unvoiced_cost =uc)\
del x\
del raw_x\
del int_x\
gc.collect()" -> script
    Residual_symetry_string <- reticulate::py_suppress_warnings( reticulate::py_run_string(script))
    
    inTable <- data.frame( "time" = py$f0_times,
                           "f0" = py$f0,
                           "corr"=py$corr)
    
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16", "INT16")
    #Use the time separation between second and pitch measurement time stamps to compute a sample frequency.
    
    sampleRate <-  1/ windowShift * 1000
    attr(outDataObj, "sampleRate") <- sampleRate
    
    attr(outDataObj, "origFreq") <-  as.numeric(py$fs) 
    startTime <- 1/sampleRate
    attr(outDataObj, "startTime") <- as.numeric(py$f0_times[[1]])
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
    class(outDataObj) = "AsspDataObj"
    
    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    # f0 track
    f0Table <- inTable %>%
      dplyr::select(f0) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    nof0Values <- nrow(f0Table)
    names(f0Table) <- NULL
    outDataObj = wrassp::addTrack(outDataObj, "f0", as.matrix(f0Table[,1]), "INT16")
    
    # Correlation track 
    corrTable <- inTable %>%
      dplyr::select(corr) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    noCorrValues <- nrow(corrTable)
    names(corrTable) <- NULL
    outDataObj = wrassp::addTrack(outDataObj, "corr", as.matrix(corrTable[,1]), "INT16")
    

  
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
      missing_corr_vals = matrix(0,
                                  nrow = nr_of_missing_samples,
                                  ncol = ncol(outDataObj$corr))

      
      # prepend values
      outDataObj$f0 = rbind(missing_f0_vals, outDataObj$f0)
      outDataObj$corr = rbind(missing_pitch_vals, outDataObj$corr)

      
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



attr(reaper,"ext") <-  c("rp0") 
attr(reaper,"tracks") <-  c("f0","corr")
attr(reaper,"outputType") <-  c("SSFF")

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

pm_times, pm, f0_times, f0, corr = pyreaper.reaper(x=int_x, fs=fs, minf0 = fMin, maxf0 = fMax, do_high_pass=hp, do_hilbert_transform= ht,  frame_period=ws, inter_pulse=ws, unvoiced_cost =uc)\
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
    
    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    # f0 track
    pmTable <- inTable %>%
      dplyr::select(pm) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    noPMValues <- nrow(pmTable)
    names(pmTable) <- NULL
    outDataObj = wrassp::addTrack(outDataObj, "pm", as.matrix(pmTable[,1]), "INT16")
    
    
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



attr(reaper_pm,"ext") <-  c("rpm") 
attr(reaper_pm,"tracks") <-  c("pm")
attr(reaper_pm,"outputType") <-  c("SSFF")


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
    
    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    # excitation track
    extTable <- inTable %>%
      dplyr::select(excitation) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    noExtValues <- nrow(extTable)
    names(extTable) <- NULL
    outDataObj = wrassp::addTrack(outDataObj, "excitation", as.matrix(extTable[,1]), "INT16")
    
    
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

attr(excite,"ext") <-  c("xte") 
attr(excite,"tracks") <-  c("excitation")
attr(excite,"outputType") <-  c("SSFF")



#' Estimate pitch using the Kaldi modifies version of RAPT
#'
#' The algorithm used is a version of the [RAPT][rapt] algorithm
#' that considers voicing also in voiceless frames and conputes a
#' Normalized Cross Correlation Function (NCCF) that can be used to
#' estimate the probability of voicing \insertCite{Ghahremani.2014.10.1109/icassp.2014.6854049}{superassp}.
#'
#' The function calls the [torchaudio](https://github.com/pytorch/audio) \insertCite{yang2021torchaudio}{superassp} library to do the pitch estimates and therefore
#' relies on it being present in a properly set up python environment to work. Please refer to the [torchaudio manual](https://pytorch.org/audio/main/generated/torchaudio.functional.compute_kaldi_pitch.html) for further information.
#'
#'
#' @inheritParams rapt
#' @param soft_min_f0 (float, optional) – Minimum f0, applied in soft way, must not exceed min-f0 (default: 10.0)
#' @param penalty_factor Cost factor for fO change. (default: 0.1)
#' @param lowpass_cutoff Cutoff frequency for LowPass filter (Hz) (default: 1000)
#' @param resample_frequency Frequency that we down-sample the signal to. Must be more than twice `lowpass_cutoff`. (default: 4000)
#' @param delta_pitch  Smallest relative change in pitch that our algorithm measures. (default: 0.005)
#' @param nccf_ballast Increasing this factor reduces NCCF for quiet frames (default: 7000)
#' @param lowpass_filter_width  Integer that determines filter width of lowpass filter, more gives sharper filter. (default: 1)
#' @param psample_filter_width  Integer that determines filter width when upsampling NCCF. (default: 5)
#' @param max_frames_latency  Maximum number of frames of latency that we allow pitch tracking to introduce into the feature processing (affects output only if `frames_per_chunk` > 0 and `simulate_first_pass_online`=`TRUE`) (default: 0)
#' @param frames_per_chunk The number of frames used for energy normalization. (default: 0)
#' @param simulate_first_pass_online  If true, the function will output features that correspond to what an online decoder would see in the first pass of decoding – not the final version of the features, which is the default. (default: `FALSE`) Relevant if `frames_per_chunk > 0`.
#' @param recompute_frame  Only relevant for compatibility with online pitch extraction. A non-critical parameter; the frame at which we recompute some of the forward pointers, after revising our estimate of the signal energy. Relevant if `frames_per_chunk > 0`. (default: 500)
#' @param snip_edges  If this is set to false, the incomplete frames near the ending edge won’t be snipped, so that the number of frames is the file size divided by the `windowShift`. This makes different types of features give the same number of frames. (default: True)
#'
#' @return An SSFF track object containing two tracks (f0 and nccf) that are
#'   either returned (toFile == FALSE) or stored on disk.
#'
#' @seealso rapt
#'
#' @references \insertAllCited{}
#' @export


kaldi_pitch <- function(listOfFiles,
                 beginTime=0,
                 endTime=0,
                 windowShift=5,
                 windowSize=25,
                 minF=70,
                 maxF=200,
                 softMinF0 = 10.0,
                 voiced_voiceless_cost = 0.10,
                 owpass_cutoff =1000.0,
                 resample_frequency = 4000.0,
                 deltaChange = 0.005,
                 nccfBallast = 7000,
                 lowpass_cutoff =1000,
                 lowpass_filter_width = 1,
                 upsample_filter_width = 5,
                 max_frames_latency = 0,
                 frames_per_chunk = 0,
                 simulate_first_pass_online = FALSE,
                 recompute_frame=500,
                 snip_edges = TRUE,
                 explicitExt="kap",
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
    py$windowShift <- reticulate::r_to_py(as.double(windowShift))
    py$windowSize <- reticulate::r_to_py(as.double(windowSize))
    py$maxF <- reticulate::r_to_py(as.double(maxF))
    py$minF <- reticulate::r_to_py(as.double(minF))
    py$beginTime <- reticulate::r_to_py(as.double(beginTime))
    py$endTime <- reticulate::r_to_py(as.double(endTime))
    py$softMinF0 <- reticulate::r_to_py(as.double(softMinF0))
    py$voiced_voiceless_cost <- reticulate::r_to_py(as.double(voiced_voiceless_cost))
    py$lowpass_cutoff <- reticulate::r_to_py(as.double(lowpass_cutoff))
    py$resample_frequency <- reticulate::r_to_py(as.double(resample_frequency))
    py$deltaChange <- reticulate::r_to_py(as.double(deltaChange))
    py$nccfBallast <- reticulate::r_to_py(as.double(nccfBallast))
    py$lowpass_filter_width <- reticulate::r_to_py(as.integer(lowpass_filter_width))
    py$upsample_filter_width <- reticulate::r_to_py(as.integer(upsample_filter_width))
    py$max_frames_latency <- reticulate::r_to_py(as.integer(max_frames_latency))
    py$frames_per_chunk <- reticulate::r_to_py(as.integer(frames_per_chunk))
    py$simulate_first_pass_online <- reticulate::r_to_py(as.logical(simulate_first_pass_online))
    py$recompute_frame <- reticulate::r_to_py(as.integer(recompute_frame))
    py$snip_edges <- reticulate::r_to_py(as.logical(snip_edges))


    reticulate::py_run_string("import torch\
import gc\
import torchaudio\
import torchaudio.functional as F \
import torchaudio.transforms as T \
import math \
\
metadata = torchaudio.info(soundFile) \
\
if  beginTime > 0 and endTime > beginTime : \
	startSample = math.floor(beginTime * metadata.sample_rate) \
else: \
	startSample = 0 \
\
if  endTime > 0 and endTime > beginTime :\
	endSample = math.ceil(endTime * metadata.sample_rate)\
	nSamples = endSample - startSample \
else: \
	nSamples = -1 \
\
SPEECH_WAVEFORM, SAMPLE_RATE = torchaudio.load(soundFile,frame_offset=startSample, num_frames= nSamples) \
\
pitch_feature = F.compute_kaldi_pitch(waveform=SPEECH_WAVEFORM,\
	sample_rate=SAMPLE_RATE, \
	frame_length= windowSize, \
	frame_shift= windowShift, \
	min_f0= minF, \
	max_f0= maxF, \
	soft_min_f0= softMinF0, \
	penalty_factor= voiced_voiceless_cost, \
	lowpass_cutoff= lowpass_cutoff, \
	resample_frequency= resample_frequency, \
	delta_pitch= deltaChange, \
	nccf_ballast= nccfBallast, \
	lowpass_filter_width= lowpass_filter_width, \
	upsample_filter_width= upsample_filter_width, \
	max_frames_latency= max_frames_latency, \
	frames_per_chunk= frames_per_chunk, \
	simulate_first_pass_online= simulate_first_pass_online, \
	recompute_frame= recompute_frame, \
	snip_edges=snip_edges)\
pitch, nccf = pitch_feature[..., 1], pitch_feature[..., 0] \
nppitch = pitch.numpy()\
npnccf = nccf.numpy()\
end_time = SPEECH_WAVEFORM.shape[1] / SAMPLE_RATE\
del pitch\
del nccf\
del pitch_feature\
del SPEECH_WAVEFORM\
gc.collect()")


    inTable <- data.frame( "f0" = as.vector(py$nppitch),
                           "nccf"=as.vector(py$npnccf))

    startTime = windowShift

    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16", "REAL32")
    #Use the time separation between second and pitch measurement time stamps to compute a sample frequency.

    sampleRate <-  1/ windowShift * 1000
    attr(outDataObj, "sampleRate") <- sampleRate

    attr(outDataObj, "origFreq") <-  as.numeric(py$SAMPLE_RATE)
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

    # Normalized Cross Correlation Function
    nccfTable <- inTable %>%
      dplyr::select(nccf) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))

    noNCCFValues <- nrow(nccfTable)
    names(nccfTable) <- NULL
    outDataObj = wrassp::addTrack(outDataObj, "nccf", as.matrix(nccfTable[,1]), "REAL32")



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
      missing_nccf_vals = matrix(0,
                                  nrow = nr_of_missing_samples,
                                  ncol = ncol(outDataObj$nccf))

      # prepend values
      outDataObj$f0 = rbind(missing_f0_vals, outDataObj$f0)
      outDataObj$nccf = rbind(missing_nccf_vals, outDataObj$nccf)


      # fix start time
      attr(outDataObj, "startTime") = startTime - nr_of_missing_samples * (1/sampleRate)
    }

    assertthat::assert_that(wrassp::is.AsspDataObj(outDataObj),
                            msg = "The AsspDataObj created by the kaldi_pitch function is invalid.")

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



attr(kaldi_pitch,"ext") <-  c("kap")
attr(kaldi_pitch,"tracks") <-  c("f0","nccf")
attr(kaldi_pitch,"outputType") <-  c("SSFF")




#' Compute pitch and periodicity using the CREPE pitch tracker
#' 
#' The CREPE \insertCite{Kim.2018.10.1109/icassp.2018.8461329}{superassp} applies a deep convolutional neural network directly on the time-domain waveform to find 
#' the fundamental frequency in a speech signal. Two versions of the models have been trained, one smaller yielding quicker results, and the full model which can be considerably 
#' more computationally intensive to apply.
#' 
#' T
#' 
#' @param listOfFiles A vector of file paths to wav files.
#' @param beginTime (Not implemented) The start time of the section of the sound file that should be processed.
#' @param endTime (Not implemented) The end time of the section of the sound file that should be processed.
#' @param windowShift  The measurement interval (frame duration), in seconds.
#' @param minF Candidate f0 frequencies below this frequency will not be considered. 
#' @param maxF Candidates above this frequency will be ignored.
#' @param voicing.threshold Voice/unvoiced threshold. Default is 0.21.
#' @param silence.threshold Frames that do not contain amplitudes above this threshold (relative to the global maximum amplitude), are probably silent.
#' @param model Use a fast ("tiny") model, or a more complete ("full") model to find pitch. The more complete model will take approximately 9-11 times longer to process the file. 
#' @param conda.env The name of the conda environment in which Python and its
#'   required packages are stored. Please make sure that you know what you are
#'   doing if you change this. Defaults to `NULL`, which means that the default enviroment or the environment set in the 
#'   `RETICULATE_PYTHON` environment variable will be used.
#' @inheritParams praat_formant_burg
#'
#' @return
#'  An SSFF track object containing two tracks (f0 and periodicity) that are either returned (toFile == FALSE) or stored on disk.
#' 
#' @export
#' 
#' 
#' @references 
#'   \insertAllCited{}
#'
crepe <- function(listOfFiles,
                 beginTime=0,
                 endTime=0,
                 windowShift=5,
                 windowSize=15,
                 minF=70, 
                 maxF=200, 
                 voicing.threshold=0.21,
                 silence.threshold=-60.0,
                 model=c("tiny","full"),
                 explicitExt="crp",
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
    py$windowShift <- reticulate::r_to_py(windowShift)
    py$windowSize <- reticulate::r_to_py(windowSize)
    py$model <- reticulate::r_to_py("tiny")    
    py$fmax <- reticulate::r_to_py(maxF)
    py$fmin <- reticulate::r_to_py(minF)
   # py$bt <- reticulate::r_to_py(beginTime)
    #py$et <- reticulate::r_to_py(endTime)
    py$silence_threshold <- reticulate::r_to_py(silence.threshold)
    py$voicing_threshold <- reticulate::r_to_py(voicing.threshold)   
    
    reticulate::py_run_string("import torchcrepe\
import gc\
import math\
audio, sr = torchcrepe.load.audio( soundFile )\
\
hop_length = int(sr / (1000.0 / windowShift))\
\
pitch, periodicity = torchcrepe.predict(audio, \
                           sr, \
                           hop_length, \
                           fmin, \
                           fmax, \
                           model, \
                           return_periodicity=True)\
\
win_length = math.ceil(windowSize / windowShift) \
\
periodicity = torchcrepe.filter.median(periodicity, win_length) \
\
periodicity = torchcrepe.threshold.Silence(silence_threshold)(periodicity, \
                                                 audio, \
                                                 sr, \
                                                 hop_length)\
\
pitch = torchcrepe.threshold.At(voicing_threshold)(pitch, periodicity)\
\
nppitch = torchcrepe.filter.mean(pitch, win_length).numpy()\
npperiodicity = periodicity.numpy()\
del pitch\
del periodicity\
del audio\
gc.collect()")
    
    inTable <- data.frame( "f0" = as.vector(py$nppitch),
                           "periodicity"=as.vector(py$npperiodicity))

    startTime = windowShift
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16", "REAL32")
    #Use the time separation between second and pitch measurement time stamps to compute a sample frequency.
    
    sampleRate <-  1/ windowShift * 1000
    attr(outDataObj, "sampleRate") <- sampleRate
    
    attr(outDataObj, "origFreq") <-  as.numeric(py$sr) 
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
    periodicityTable <- inTable %>%
      dplyr::select(periodicity) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    noPeriodicityValues <- nrow(periodicityTable)
    names(periodicityTable) <- NULL
    outDataObj = wrassp::addTrack(outDataObj, "periodicity", as.matrix(periodicityTable[,1]), "REAL32")
    
    
    
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
      
      missing_periodicity_vals = matrix(0,
                                  nrow = nr_of_missing_samples,
                                  ncol = ncol(outDataObj$periodicity))
      
      # prepend values
      outDataObj$f0 = rbind(missing_f0_vals, outDataObj$f0)
      outDataObj$periodicity = rbind(missing_periodicity_vals, outDataObj$periodicity)
      
      
      # fix start time
      attr(outDataObj, "startTime") = startTime - nr_of_missing_samples * (1/sampleRate)
    }
    
    assertthat::assert_that(wrassp::is.AsspDataObj(outDataObj),
                            msg = "The AsspDataObj created by the crepe function is invalid.")
    
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



attr(crepe,"ext") <-  c("crp") 
attr(crepe,"tracks") <-  c("f0","periodicity")
attr(crepe,"outputType") <-  c("SSFF")


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
yin <- function(listOfFiles,
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
pitch = librosa.yin(waveform,\
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



attr(yin,"ext") <-  c("yip") 
attr(yin,"tracks") <-  c("f0")
attr(yin,"outputType") <-  c("SSFF")




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
    
    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    # Pitch track
    f0Table <- inTable %>%
      dplyr::select(f0) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    nof0Values <- nrow(f0Table)
    names(f0Table) <- NULL

    outDataObj = wrassp::addTrack(outDataObj, "f0", as.matrix(f0Table[,1]), "INT16")
    
    # Voiced / unvoiced (1,0) track
    voicedTable <- inTable %>%
      dplyr::select(voiced) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    noVoicedValues <- nrow(voicedTable)
    names(voicedTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "voiced", as.matrix(voicedTable[,1]), "INT16")
    
    # Voicing probability track
    
    vprobTable <- inTable %>%
      dplyr::select(vprob) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    noVprobValues <- nrow(vprobTable)
    names(vprobTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "vprob", as.matrix(vprobTable[,1]), "REAL32") 
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



attr(pyin,"ext") <-  c("pyp") 
attr(pyin,"tracks") <-  c("f0","voiced","vprob")
attr(pyin,"outputType") <-  c("SSFF")

#' Compute f0 using the Harvest algorithm
#'
#' The DIO algorithm \insertCite{morise2010rapid}{superassp} was developed for
#' the WORLD vocoder
#' \insertCite{MORISE.2016.10.1587/transinf.2015edp7457}{superassp} and aims
#' provide a fast estimate of the f0 contour.
#'
#'
#' @inheritParams swipe
#' @param voiced.voiceless.threshold Threshold for voiced/unvoiced decision. Can
#'   be any value >= 0, but 0.02 to 0.2 is a reasonable range. Lower values will
#'   cause more frames to be considered unvoiced (in the extreme case of
#'   `threshold=0`, almost all frames will be unvoiced).
#'
#' @return An SSFF track object containing two tracks (f0 and corr) that are
#'   either returned (toFile == FALSE) or stored on disk.
#' @references \insertAllCited{}
#'
#' 
dio<- function(listOfFiles,
                 beginTime=0,
                 endTime=0,
                 windowShift=5,
                 minF=70, 
                 maxF=200, 
                 voiced_voiceless_threshold = 0.01,
                 explicitExt="wd0",
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
del _f0\
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
    
    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    # Pitch track
    f0Table <- inTable %>%
      dplyr::select(f0) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    nof0Values <- nrow(f0Table)
    names(f0Table) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "f0", as.matrix(f0Table[,1]), "INT16")
    

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



attr(dio,"ext") <-  c("wd0") 
attr(dio,"tracks") <-  c("f0")
attr(dio,"outputType") <-  c("SSFF")

#' Compute f0 using the Harvest algorithm
#' 
#' The Harvest algorithm \insertCite{Morise.2017.10.21437/interspeech.2017-68}{superassp} was developed for the WORLD vocoder \insertCite{MORISE.2016.10.1587/transinf.2015edp7457}{superassp} and aims to obtain a reliable F0 contour and reduce erroneously identified voice frames. 
#' 
#' The algorithm consists of two steps. In the first step, the algorithm uses fundamental component extraction by many band-pass filters with different center frequencies and obtains the basic f0 candidates from filtered signals. In the second step, basic f0 candidates are refined and scored by using the instantaneous frequency, and then several f0 candidates in each frame are estimated.
#'  
#' @inheritParams swipe
#'
#' @return An SSFF track object containing two tracks (f0 and corr) that are
#'   either returned (toFile == FALSE) or stored on disk.
#' @references 
#'   \insertAllCited{}
#'
#'
harvest<- function(listOfFiles,
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
f0, t = pw.harvest(x,\
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
    
    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    # Pitch track
    f0Table <- inTable %>%
      dplyr::select(f0) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    nof0Values <- nrow(f0Table)
    names(f0Table) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "f0", as.matrix(f0Table[,1]), "INT16")
    
    
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

attr(harvest,"ext") <-  c("hf0") 
attr(harvest,"tracks") <-  c("f0")
attr(harvest,"outputType") <-  c("SSFF")


npy_import <- function(listOfFiles,
               #inputExt="npy",
               explicitExt="npf",
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

    
    reticulate::py_run_string("
import numpy as np\
pitch = np.load(soundFile)")
    
    inTable <- data.frame( "f0" = py$pitch)
    
    startTime = windowSize/2
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16")
    #Use the time separation between second and pitch measurement time stamps to compute a sample frequency.
    
    sampleRate <-  1/ windowShift * 1000
    attr(outDataObj, "sampleRate") <- sampleRate
    
    attr(outDataObj, "origFreq") <-  as.numeric(py$fs) 
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
    class(outDataObj) = "AsspDataObj"
    
    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    # Pitch track
    f0Table <- inTable %>%
      dplyr::select(f0) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    nof0Values <- nrow(f0Table)
    names(f0Table) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "f0", as.matrix(f0Table[,1]), "INT16")
    
    
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



attr(npy_import,"ext") <-  c("npf") 
attr(npy_import,"tracks") <-  c("f0")
attr(npy_import,"outputType") <-  c("SSFF")


#' Compute f0 using the algorithm named Yet Another Algorithm for Pitch Tracking
#'
#' The Yet Another Algorithm for Pitch Tracking algorithm
#' \insertCite{Kasi.2002.10.1109/icassp.2002.5743729}{superassp} that computes
#' f0 using Normalized Cross Correlation (NCCF) and the work of Talkin
#' \insertCite{talkin1995robust}{superassp} in developing the [RAPT][rapt]
#' algorithm.
#'
#' The YAAPT algorithm processes the original acoustic signal and a non-linearly
#' processed version of the signal to partially restore very weak f0 components.
#' Intelligent peak picking to select multiple f0 candidates and assign merit
#' factors; and, incorporation of highly robust pitch contours obtained from
#' smoothed versions of low frequency portions of spectrograms. Dynamic
#' programming is used to find the “best” pitch track among all the candidates,
#' using both local and transition costs.
#'
#'
#'
#' @inheritParams rapt
#' @param windowSize length of each analysis frame (default: 35 ms)
#' @param tda_frame_length The frame length employed in the time domain analysis
#'   (defaults to the same as windowSize 35 ms).
#' @param fft_length FFT length (default: 8192 samples)
#' @param bp_forder order of band-pass filter (default: 150)
#' @param bp_low low frequency of filter passband (default: 50 Hz)
#' @param bp_high high frequency of filter passband (default: 1500 Hz)
#' @param nlfer_thresh1 NLFER (Normalized Low Frequency Energy Ratio) boundary
#'   for voiced/unvoiced decisions (default: 0.75)
#' @param nlfer_thresh2 threshold for NLFER definitely unvoiced (default: 0.1)
#' @param shc_numharms number of harmonics in SHC (Spectral Harmonics
#'   Correlation) calculation (default: 3)
#' @param shc_window SHC window length (default: 40 Hz)
#' @param shc_maxpeaks maximum number of SHC peaks to be found (default: 4)
#' @param shc_pwidth window width in SHC peak picking (default: 50 Hz)
#' @param shc_thresh1 threshold 1 for SHC peak picking (default: 5)
#' @param shc_thresh2 threshold 2 for SHC peak picking (default: 1.25)
#' @param f0_double pitch doubling decision threshold (default: 150 Hz)
#' @param f0_half pitch halving decision threshold (default: 150 Hz)
#' @param dp5_k1 weight used in dynamic program (default: 11)
#' @param dec_factor factor for signal resampling (default: 1)
#' @param nccf_thresh1 threshold for considering a peak in NCCF (Normalized
#'   Cross Correlation Function) (default: 0.3)
#' @param nccf_thresh2 threshold for terminating search in NCCF (default: 0.9)
#' @param nccf_maxcands maximum number of candidates found (default: 3)
#' @param nccf_pwidth window width in NCCF peak picking (default: 5)
#' @param merit_boost boost merit (default. 0.20)
#' @param merit_pivot merit assigned to unvoiced candidates in definitely
#'   unvoiced frames (default: 0.99)
#' @param merit_extra merit assigned to extra candidates in reducing pitch
#'   doubling/halving errors (default: 0.4)
#' @param median_value order of medial filter (default: 7)
#' @param dp_w1 DP (Dynamic Programming) weight factor for voiced-voiced
#'   transitions (default: 0.15)
#' @param dp_w2 DP weight factor for voiced-unvoiced or unvoiced-voiced
#'   transitions (default: 0.5)
#' @param dp_w3 DP weight factor of unvoiced-unvoiced transitions (default: 0.1)
#' @param dp_w4 Weight factor for local costs (default: 0.9)
#'
#' @return An SSFF track object containing two tracks ("f0" and "voiced") which
#'   contains the computed pitch values, and a binary (0 or 1) indication of
#'   whether the frame was considered "voiced" (1)  or not (0). The tracks are
#'   either returned (toFile == FALSE) or stored on disk.
#' @references \insertAllCited{}
#'
#'   
yaapt<- function(listOfFiles,
                   beginTime=0,
                   endTime=0,
                   windowShift=5,
                 windowSize=35,
                   minF=70, 
                   maxF=200,
                   tda_frame_length = 35,
                   fft_length = 8192,
                   bp_forder=  150,
                   bp_low= 50 ,
                   bp_high=  1500 ,
                   nlfer_thresh1= 0.75,
                   nlfer_thresh2=  0.1,
                   shc_numharms= 3,
                   shc_window= 40 ,
                   shc_maxpeaks=  4,
                   shc_pwidth=  50,
                   shc_thresh1=  5,
                   shc_thresh2=  1.25,
                   f0_double= 150, 
                   f0_half= 150, 
                   dp5_k1= 11,
                   dec_factor=  1,
                   nccf_thresh1=  0.3,
                   nccf_thresh2= 0.9,
                   nccf_maxcands=  3,
                   nccf_pwidth=  5,
                   merit_boost=  0.20,
                   merit_pivot=  0.99,
                   merit_extra=  0.4,
                   median_value=  7,
                   dp_w1=  0.15,
                   dp_w2=  0.5,
                   dp_w3= 0.1,
                   dp_w4= 0.9,
                   explicitExt="yf0",
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
    py$windowShift <- reticulate::r_to_py(as.integer(windowShift))
    py$windowSize <- reticulate::r_to_py(as.integer(windowSize))
    py$maxF <- reticulate::r_to_py(as.integer(maxF))
    py$minF <- reticulate::r_to_py(as.integer(minF))
    py$beginTime <- reticulate::r_to_py(as.double(beginTime))
    py$endTime <- reticulate::r_to_py(as.double(endTime))
    py$tda_frame_length <- reticulate::r_to_py(as.integer(tda_frame_length))
    py$fft_length <- reticulate::r_to_py(as.integer(fft_length ))
    py$bp_forder <- reticulate::r_to_py(as.integer(bp_forder))
    py$bp_low <- reticulate::r_to_py(as.integer(bp_low))
    py$bp_high <- reticulate::r_to_py(as.integer(bp_high))
    py$nlfer_thresh1 <- reticulate::r_to_py(as.double(nlfer_thresh1))
    py$nlfer_thresh2 <- reticulate::r_to_py(as.double(nlfer_thresh2))
    py$shc_numharms <- reticulate::r_to_py(as.integer(shc_numharms))
    py$shc_window <- reticulate::r_to_py(as.integer(shc_window))
    py$shc_maxpeaks <- reticulate::r_to_py(as.integer(shc_maxpeaks))
    py$shc_pwidth <- reticulate::r_to_py(as.integer(shc_pwidth))
    py$shc_thresh1 <- reticulate::r_to_py(as.integer(shc_thresh1))
    py$shc_thresh2 <- reticulate::r_to_py(as.double(shc_thresh2))
    py$f0_double <- reticulate::r_to_py(as.integer(f0_double))
    py$f0_half <- reticulate::r_to_py(as.integer(f0_half))
    py$dp5_k1 <- reticulate::r_to_py(as.integer(dp5_k1))
    py$dec_factor <- reticulate::r_to_py(as.integer(dec_factor))
    py$nccf_thresh1 <- reticulate::r_to_py(as.double(nccf_thresh1))
    py$nccf_thresh2 <- reticulate::r_to_py(as.double(nccf_thresh2))
    py$nccf_maxcands <- reticulate::r_to_py(as.integer(nccf_maxcands))
    py$nccf_pwidth <- reticulate::r_to_py(as.integer(nccf_pwidth))
    py$merit_boost <- reticulate::r_to_py(as.double(merit_boost))
    py$merit_pivot <- reticulate::r_to_py(as.double(merit_pivot))
    py$merit_extra <- reticulate::r_to_py(as.double(merit_extra))
    py$median_value <- reticulate::r_to_py(as.integer(median_value))
    py$dp_w1 <- reticulate::r_to_py(as.double(dp_w1))
    py$dp_w2 <- reticulate::r_to_py(as.double(dp_w2))
    py$dp_w3 <- reticulate::r_to_py(as.double(dp_w3))
    py$dp_w4 <- reticulate::r_to_py(as.double(dp_w4))
    
    reticulate::py_run_string("import amfm_decompy.pYAAPT as pYAAPT\
import gc\
import amfm_decompy.basic_tools as basic\
import math as m\
\
signal = basic.SignalObj(soundFile)\
fs = signal.fs\
if endTime > 0.0 or beginTime > 0.0:\
	startSample = m.floor(beginTime * signal.fs)\
	endSample = m.ceil(endTime * signal.fs)\
	subsignal = basic.SignalObj(signal.data[startSample:endSample],signal.fs)\
else:\
	subsignal = signal\
\
pitch = pYAAPT.yaapt(subsignal, **{'f0_min' : minF, \
	'tda_frame_length' : tda_frame_length,\
	'f0_max' : maxF, \
	'frame_length' : windowSize, \
	'frame_space' : windowShift,\
	'tda_frame_length' : tda_frame_length,\
	'fft_length' : fft_length,\
	'bp_forder': bp_forder,\
	'bp_low' : bp_low,\
	'bp_high' : bp_high,\
	'nlfer_thresh1' : nlfer_thresh1,\
	'nlfer_thresh2' : nlfer_thresh2,\
	'shc_numharms' : shc_numharms,\
	'shc_window' : shc_window,\
	'shc_maxpeaks' : shc_maxpeaks,\
	'shc_pwidth' : shc_pwidth,\
	'shc_thresh1' : shc_thresh1,\
	'shc_thresh2' : shc_thresh2,\
	'f0_double' : f0_double,\
	'f0_half' : f0_half,\
	'dp5_k1' : dp5_k1,\
	'dec_factor' : dec_factor,\
	'nccf_thresh1' : nccf_thresh1,\
	'nccf_thresh2' : nccf_thresh2,\
	'nccf_maxcands' : nccf_maxcands,\
	'nccf_pwidth' : nccf_pwidth,\
	'merit_boost' : merit_boost,\
	'merit_pivot' : merit_pivot,\
	'merit_extra'  : merit_extra,\
	'median_value' : median_value,\
	'dp_w1' : dp_w1,\
	'dp_w2' : dp_w2,\
	'dp_w3' : dp_w3,\
	'dp_w4' : dp_w4 })\
\
f0 =  pitch.samp_values\
vuv =   pitch.vuv\
del signal\
gc.collect()")
    
    inTable <- data.frame( "f0" = as.integer(py$f0),
                           "voiced"=ifelse(py$vuv,1,0))
    #return(inTable)
    startTime = windowSize/2
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16","INT16")
    #Use the time separation between second and pitch measurement time stamps to compute a sample frequency.
    
    sampleRate <-  1/ windowShift * 1000
    attr(outDataObj, "sampleRate") <- sampleRate
    
    attr(outDataObj, "origFreq") <-  as.numeric(py$fs) 
    #startTime <- 1/sampleRate
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
    class(outDataObj) = "AsspDataObj"
    
    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    # Pitch track
    f0Table <- inTable %>%
      dplyr::select(f0) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    nof0Values <- nrow(f0Table)
    names(f0Table) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "f0", as.matrix(f0Table[,1]), "INT16")
    
    # Voiced voiceless track
    voicedTable <- inTable %>%
      dplyr::select(voiced) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    noVoicedValues <- nrow(voicedTable)
    names(voicedTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "voiced", as.matrix(voicedTable[,1]), "INT16")
    
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
      
      # prepend values
      outDataObj$f0 = rbind(missing_f0_vals, outDataObj$f0)   
      outDataObj$voiced = rbind(missing_voiced_vals, outDataObj$voiced)   
      
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

attr(yaapt,"ext") <-  c("yf0") 
attr(yaapt,"tracks") <-  c("f0","voiced")
attr(yaapt,"outputType") <-  c("SSFF")


#' Estimates aperiodicity of a speech signal
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
    
    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    
    
    noAperiodicitiesValues <- nrow(aperiodicityTable)
    names(aperiodicityTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "aperiod", as.matrix(aperiodicityTable), "REAL32")
    
    
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



attr(aperiodicities,"ext") <-  c("wap") 
attr(aperiodicities,"tracks") <-  c("aperiod")
attr(aperiodicities,"outputType") <-  c("SSFF")


#' (En)coded spectral envelope estimation
#'
#' The Spectral envelope is estimated using function for assessing band-aperiodicities
#' using the CheapTrick algorithm
#' \insertCite{Morise:2015ia}{superassp} implemented
#' in the WORLD vocoder
#' \insertCite{MORISE.2016.10.1587/transinf.2015edp7457}{superassp}. The the [harvest]
#' pitch algorithm \insertCite{Morise.2017.10.21437/interspeech.2017-68}{superassp} is used to calculate
#' the periodic component.
#'
#'
#' @inheritParams harvest
#' @param dimensions Number of dimensions of coded spectral envelope
#'
#' @return An SSFF track object containing two tracks (f0 and corr) that are
#'   either returned (toFile == FALSE) or stored on disk.
#' @references \insertAllCited{}
#' 
#' @export
#'
#'   
seenc<- function(listOfFiles,
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

f0, t = pw.harvest(x,\
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
    
    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    
    
    noSeencValues <- nrow(seencTable)
    names(seencTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "seenc", as.matrix(seencTable), "REAL32")
    
    
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



attr(seenc,"ext") <-  c("sec") 
attr(seenc,"tracks") <-  c("seenc")
attr(seenc,"outputType") <-  c("SSFF")

## FOR INTERACTIVE TESTING
#library(superassp)
#library(reticulate)
#library(dplyr)
#f <- "/Users/frkkan96/Desktop/a1.wav"
#swipe(f,beginTime=0,endTime=1,toFile=FALSE) -> outportion
#swipe(f,beginTime=0,endTime=0,toFile=FALSE) -> outsa
#reaper(f,beginTime=0,endTime=1,toFile=FALSE) -> outportion
#excite(f,beginTime=0,endTime=0,toFile=FALSE) -> outex
#reaper_pm(f,beginTime=0,endTime=0,toFile=FALSE) -> outpm
#kaldi_pitch(f,beginTime=0,endTime=0,toFile=FALSE) -> outkaldi
#kaldi_pitch(f,beginTime=0,endTime=0,toFile=TRUE)
#dio(f,beginTime=0,endTime=0,toFile=FALSE) -> out
#dio(f,beginTime=0,endTime=0,toFile=TRUE)
#kaldi_pitch(f,beginTime=0,endTime=0,toFile=FALSE) -> out
#harvest(f,beginTime=0,endTime=0,toFile=FALSE) -> out1
#seenc(f,beginTime=0,endTime=0,toFile=FALSE) -> out2

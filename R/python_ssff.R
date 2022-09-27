

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
    origSoundFile <- fileBeginEnd[i, "listOfFiles"]
    
    beginTime <- fileBeginEnd[i, "beginTime"]
    endTime <- fileBeginEnd[i, "endTime"]
    
    py$soundFile <- reticulate::r_to_py(origSoundFile)
    py$ws <- reticulate::r_to_py(windowShift)
    py$fMax <- reticulate::r_to_py(maxF)
    py$fMin <- reticulate::r_to_py(minF)
    py$bt <- reticulate::r_to_py(beginTime)
    py$et <- reticulate::r_to_py(endTime)
    
    
    reticulate::py_run_string("import numpy as np\
import pysptk as sp\
import librosa as lr\
if et > 0:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt, duration = et - bt)\
else:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt)\

f0_swipe = sp.swipe(x.astype(np.float64), fs=fs, hopsize=ws / 1000 * fs, min=fMin, max=fMax, otype=\"f0\")\
pitch_swipe = sp.swipe(x.astype(np.float64), fs=fs, hopsize=ws / 1000 * fs, min=fMin, max=fMax, otype=\"pitch\")")
    
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
    origSoundFile <- fileBeginEnd[i, "listOfFiles"]
    
    beginTime <- fileBeginEnd[i, "beginTime"]
    endTime <- fileBeginEnd[i, "endTime"]
    
    
    py$soundFile <- reticulate::r_to_py(origSoundFile)
    py$ws <- reticulate::r_to_py(windowShift)
    py$fMax <- reticulate::r_to_py(maxF)
    py$fMin <- reticulate::r_to_py(minF)
    py$bt <- reticulate::r_to_py(beginTime)
    py$et <- reticulate::r_to_py(endTime)
    
    
    reticulate::py_run_string("import numpy as np\
import pysptk as sp\
import librosa as lr\
if et > 0:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt, duration = et - bt)\
else:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt)\

f0_rapt = sp.rapt(x.astype(np.float64), fs=fs, hopsize=ws / 1000 * fs, min=fMin, max=fMax, otype=\"f0\")\
pitch_rapt = sp.rapt(x.astype(np.float64), fs=fs, hopsize=ws / 1000 * fs, min=fMin, max=fMax, otype=\"pitch\")")
    
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
#' @references \insertAllCited{}
#'
#' @param unvoiced_cost Set the cost for unvoiced segments. Default is 0.9, the higher the value the more f0 estimates in noise.
#' @param high.pass Perform high-pass filtering to remove DC and low-frequency noise?
#' @param hilbert.transform Remove phase distortion using Hilbert transform?
#' @inheritParams swipe
#'
#'
#' @return An SSFF track object containing two tracks (f0 and corr) that are
#'   either returned (toFile == FALSE) or stored on disk.
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
    origSoundFile <- fileBeginEnd[i, "listOfFiles"]
    
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
import pysptk as sp\
import librosa as lr\
import pyreaper\
if et > 0:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt, duration = et - bt)\
else:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt)\

raw_x = x * 2**15\
int_x = raw_x.astype(np.int16)\

pm_times, pm, f0_times, f0, corr = pyreaper.reaper(x=int_x, fs=fs, minf0 = fMin, maxf0 = fMax, do_high_pass=hp, do_hilbert_transform= ht,  frame_period=ws, inter_pulse=ws, unvoiced_cost =uc)" -> script
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
    origSoundFile <- fileBeginEnd[i, "listOfFiles"]
    
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
import pysptk as sp\
import librosa as lr\
import pyreaper\
if et > 0:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt, duration = et - bt)\
else:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt)\

raw_x = x * 2**15\
int_x = raw_x.astype(np.int16)\

pm_times, pm, f0_times, f0, corr = pyreaper.reaper(x=int_x, fs=fs, minf0 = fMin, maxf0 = fMax, do_high_pass=hp, do_hilbert_transform= ht,  frame_period=ws, inter_pulse=ws, unvoiced_cost =uc)" -> script
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
    outDataObj = wrassp::addTrack(outDataObj, "f0", as.matrix(pmTable[,1]), "INT16")
    
    
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
    origSoundFile <- fileBeginEnd[i, "listOfFiles"]
    
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
    # import numpy as np
    # import pysptk as sp
    # #from scipy.io import wavfile
    # import librosa as lr
    # x, fs = lr.load("/Users/frkkan96/Desktop/a1.wav",dtype=np.float64, offset= 0.5, duration = 15.0)
    # 
    # #fs, x = wavfile.read(pysptk.util.example_audio_file())
    # #x, f = lr.lood(pysptk.util.example_audio_file())
    # pitch_swipe = sp.swipe(x, fs=fs, hopsize=80, min=60, max=200, otype="pitch")
    # 
    # exct = sp.excite(pitch_swipe, hopsize=80)
    # 

        
    "import numpy as np\
import pysptk as sp\
import librosa as lr\
import pyreaper\
if et > 0:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt, duration = et - bt)\
else:\
  x, fs = lr.load(soundFile,dtype=np.float64, offset= bt)\

hs = ws * fs\
pitch_swipe = sp.swipe(x.astype(np.float64), fs=fs, hopsize=hs, min=fMin, max=fMax, otype=\"pitch\",threshold=vt)\
exct = sp.excite(pitch_swipe, hopsize=hs)" -> script
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



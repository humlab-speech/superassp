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
trk_yaapt <- function(listOfFiles,
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
pitch = pYAAPT.trk_yaapt(subsignal, **{'f0_min' : minF, \
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
    
    # Voiced voiceless track
    voicedTable <- inTable %>%
      dplyr::select(voiced) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    noVoicedValues <- nrow(voicedTable)
    names(voicedTable) <- NULL
    
    outDataObj = addTrack(outDataObj, "voiced", as.matrix(voicedTable[,1]), "INT16")
    
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

attr(trk_yaapt,"ext") <-  c("yf0") 
attr(trk_yaapt,"tracks") <-  c("f0","voiced")
attr(trk_yaapt,"outputType") <-  c("SSFF")





#' Compute pitch and periodicity using the CREPE pitch tracker
#'
#' The CREPE \insertCite{Kim.2018.10.1109/icassp.2018.8461329}{superassp} applies a deep convolutional neural network directly on the time-domain waveform to find
#' the fundamental frequency in a speech signal. Two versions of the models have been trained, one smaller yielding quicker results, and the full model which can be considerably
#' more computationally intensive to apply.
#'
#' This function uses the \code{av} package to load audio files, supporting all media formats
#' (WAV, MP3, MP4, video files, etc.).
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
#' @inheritParams trk_formantp
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
trk_crepe <- function(listOfFiles,
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

    # Initialize Python environment
    py <- reticulate::import_main()

    # Check if torchcrepe is available
    if (!reticulate::py_module_available("torchcrepe")) {
      stop(
        "trk_crepe() requires the torchcrepe Python module.\n\n",
        "Install with: pip install torchcrepe\n\n",
        "Note: torchcrepe requires PyTorch and CUDA (for GPU acceleration).\n",
        "For CPU-only usage, you may need to install PyTorch CPU version first:\n",
        "  pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu\n",
        "  pip install torchcrepe",
        call. = FALSE
      )
    }

    # Load audio using av package (supports all media formats)
    audio_data <- av::read_audio_bin(
      audio = origSoundFile,
      start_time = if (beginTime > 0) beginTime else NULL,
      end_time = if (endTime > 0) endTime else NULL,
      channels = 1
    )

    # Get sample rate
    sr <- attr(audio_data, "sample_rate")

    # Convert to float32 for Python/PyTorch
    audio_float <- as.numeric(audio_data) / 2147483647.0  # INT32_MAX

    # Import torch and create tensor
    torch <- reticulate::import("torch", convert = FALSE)
    audio_tensor <- torch$from_numpy(
      reticulate::np_array(audio_float, dtype = "float32")
    )

    # Pass parameters to Python
    py$audio <- audio_tensor
    py$sr <- reticulate::r_to_py(as.integer(sr))
    py$windowShift <- reticulate::r_to_py(windowShift)
    py$windowSize <- reticulate::r_to_py(windowSize)
    py$model <- reticulate::r_to_py("tiny")
    py$fmax <- reticulate::r_to_py(maxF)
    py$fmin <- reticulate::r_to_py(minF)
    py$silence_threshold <- reticulate::r_to_py(silence.threshold)
    py$voicing_threshold <- reticulate::r_to_py(voicing.threshold)

    reticulate::py_run_string("import torchcrepe\
import gc\
import math\
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
    
    # Auto-correlation track
    periodicityTable <- inTable %>%
      dplyr::select(periodicity) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    noPeriodicityValues <- nrow(periodicityTable)
    names(periodicityTable) <- NULL
    outDataObj = addTrack(outDataObj, "periodicity", as.matrix(periodicityTable[,1]), "REAL32")
    
    
    
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
    
    assertthat::assert_that(is.AsspDataObj(outDataObj),
                            msg = "The AsspDataObj created by the crepe function is invalid.")
    
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



attr(trk_crepe,"ext") <-  c("crp") 
attr(trk_crepe,"tracks") <-  c("f0","periodicity")
attr(trk_crepe,"outputType") <-  c("SSFF")



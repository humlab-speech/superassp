#' Pitch tracking using the torch pitch tracker
#'
#' This function estimates pitch by normalized cross-correlation function (NCCF) and
#' median smoothing, as implemented in the torchaudio
#' \insertCite{yang2021torchaudio}{superassp} library. The exact algorithm is undisclosed by the implementing library
#' but approach likely builds on earlier implementations that use NCCFs \insertCite{talkin1995robust,Kasi.2002.10.1109/icassp.2002.5743729}{superassp} 
#' including the \link[superassp::rapt]{RAPT} algorithm. 
#'
#' @inheritParams swipe
#'
#' @return
#'  An SSFF track object containing two tracks (f0 and pitch) that are either returned (toFile == FALSE) or stored on disk.
#' 
#' @seealso rapt
#' @references \insertAllCited{}
#' 
torch_pitch <- function(listOfFiles,
                  beginTime=0,
                  endTime=0,
                  windowShift=10,
                  windowSize=30,
                  minF=70, 
                  maxF=200, 
                  explicitExt="tpi",
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


    audio = transform_to_tensor(audiofile_loader(filepath=origSoundFile,
                                                # offset=beginTime,
                                                # duration=(endTime - beginTime), #A duration of 0 seems to be interpreted as the complete file
                                                 unit="time"))
    waveform <- audio[[1]]
    sample_rate <- audio[[2]]
    
    pitch <- functional_detect_pitch_frequency(waveform,
                                             sample_rate = sample_rate,
                                             frame_time = windowShift /1000,
                                             win_length = windowSize,
                                             freq_low=minF,
                                             freq_high=maxF) # Expects seconds
    
    # nfcc <- functional__compute_nccf(waveform,
    #                                sample_rate = sample_rate,
    #                                frame_time = windowShift/1000,
    #                                freq_low = minF)
    # return(c(pitch,nfcc))
    inTable <- data.frame( "f0" = as.integer(pitch)#,
                           #"nfcc"=as.numeric(nfcc)
                           )
    
    startTime = windowShift
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16")
    #Use the time separation between second and pitch measurement time stamps to compute a sample frequency.
    
    sampleRate <-  1/ windowShift * 1000
    attr(outDataObj, "sampleRate") <- sampleRate
    
    attr(outDataObj, "origFreq") <-  as.numeric(sample_rate) 
    startTime <- 1/sampleRate
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
    class(outDataObj) = "AsspDataObj"
    
    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    return(inTable)
    # Cross-correlation track
    f0Table <- inTable %>%
      dplyr::select(all_of(f0)) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    
    nof0Values <- nrow(f0Table)
    names(f0Table) <- NULL
    outDataObj = wrassp::addTrack(outDataObj, "f0", as.matrix(f0Table[,1]), "INT16")
    
    # # Auto-correlation track
    # nfccTable <- inTable %>%
    #   dplyr::select(nfcc) %>%
    #   replace(is.na(.), 0) %>%
    #   dplyr::mutate(
    #     dplyr::across(
    #       tidyselect::everything(),as.integer))
    # 
    # nfccTable <- nrow(nfccTable)
    # names(pitchTable) <- NULL
    # outDataObj = wrassp::addTrack(outDataObj, "nfcc", as.matrix(nfccTable[,1]), "INT16")
    # 
    # 
    
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
      # missing_nfcc_vals = matrix(0,
      #                             nrow = nr_of_missing_samples,
      #                             ncol = ncol(outDataObj$nfcc))
      
      # prepend values
      outDataObj$f0 = rbind(missing_f0_vals, outDataObj$f0)
      # outDataObj$nfcc = rbind(missing_nfcc_vals, outDataObj$nfcc)
      # 
      
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

attr(torch_pitch,"ext") <-  c("tpi") 
attr(torch_pitch,"tracks") <-  c("f0")
attr(torch_pitch,"outputType") <-  c("SSFF")


## FOR INTERACTIVE TESTING
#library(superassp)
#library(torch)
#library(torchaudio)
#library(dplyr)
#f <- "/Users/frkkan96/Desktop/a1.wav"
    
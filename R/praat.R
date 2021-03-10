add_trackDefinition <- function(
  emuDBhandle,
  name,
  columnName = NULL,
  fileExtension = NULL,
  onTheFlyFunctionName = NULL,
  onTheFlyParams = NULL,
  onTheFlyOptLogFilePath = NULL,
  verbose = TRUE,
  interactive = TRUE){
  
  # If the function extists in wrassp, just call that function
  if(!is.null(wrassp::wrasspOutputInfos[[onTheFlyFunctionName]])){
    emuR::add_ssffTrackDefinition(emuDBhandle=emuDBhandle,
                                  name=name,
                                  columnName = columnName,
                                  fileExtension = fileExtension,
                                  onTheFlyFunctionName = onTheFlyFunctionName,
                                  onTheFlyParams = onTheFlyParams,
                                  onTheFlyOptLogFilePath = onTheFlyOptLogFilePath,
                                  verbose = verbose,
                                  interactive = interactive)
    
  }else{
    
    # Check that the function extists 
    if(exists(onTheFlyFunctionName) & is.function(get(onTheFlyFunctionName))){
      
      fun <- get(onTheFlyFunctionName)
      #Check that the function has been prepared for use with this function by 
      # giving it the the required additional attributes "ext" and "tracks"
      if(!is.null(attr(fun,"ext")) & !is.null(attr(fun,"tracks")) ){
        #Set the default file extension to the one set as an attribute, if missing in the arguments
        ext <- ifelse(!is.null(fileExtension),fileExtension,attr(fun,"ext"))
        if(!columnName %in% attr(fun,"tracks") ) stop("The track ",columnName, " is not a defined output track name of the function ",onTheFlyFunctionName)
        columnName <- ifelse(!is.null(columnName),columnName,attr(fun,"tracks")[[1]])
        
      }else{
        stop("The function ",onTheFlyFunctionName," is not defined correctly. Please provide it with the attributes \"ext\" and \"tracks\".\n See ?attr for details, as well as attributes(praat_formant_burg) for an example." )
      }
      dbConfig = emuR:::load_DBconfig(emuDBhandle)
      funcFormals = formals(onTheFlyFunctionName)
      funcFormals[names(onTheFlyParams)] = onTheFlyParams
      funcFormals$optLogFilePath = onTheFlyOptLogFilePath
      fp = emuR::list_files(emuDBhandle, dbConfig$mediafileExtension)
      funcFormals$listOfFiles = paste(emuDBhandle$basePath, paste0(fp$session, emuR:::session.suffix), paste0(fp$bundle, emuR:::bundle.dir.suffix), fp$file, sep = .Platform$file.sep)
      funcFormals$explicitExt = fileExtension
      funcFormals$verbose = verbose
      do.call(onTheFlyFunctionName, funcFormals)
      #add the definition
      emuR::add_ssffTrackDefinition(emuDBhandle,name=name,columnName = columnName,fileExtension = ext)
    }else{
      
      stop("Could not find a definition of the function ",onTheFlyFunctionName,"." )
    }
  }
}


#' A simple check of a presence of a Praat executable
#' 
#' 
#'
#' @param praat_path A character string containing the path to the executable that the function was able to find (or the executable that the function was able to verify the existance of), or NULL if no Praat executable was found and verified.
#'
#' @return A boolean indicating whether the Praat executable could be found or not.
#' @export
#'
have_praat <- function(praat_path=NULL){
  return(ifelse(is.null(get_praat(praat_path)) ,FALSE,TRUE))
}



#' Utility function for getting the full path of the Praat executable.
#' 
#' This function checks the system variables and deduces where Praat is installed. On the OSX platform (Darwin) the Praat app is assumed to exist in the Applications folder, 
#' and the actual binary inside of the application package is then used. If not OSX, then the function will search the default search paths for executables set up in the OS. 
#' If an explicit path is given, then function will just check whether the executable is actualy present there.
#'
#' @param praat_path A character string containing the path to the executable that the function was able to find (or the executable that the function was able to verify the existance of), or NULL if no Praat executable was found and verified.
#'
#' @return A character string containing the path to the executable that the function was able to find (or the executable that the function was able to verify the existance of), or NULL if no Praat executable was found and verified.
#' @export
#'
get_praat <- function(praat_path=NULL){
  sysname <- as.vector(Sys.info()["sysname"])
  
  if(is.null(praat_path)){
    praat_path <- switch(sysname,
                         Darwin = "/Applications/Praat.app/Contents/MacOS/Praat",
                         Sys.which("praat")
    )    
  }
  if(! file.exists(praat_path)){
    praat_path <- NULL
    warning("The function get_praat could not find your praat binary. Please provide a full path directly.")
  }
  
  return(praat_path)
  
}


#' Use Praat to compute a formant track using the burg method.
#' 
#' This function asks Praat to compute a formant track using the burg method, and the result is converted in an SSFF file.
#' This function should be a drop-in replacement for the \code{\link[wrassp]{forest}} function, but with some additional
#' arguments. If the function cannot find the Praat binary automatically, you have to give an explicit path (e.g. "/Applications/Praat.app/Contents/MacOS/Praat" if you placed Praat in the Applications folder on your Mac). 
#' You can check whether you need to supply a explicit path using the \code{\link{get_praat}} or 
#' \code{\link{have_praat}} functions.
#' 
#'
#' @param listOfFiles a vector of wav file paths to be processed by function.
#' @param beginTime the time where processing should end (in s) The default is 0 (zero) which means that the computation of formants will start at the start of the sound file.
#' @param endTime the time where processing should end (in s) The default is 0 (zero) which means that formants will be computed up to the end of the file.
#' @param windowShift the analysis window shift length (in ms).
#' @param numFormants the number of formants that the analysis should try to find 
#' @param maxhzformant praat will try to find formants only up to this frequency in the spectrum.
#' @param windowSize the analysis window length (in ms).
#' @param preemphasis the frequency from which a preemphasis will be applied..
#' @param window the analysis window function used when extracting part of a sound file for analysis. De faults to "Hanning".
#' @param relativeWidth the relative width of the windowing function used.
#' @param toFile write the output to a file? The file will be written in  `outputDirectory`, if defined, or in the same directory as the soundfile. 
#' @param explicitExt the file extension that should be used.
#' @param outputDirectory set an explicit directory for where the signal file will be written. If not defined, the file will be written to the same directory as the sound file.
#' @param verbose Not implemented. Only included here for compatibility.  
#' @param praat_path give an explicit path for Praat. If the praat 
#'
#' @return a list of 
#' @export
#'
#' 
#' 
praat_formant_burg <- function(listOfFiles,beginTime=0,endTime=0,windowShift=0.0,numFormants=4.0,maxhzformant=5500.0,windowSize=0.025,preemphasis=50.0,window="hanning",relativeWidth=1.0,toFile=TRUE,explicitExt="fms",outputDirectory=NULL,verbose=FALSE,praat_path=NULL){


  
  if(! has_praat(praat_path)){
    stop("Could not find praat. Please specify a full path.")
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
  

  
  praat_script <- ifelse(PRAAT_DEVEL== TRUE,
                     file.path("inst","praat","formant_burg.praat"),
                     file.path(system.file(package = "superassp",mustWork = TRUE),"praat","formant_burg.praat")
                     )
  
  formant_burg <- tjm.praat::wrap_praat_script(praat_location = get_praat(),
                                    script_code_to_run = readLines(praat_script)
                                    ,return="last-argument")
  
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
    
    formantTabFile <- tempfile(fileext = ".csv")

    #Required for preventing errors in the handoff of file names containing spaces and () characters
    # to Praat
    soundFile <- tempfile(fileext = ".wav")
    R.utils::createLink(soundFile,origSoundFile)
    #Alternative route - much slower
    #file.copy(origSoundFile,soundFile)
    

    outFormantTabFile <- formant_burg(soundFile,
                            beginTime,
                            endTime,
                            windowShift,
                            numFormants,
                            maxhzformant,
                            windowSize,
                            preemphasis,
                            window,
                            relativeWidth,
                            formantTabFile)
    
    inTable <- read.csv(file=outFormantTabFile
                               ,header=TRUE
                               ,na.strings =c("--undefined--","NA"),
                        sep = ",")
      

    
    # We need the sound file to extract some information
    origSound <- wrassp::read.AsspDataObj(soundFile)

    starTime = inTable[1,"time.s."]
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16", "INT16")
    #Use the time separation between second and first formant measurement time stamps to compute a sample frequency.
    sampleRate <-  as.numeric(1 / (inTable[2,"time.s."] - inTable[1,"time.s."]))
    attr(outDataObj, "sampleRate") <- sampleRate

    attr(outDataObj, "origFreq") <-  as.numeric(attr(origSound, "sampleRate"))
    startTime <- as.numeric(inTable[1,"time.s."])
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
    class(outDataObj) = "AsspDataObj"

    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary

    fmTable <- inTable %>%
      dplyr::select(tidyselect::starts_with("F",ignore.case = FALSE)) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))

    noFormantsValues <- nrow(fmTable)
    noFormants <- ncol(fmTable)
    
    names(fmTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "fm", as.matrix(fmTable), "INT16")

    bwTable <- inTable %>%
      dplyr::select(tidyselect::starts_with("B",ignore.case = FALSE)) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))

    names(bwTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "bw", as.matrix(bwTable), "INT16")


    ## Apply fix from Emu-SDMS manual
    ##https://raw.githubusercontent.com/IPS-LMU/The-EMU-SDMS-Manual/master/R/praatToFormants2AsspDataObj.R

    # add missing values at the start as Praat sometimes
    # has very late start values which causes issues
    # in the SSFF file format as this sets the startRecord
    # depending on the start time of the first sample
    if( startTime > (1/sampleRate) ){

      nr_of_missing_samples = as.integer(floor(startTime / (1/sampleRate)))

      missing_fm_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$fm))


      missing_bw_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$bw))

      # prepend values
      outDataObj$fm = rbind(missing_fm_vals, outDataObj$fm)
      outDataObj$bw = rbind(missing_fm_vals, outDataObj$bw)

      # fix start time
      attr(outDataObj, "startTime") = startTime - nr_of_missing_samples * (1/sampleRate)
    }

    assertthat::assert_that(wrassp::is.AsspDataObj(outDataObj),
                            msg = paste("The AsspDataObj created by the praat_formant_burg function is invalid.\nPlease check the table file '",tabfile,"' for errors.",sep=""))

    ssff_file <- gsub("wav$",explicitExt,origSoundFile)
    if(!is.null(outputDirectory)){
      ssff_file <- file.path(outputDirectory,basename(ssff_file))
    }
    
    attr(outDataObj,"filePath") <- as.character(ssff_file)
    if(toFile){
      wrassp::write.AsspDataObj(dobj=outDataObj,file=ssff_file)
      #Here we can be sure that the list is a valid SSFF object, so the
      # so we add TRUE to the out vector
      outListOfFiles <- c(outListOfFiles,TRUE)
    }

  }

  if(toFile){
    return(length(outListOfFiles))
  }else{
    return(outDataObj)
  }
    
}

attr(praat_formant_burg,"ext") <-  c("fms") 
attr(praat_formant_burg,"tracks") <-  c("fm", "bw")

# fm
# F1: The frequency of the first formant
# F2: The frequency of the second formant
# F3: The frequency of the third formant
# bw
# B1: The bandwidth of the first formant
# B2: The bandwidth of the second formant
# B3: The bandwidth of the third formant
# 
# uncorrected_harmonics
# H1u: The uncorrected amplitude of the first harmonic
# H2u: The uncorrected amplitude of the second harmonic
# H4u: The uncorrected amplitude of the fourth harmonic
# 
# uncorrected_harmonics_k
# H2Ku: The uncorrected amplitude of the harmonic closest to 2000Hz
# H5Ku: The uncorrected amplitude of the harmonic closest to 5000Hz
# 
# uncorrected_harmonics_formants
# A1u: The uncorrected amplitude of the harmonic closest to the first formant
# A2u: The uncorrected amplitude of the harmonic closest to the second formant
# A3u: The uncorrected amplitude of the harmonic closest to the third formant
# 
# H1H2u: The difference in amplitudes of the uncorrected first and second harmonics
# 
# H2H4u: The difference in amplitudes of the uncorrected first and fourth harmonics
# 
# H1A1u: The difference in the uncorrected amplitudes of the first harmonic and the harmonic closes to the first formant
# 
# H1A2u: The difference in the uncorrected amplitudes of the first harmonic and the harmonic closes to the second formant
# 
# H1A3u: The difference in the uncorrected amplitudes of the first harmonic and the harmonic closes to the third formant
# 
# H2KH5Ku: The difference in the uncorrected amplitudes of the harmonics closest to 2kHz and 5kHz respectivelly
# 
# 
# corrected_harmonics
# H1c: The corrected amplitude of the first harmonic
# H2c: The corrected amplitude of the second harmonic
# H4c: The corrected amplitude of the fourth harmonic
# 
# corrected_harmonics_formants
# A1c: The corrected amplitude of the harmonic closest to the first formant
# A2c: The corrected amplitude of the harmonic closest to the first formant
# A3c: The corrected amplitude of the harmonic closest to the first formant
# 
# 
# H1H2c: The difference in amplitudes of the corrected first and second harmonics
# 
# H2H4c: The difference in amplitudes of the corrected first and second harmonics
# 
# H1A1c: The difference in the corrected amplitudes of the first harmonic and the harmonic closest to the first formant
# 
# H1A2c: The difference in the corrected amplitudes of the first harmonic and the harmonic closest to the first formant
# 
# H1A3c: The difference in the corrected amplitudes of the first harmonic and the harmonic closest to the first formant
# 
# CPP: The cepstral peak prominence (Hillenbrand et al., 1994)
# 
# HNR05: The Harmonic-to-noise ratio as defined by de Krom (1993), measured from 0-500 Hz
# 
# HNR15: The Harmonic-to-noise ratio as defined by de Krom (1993), measured from 0-1500 Hz
# 
# HNR25: The Harmonic-to-noise ratio as defined by de Krom (1993), measured from 0-2500 Hz
# 
# HNR35: The Harmonic-to-noise ratio as defined by de Krom (1993), measured from 0-3500 Hz



praat_sauce <- function(listOfFiles,beginTime=0,endTime=0,channel=1,measure=2,points=5,resample_to_16k=TRUE,pitchTracking=TRUE,formantMeasures=TRUE,spectralMeasures=TRUE,windowLength=0.025,windowPosition=0.5,maxFormantHz=5000,spectrogramWindow=0.005,f0min=50,f0max=300,timeStep=0,maxNumFormants=5,preEmphFrom=50,formantTracking=1,F1ref=500,F2ref=1500,F3ref=2500,useBandwidthFormula=FALSE,toFile=TRUE,explicitExt="psa",outputDirectory=NULL,verbose=FALSE,praat_path=NULL){
  
  if( ! (pitchTracking|formantMeasures|spectralMeasures ) ){
    stop("Calling the praat_sauce function without wanting some acoustic measurements in return makes no sense.\n",
         "Please set either pitchTracking, formantMeasures or spectralMeasures to TRUE.1") 
  }
  
  #Right now, you cannot get spectral measures from praatsauce without computing formants and pitch too
  # so the user will get them too
  if(spectralMeasures){
    formantMeasures = TRUE
    pitchTracking = TRUE
  }
  
  
  
  if(! has_praat(praat_path)){
    stop("Could not find praat. Please specify a full path.")
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
  
  
  
  praat_script <- ifelse(PRAAT_DEVEL== TRUE,
                         file.path("inst","praat","praatsauce.praat"),
                         file.path(system.file(package = "superassp",mustWork = TRUE),"praat","praatsauce.praat")
  )
  
  additional_script_names <- c("formantMeasures.praat","pitchTracking.praat","spectralMeasures.praat","correct_iseli_z.praat","getbw_HawksMiller.praat") 
  
  additional_scripts <- c()
  for(scriptname in additional_script_names){
    additional_scripts <- c(additional_scripts,
                            ifelse(PRAAT_DEVEL== TRUE,
                                   file.path("inst","praat",scriptname),
                                   file.path(system.file(package = "superassp",mustWork = TRUE),"praat",scriptname))
                            )
                            
  }
  
  praatsauce <- tjm.praat::wrap_praat_script(praat_location = get_praat(),
                                               script_code_to_run = readLines(praat_script)
                                               ,return="last-argument")
  #Copy additional files
  copied <- file.copy(additional_scripts,tempdir(),overwrite = TRUE)

  #Check that all files exists before we begin
  filesEx <- file.exists(listOfFiles)
  if(!all(filesEx) ){
    filedNotExists <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ",paste(filedNotExists, collapse = ", "))
  }
  if(!all(copied)){
    stop("Not all required praat script files were copied correctly: The script files ", paste(additional_scripts[!copied],collapse = ",")," are missing.")
  }
  
  #The empty vector of file names that should be returned
  outListOfFiles <- c()
  
  for(i in 1:nrow(fileBeginEnd)){ 
    origSoundFile <- fileBeginEnd[i, "listOfFiles"]
    
    beginTime <- fileBeginEnd[i, "beginTime"]
    endTime <- fileBeginEnd[i, "endTime"]
    
    outputfile <- tempfile(fileext = ".csv")
    
    #Required for preventing errors in the handoff of file names containing spaces and () characters
    # to Praat
    soundFile <- tempfile(fileext = ".wav")
    R.utils::createLink(soundFile,origSoundFile)
    #Alternative route - much slower
    #file.copy(origSoundFile,soundFile)
    
    
    outputfile <- praatsauce(soundFile,
                                      beginTime,
                                      endTime,
                                      channel,
                                      measure,
                                      points,
                                      ifelse(resample_to_16k,1,0),
                                      ifelse(pitchTracking,1,0),
                                      ifelse(formantMeasures,1,0),
                                      ifelse(spectralMeasures,1,0),
                                      windowLength,
                                      windowPosition,
                                      maxFormantHz,
                                      spectrogramWindow,
                                      0, #ifelse(useExistingPitch,1,0), TO IMPLEMENT LATER
                                      f0min,
                                      f0max,
                                      timeStep,
                                      maxNumFormants,
                                      preEmphFrom,
                                      ifelse(formantTracking,1,0),
                                      F1ref,
                                      F2ref,
                                      F3ref,
                                      0, #ifelse(useExistingFormants,1,0), TO IMPLEMENT LATER
                                      ifelse(useBandwidthFormula,1,0),
                                      outputfile)
    
    inTable <- read.csv(file=outputfile
                        ,header=TRUE
                        ,na.strings =c("--undefined--","NA"),
                        sep = ",")
    
    

    #####
    # Create the SSFF object
    #####
    
    # We need the sound file to extract some information
    origSound <- wrassp::read.AsspDataObj(soundFile)
    
    starTime = inTable[1,"t"]
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- rep("INT16",ncol(inTable)-1) #All but the "t" column
    #Use the time separation between second and first formant measurement time stamps to compute a sample frequency.
    sampleRate <-  as.numeric(1 / (inTable[2,"t"] - inTable[1,"t"]))
    attr(outDataObj, "sampleRate") <- sampleRate
    
    attr(outDataObj, "origFreq") <-  as.numeric(attr(origSound, "sampleRate"))
    startTime <- as.numeric(inTable[1,"t"])
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
    class(outDataObj) = "AsspDataObj"
    
    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    if(pitchTracking){
      ###########
      ## f0 values are placed in the track "f0"
      ###########
      # Too elaborate for a single column manipulation, but it is consistent with subsequent processing of 
      # multiple columns
      f0Table <- inTable %>%
        dplyr::select(f0) %>%
        replace(is.na(.), 0) %>%
        dplyr::mutate(f0=as.integer(f0))
      
      names(f0Table) <- NULL
      
      outDataObj = wrassp::addTrack(outDataObj, "f0", as.matrix(f0Table), "INT16")
      
      
    }
    
    ### 
    # This properties are only extracted if the user has specified that formants should be computed,
    # OR spectral measures are to be computed
    ### 
    if(formantMeasures){
        
      ########### Something
      # Formant frequencies are placed in the track "fm"
      ###########
      
      fmTable <- inTable %>%
        dplyr::select(tidyselect::starts_with("F",ignore.case = FALSE)) %>%
        replace(is.na(.), 0) %>%
        dplyr::mutate(
          dplyr::across(
            tidyselect::everything(),as.integer))
      
      #noFormantsValues <- nrow(fmTable)
      #noFormants <- ncol(fmTable)
      
      names(fmTable) <- NULL
      
      outDataObj = wrassp::addTrack(outDataObj, "fm", as.matrix(fmTable), "INT16")
      
      ###########
      ## Formant bandwidths are placed in the track "bw"
      ###########
      
      bwTable <- inTable %>%
        dplyr::select(tidyselect::starts_with("B",ignore.case = FALSE)) %>%
        replace(is.na(.), 0) %>%
        dplyr::mutate(
          dplyr::across(
            tidyselect::everything(),as.integer))
      
      names(bwTable) <- NULL
      
      outDataObj = wrassp::addTrack(outDataObj, "bw", as.matrix(bwTable), "INT16")
    
    }
    
    
    ### 
    # Extraction of the more special voice characteristics
    ###    
    
    
    if(spectralMeasures){
      
      ###########
      ## Amplitudes (uncorrected) of harmonics  are placed in the track "H"
      ###########
      
      harmTable <- inTable %>%
        dplyr::select(H1u:H4u) %>%
        replace(is.na(.), 0) %>%
        dplyr::mutate(
          dplyr::across(
            tidyselect::everything(),as.integer))
      
      names(harmTable) <- NULL
      
      outDataObj = wrassp::addTrack(outDataObj, "H", as.matrix(harmTable), "INT16")

      ###########
      ## Corrected amplitudes  of harmonics  are placed in the track "Hc"
      ###########
      
      harmTable <- inTable %>%
        dplyr::select(H1c:H4c) %>%
        replace(is.na(.), 0) %>%
        dplyr::mutate(
          dplyr::across(
            tidyselect::everything(),as.integer))
      
      names(harmTable) <- NULL
      
      outDataObj = wrassp::addTrack(outDataObj, "Hc", as.matrix(harmTable), "INT16")
      
      
      ###########
      ## The (uncorrected) amplitudes of harmonics  closest to F1-F3 are placed in the track "A"
      ###########
      
      harmTable <- inTable %>%
        dplyr::select(A1u:A3u) %>%
        replace(is.na(.), 0) %>%
        dplyr::mutate(
          dplyr::across(
            tidyselect::everything(),as.integer))
      
      names(harmTable) <- NULL
      
      outDataObj = wrassp::addTrack(outDataObj, "A", as.matrix(harmTable), "INT16")
      ###########
      ## The corrected amplitudes of harmonics closest to F1-F3 are placed in the track "Ac"
      ###########
      
      harmTable <- inTable %>%
        dplyr::select(A1c:A3c) %>%
        replace(is.na(.), 0) %>%
        dplyr::mutate(
          dplyr::across(
            tidyselect::everything(),as.integer))
      
      names(harmTable) <- NULL
      
      outDataObj = wrassp::addTrack(outDataObj, "Ac", as.matrix(harmTable), "INT16")
      
      ###########
      ## The (uncorrected) first and second harmonics closest to 2 and 5k Hz respectively are placed in the columns in "H25K"
      ###########
      harmTable <- inTable %>%
        dplyr::select(H2Ku,H5Ku) %>%
        replace(is.na(.), 0) %>%
        dplyr::mutate(
          dplyr::across(
            tidyselect::everything(),as.integer))
      
      names(harmTable) <- NULL
      
      outDataObj = wrassp::addTrack(outDataObj, "H25K", as.matrix(harmTable), "INT16")

      
      ###########
      ## The  differences between the (uncorrected) amplitudes of the first and second harmonics and the
      ## second and fourth are stored in the "HH" field
      ###########

      harmTable <- inTable %>%
        dplyr::select(H1H2u,H2H4u) %>%
        replace(is.na(.), 0) %>%
        dplyr::mutate(
          dplyr::across(
            tidyselect::everything(),as.integer))
      
      names(harmTable) <- NULL
      
      outDataObj = wrassp::addTrack(outDataObj, "HH", as.matrix(harmTable), "INT16")
      
      ###########
      ## The  differences between the corrected amplitudes of the first and second harmonics and the
      ## second and fourth are stored in the "HHc" field
      ###########
      
      harmTable <- inTable %>%
        dplyr::select(H1H2c,H2H4c) %>%
        replace(is.na(.), 0) %>%
        dplyr::mutate(
          dplyr::across(
            tidyselect::everything(),as.integer))
      
      names(harmTable) <- NULL
      
      outDataObj = wrassp::addTrack(outDataObj, "HHc", as.matrix(harmTable), "INT16")

      ###########
      ## The  differences between the (uncorrected) amplitudes of the first harmonic and the harmonics
      ## closest to F1,F2, and F3 are placed as columns in the "HAd" field
      ###########
      
      harmTable <- inTable %>%
        dplyr::select(H1A1u,H1A2u,H1A3u) %>%
        replace(is.na(.), 0) %>%
        dplyr::mutate(
          dplyr::across(
            tidyselect::everything(),as.integer))
      
      names(harmTable) <- NULL
      
      outDataObj = wrassp::addTrack(outDataObj, "HA", as.matrix(harmTable), "INT16")
      
      ###########
      ## The  differences between the corrected amplitudes of the first harmonic and the harmonics
      ## closest to F1,F2, and F3 are placed as columns in the "HAc" field
      ###########
      
      harmTable <- inTable %>%
        dplyr::select(H1A1c,H1A2c,H1A3c) %>%
        replace(is.na(.), 0) %>%
        dplyr::mutate(
          dplyr::across(
            tidyselect::everything(),as.integer))
      
      names(harmTable) <- NULL
      
      outDataObj = wrassp::addTrack(outDataObj, "HAc", as.matrix(harmTable), "INT16")

      ###########
      ## The cepstral peak prominence is inserted into the "cpp" field
      ###########
      
      cppTable <- inTable %>%
        dplyr::select(CPP) %>%
        replace(is.na(.), 0) %>%
        dplyr::mutate(
          dplyr::across(
            tidyselect::everything(),as.integer))
      
      names(cppTable) <- NULL
      
      outDataObj = wrassp::addTrack(outDataObj, "cpp", as.matrix(cppTable), "INT16")
      
      ###########
      ## The harmonic-to-noise ratio measured from 0 to 500, 1500, 2500 and 3500 Hz respectively 
      ## is inserted into columns of the field "hnr"
      ###########
      
      hnrTable <- inTable %>%
        dplyr::select(HNR05,HNR15,HNR25,HNR35) %>%
        replace(is.na(.), 0) %>%
        dplyr::mutate(
          dplyr::across(
            tidyselect::everything(),as.integer))
      
      names(hnrTable) <- NULL
      
      outDataObj = wrassp::addTrack(outDataObj, "hnr", as.matrix(hnrTable), "INT16")     

    }
    
    ## Apply fix from Emu-SDMS manual
    ##https://raw.githubusercontent.com/IPS-LMU/The-EMU-SDMS-Manual/master/R/praatToFormants2AsspDataObj.R
    
    # add missing values at the start as Praat sometimes
    # has very late start values which causes issues
    # in the SSFF file format as this sets the startRecord
    # depending on the start time of the first sample
    if( startTime > (1/sampleRate) ){
      
      nr_of_missing_samples = as.integer(floor(startTime / (1/sampleRate)))
      
      missing_fm_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$fm))
      
      
      missing_bw_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$bw))
      
      # prepend values
      outDataObj$fm = rbind(missing_fm_vals, outDataObj$fm)
      outDataObj$bw = rbind(missing_fm_vals, outDataObj$bw)
      
      # fix start time
      attr(outDataObj, "startTime") = startTime - nr_of_missing_samples * (1/sampleRate)
    }
    
    assertthat::assert_that(wrassp::is.AsspDataObj(outDataObj),
                            msg = paste("The AsspDataObj created by the praat_sauce function is invalid.\nPlease check the table file '",tabfile,"' for errors.",sep=""))
    
    ssff_file <- gsub("wav$",explicitExt,origSoundFile)
    if(!is.null(outputDirectory)){
      ssff_file <- file.path(outputDirectory,basename(ssff_file))
    }
    
    attr(outDataObj,"filePath") <- as.character(ssff_file)
    if(toFile){
      wrassp::write.AsspDataObj(dobj=outDataObj,file=ssff_file)
      #Here we can be sure that the list is a valid SSFF object, so the
      # so we add TRUE to the out vector
      outListOfFiles <- c(outListOfFiles,TRUE)
    }
    
  }
  
  if(toFile){
    return(length(outListOfFiles))
  }else{
    return(outDataObj)
  }
  
}
attr(praat_sauce,"ext") <-  c("psa") 
attr(praat_sauce,"tracks") <-  c("f0","fm", "bw","H","Hc","A","Ac","H25K","HH","HHc","HA","HAc","cpp","hnr")



# FOR INTERACTIVE TESTING
#Use this to mark that a Praat script is being developed (and the version in the installed copy of the library can not be used)
# PRAAT_DEVEL = TRUE
# 
# path2demoData = file.path(tempdir(),"emuR_demoData")
# unlink(path2demoData, recursive = TRUE)
# 
# emuR::create_emuRdemoData()
# 
# ae <- emuR::load_emuDB(file.path(path2demoData,"ae_emuDB"))
# 
# 
# add_trackDefinition(ae,
#                      name="cpp",
#                      fileExtension ="cpp",columnName = "cpp",
#                      onTheFlyFunctionName = "praat_sauce",onTheFlyParams = list(pitchTracking=TRUE,formantMeasures=TRUE,spectralMeasures=TRUE))
# add_perspective(ae,"Praat")
# 
# set_levelCanvasesOrder(ae,"Praat",c("Phonetic","Tone"))
# set_signalCanvasesOrder(ae,"Praat",c("OSCI","SPEC","fm","cpp"))

#rm(dbConfig)
#navigateToFile(file.path(ae$basePath,"ae_DBconfig.json"),line=395)
#httpuv::stopAllServers()
#serve(ae,autoOpenURL = NULL)


# library('testthat')
# test_file('tests/testthat/test_aaa_initDemoDatabase.R')
# test_file('tests/testthat/test_praat.R')




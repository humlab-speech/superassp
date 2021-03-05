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

# wrassp::wrasspOutputInfos -> wrasspOutputInfos
# 
# wrasspOutputInfos[["praat_formant_burg"]] <- wrasspOutputInfos[["forest"]]

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

  #Use this to mark that the Praat script is being developed
  PRAAT_DEVEL = TRUE
  
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
        columnName <- ifelse(is.null(columnName),columnName,attr(fun,"tracks")[[1]])
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
      add_ssffTrackDefinition(emuDBhandle,name=name,columnName = columnName,fileExtension = ext)
    }else{
      
      stop("Could not find a definition of the function ",onTheFlyFunctionName,"." )
    }
  }
}


enable_formantOverlay <- function(emuDBhandle,perspective){
  perspectiveNames <- list_perspectives(emuDBhandle)$name
  trackNames <- list_ssffTrackDefinitions(emuDBhandle)$name
  
  #Stop processing if the perspective is not defined in the database
  if(! perspective %in% perspectiveNames) {stop("The perspective  ",perspective," is not defined in the database ", emuDBhandle$dbName,"!")}
  
  #Stop processing if a track FORMANTS is not defined in the database
  if(! "FORMANTS" %in% trackNames) {stop("In order to enable formant overlays, a track named 'FORMANTS' must be defined in the database !")}
  
  which(grepl(perspective,perspectiveNames)) -> perspid
  dbConfig = emuR:::load_DBconfig(ae)
  
  dbConfig$EMUwebAppConfig$perspectives[[perspid]]$signalCanvases$assign[[1]] <- list("signalCanvasName"="SPEC","ssffTrackName"="FORMANTS")
  res <- emuR:::store_DBconfig(emuDBhandle,dbConfig = dbConfig)
  return(res)
}

set_overlayTrack <- function(emuDBhandle,perspective,trackname, overlay.on="SPEC",overwrite=FALSE){
  perspectiveNames <- list_perspectives(emuDBhandle)$name
  trackNames <- list_ssffTrackDefinitions(emuDBhandle)$name
  
  #Stop processing if the perspective is not defined in the database
  if(! perspective %in% perspectiveNames) {stop("The perspective  ",perspective," is not defined in the database ", emuDBhandle$dbName,"!")}
  
  #Stop processing if the track is not defined in the database
  if(! trackname %in% trackNames) {stop("The track  ",trackname," is not defined in the database ", emuDBhandle$dbName,"!")}
  
  which(grepl(perspective,perspectiveNames)) -> perspid
  dbConfig = emuR:::load_DBconfig(ae)
  
  #Check and stop processing if an overlay is alrady set 
  overlay <- dbConfig$EMUwebAppConfig$perspectives[[perspid]]$signalCanvases$assign
  
  if(length(overlay) > 0 ){
     for(ov in 1:length(overlay)){
       if(overlay[[ov]]$signalCanvasName == overlay.on ){
         if(! overwrite) {stop("Cannot set an overlay on ", overlay.on, " as one is already defined in the database.\nPlease set overwrite=TRUE if you wish to overwrite the previous setting.")}
         
         dbConfig$EMUwebAppConfig$perspectives[[perspid]]$signalCanvases$assign[[ov]] <- list("signalCanvasName"=overlay.on,"ssffTrackName"= trackname)
       }else{
         #In this case, existing overlay settings do not exist
         dbConfig$EMUwebAppConfig$perspectives[[perspid]]$signalCanvases$assign[[ov+1]] <- list("signalCanvasName"=overlay.on,"ssffTrackName"= trackname)
       }
     }
  }

  res <- emuR:::store_DBconfig(emuDBhandle,dbConfig = dbConfig)
  return(res)
}


# FOR INTERACTIVE TESTING

path2demoData = file.path(tempdir(),"emuR_demoData")
unlink(path2demoData, recursive = TRUE)

emuR::create_emuRdemoData()

ae <- emuR::load_emuDB(file.path(path2demoData,"ae_emuDB"))

add_trackDefinition(ae,
                    name="FORMANTS",
                    fileExtension ="pfm",columnName = "fm",
                    onTheFlyFunctionName = "praat_formant_burg",onTheFlyParams = list(maxhzformant=5000.0))
add_perspective(ae,"Praat")
add_ssffTrackDefinition(ae,"pitch",fileExtension = "pitch",columnName = "pitch",onTheFlyFunctionName = "mhsF0")
set_levelCanvasesOrder(ae,"Praat",c("Phonetic","Tone"))
set_signalCanvasesOrder(ae,"Praat",c("OSCI","SPEC","fm","FORMANTS"))

enable_formantOverlay(ae,"Praat")
set_overlayTrack(ae,"Praat","pitch","OSCI")

#rm(dbConfig)
#navigateToFile(file.path(ae$basePath,"ae_DBconfig.json"),line=395)
httpuv::stopAllServers()
#serve(ae,autoOpenURL = NULL)

# library('testthat')
# test_file('tests/testthat/test_aaa_initDemoDatabase.R')
# test_file('tests/testthat/test_praat.R')




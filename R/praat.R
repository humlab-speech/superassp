

#' A simple check of a presence of a Praat executable
#' 
#' 
#'
#' @param praat_path A character string containing the path to the executable that the function was able to find (or the executable that the function was able to verify the existance of), or NULL if no Praat executable was found and verified.
#'
#' @return A boolean indicating whether the Praat executable could be found or not.
#' @export
#'
has_praat <- function(praat_path=NULL){
  return(ifelse(is.null(get_praat(praat_path)) ,FALSE,TRUE))
}



#' Utility function for getting the full path of the Praat executable.
#' 
#' This function checks the system variables and deduces where Praat is installed. On the OSX platform (Darwin) the Praat app is assumed to exist in the Applications folder, 
#' and the actual binary inside of the application package is then used. If not OSX, then the function will search the default search paths for executables set up in the OS. 
#' If an explicit path is given, then function will just check whether the executable is actualy present there.
#'
#' @param praat_path 
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
#' \code{\link{has_praat}} functions.
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
#' @param praat_path give an explicit path for Praat. If the praat 
#'
#' @return a list of 
#' @export
#'
#' @examples
praat_formant_burg <- function(listOfFiles,beginTime=0,endTime=0,windowShift=0.0,numFormants=5.0,maxhzformant=5500.0,windowSize=0.025,preemphasis=50.0,window="hanning",relativeWidth=1.0,toFile=TRUE,explicitExt="fms",outputDirectory=NULL,praat_path=NULL){
  
  if(! has_praat(praat_path)){
    stop("Could not find praat. Please specify a full path.")
  }
  
  formant_burg <- tjm.praat::wrap_praat_script(praat_location = get_praat(),
                                    script_code_to_run = readLines(file.path(
                                      system.file(package = "superassp",mustWork = TRUE),"praat","formant_burg.praat"))
                                    ,return="last-argument")
  
  #Check that all files exists before we begin
  filesEx <- file.exists(listOfFiles)
  if(!all(filesEx)){
    filedNotExists <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ",paste(filedNotExists, collapse = ", "))
  }
  #The empty vector of file names that should be returned
  outListOfFiles <- c()
  
  for(soundFile in listOfFiles){ 
    
    outTabFile <- tempfile(fileext = ".Table")
    
    tabfile <- formant_burg(soundFile,
                            beginTime,
                            endTime,
                            windowShift,
                            numFormants,
                            maxhzformant,
                            windowSize,
                            preemphasis,
                            window,
                            relativeWidth,
                            outTabFile)
    
    # We need the sound file to extract some information
    
    origSound <- wrassp::read.AsspDataObj(soundFile)
    
    inTable  <- readr::read_csv(tabfile,na=c("--undefined--",""),col_types=list(col_double(),
                                                                           readr::col_double(),
                                                                           readr::col_double(),
                                                                           readr::col_double(),
                                                                           readr::col_double(),
                                                                           readr::col_double(),
                                                                           readr::col_double(),
                                                                           readr::col_double(),
                                                                           readr::col_double(),
                                                                           readr::col_double(),
                                                                           readr::col_double(),
                                                                           readr::col_double(),
                                                                           readr::col_double(),
                                                                           readr::col_double()))
  #  if(is.null(columns)){columns <- seq(5,ncol(inTable),2)}
    
    starTime = inTable[1,"time(s)"]
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16", "INT16")
    #Use the time separation between second and first formant measurement time stamps to compute a sample frequency.
    sampleRate <-  as.numeric(1 / (inTable[2,"time(s)"] - inTable[1,"time(s)"]))
    attr(outDataObj, "sampleRate") <- sampleRate
      
    attr(outDataObj, "origFreq") <-  as.numeric(attr(origSound, "sampleRate"))
    startTime <- as.numeric(inTable[1,"time(s)"])
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
    class(outDataObj) = "AsspDataObj"
    
    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    fmTable <- inTable %>%
      dplyr::select(starts_with("F",ignore.case = FALSE)) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer)) 
    
    noFormantsValues <- nrow(fmTable)
    
    outDataObj = wrassp::addTrack(outDataObj, "fm", as.matrix(fmTable), "INT16")
    
    bwTable <- inTable %>%
      dplyr::select(starts_with("F",ignore.case = FALSE)) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer)) 
    
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
                               ncol = 5) #ncol(outDataObj$fm)
  
  
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
    #Here we can be sure that the list is a valid SSFF object, so the 
    # so we add TRUE to the out vector
    outListOfFiles <- c(outListOfFiles,TRUE)
    if(toFile){
      ssff_file <- gsub("wav$",explicitExt,soundFile)
      if(!is.null(outputDirectory)){
        ssff_file <- file.path(outputDirectory,basename(ssff_file))
      }
      
      attr(outDataObj,"filePath") <- as.character(ssff_file)
      wrassp::write.AsspDataObj(dobj=outDataObj,file=ssff_file)
    }
    
  }
  
  #Return a summary indicating whether all files were successfully created
  return(all(outListOfFiles))
}


# FOR INTERACTIVE TESTING
library('testthat')
test_file('tests/testthat/test_aaa_initDemoDatabase.R')
test_file('tests/testthat/test_praat.R')




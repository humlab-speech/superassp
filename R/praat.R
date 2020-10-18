

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
#' 
#' 
praat_formant_burg <- function(listOfFiles,beginTime=0,endTime=0,windowShift=0.0,numFormants=5.0,maxhzformant=5500.0,windowSize=0.025,preemphasis=50.0,window="hanning",relativeWidth=1.0,toFile=TRUE,explicitExt="fms",outputDirectory=NULL,praat_path=NULL){
  
  if(! has_praat(praat_path)){
    stop("Could not find praat. Please specify a full path.")
  }
  
  tryCatch({
    fileBeginEnd <- data.frame(
      listOfFiles = listOfFiles, 
      beginTime = beginTime,
      endTime=endTime
      )
 },error=function(e){stop("The beginTime and endTime must either be a single value or the same length as listOfFiles")})
  
  
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
  
  for(i in 1:nrow(fileBeginEnd)){ 
    soundFile <- fileBeginEnd[i, "listOfFiles"]
    beginTime <- fileBeginEnd[i, "beginTime"]
    endTime <- fileBeginEnd[i, "endTime"]
    
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
    
    inTable  <- readr::read_csv(tabfile,na=c("--undefined--",""),col_types=list(readr::col_double(),
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
      dplyr::select(tidyselect::starts_with("F",ignore.case = FALSE)) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer)) 
    
    noFormantsValues <- nrow(fmTable)
    
    outDataObj = wrassp::addTrack(outDataObj, "fm", as.matrix(fmTable), "INT16")
    
    bwTable <- inTable %>%
      dplyr::select(tidyselect::starts_with("F",ignore.case = FALSE)) %>%
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


#' Preform a Praat voice report on a sound sample
#' 
#' This function sends a sounds file to Praat to perform a voice analys (Voice Report) and return the results.
#' It is possible to specify the part of sound signal that should be sent to 
#' the underlying Praat function in a manner that that is congruent with 
#' \code{wrassp} functions (e.g. \code{\link[wrassp]{forest}}, 
#' \code{\link[wrassp]{acfana}} and so on) so that it could work well wihtin for instance the \code{emuR} package. 
#' 
#' In addition, it is possible to specify where for instance a vowel measures should be made. The boundaries of the vowel are in this case given using the \code{beginTime} and \code{endTime} parameters as is usual for  \code{wrassp} functions. 
#' Then, the user may specify an offset (\code{selectionOffset}) from the \code{beginTime} time where measures should begin, and the length (\code{selectionDuration}) of the part of vowel that the user wants to analyse (2s or 3s are usual values). 
#' 
#' While it would be the equivalent to specify the start and end points of the
#' part of the signal that should be analysed directly via the \code{beginTime} 
#' and \code{endTime} parameters, this possibility of separating the "where the vowel is located" and "what part of the vowel should be extracted for analysis" makes using this function easier than calculating the start and end times of analysis directly for each vowel sample, and should work better with \code{emuR}.
#' 
#' The function will return all parameters in wide format 
#' (one column per returned voice parameter) per default, but may return wide format also if \code{returnWide=TRUE} is specified
#'
#' @param listOfFiles The sound files that should be processed.
#' @param beginTime the start time (in s) of the unit to be extracted
#' @param endTime  end end (in s) of the unit to be extracted
#' @param selectionOffset an offset that allow for starting the analysis for instance 1s into the vowel.
#' @param selectionDuration the duration (counted from the offset, in s) of the part of signal that should be analysed.
#' @param window the type of window function that should be used when extracting the signal for analysis
#' @param relativeWidth the width of the window.
#' @param returnWide whether to return a wide format \code{tibble} or a long format (see the return value below)
#' @param praat_path the full path of the Praat binary.
#'
#' @return Computed voice properties (either in wide or long format) with the 
#'   \itemize{              
#'   \item Selection start
#'   \item Selection end
#'   \item Vowel start
#'   \item Vowel end
#'   \item Median Pitch
#'   \item Mean Pitch
#'   \item Pitch SD
#'   \item Min Pitch
#'   \item Max Pitch
#'   \item Number Of Pulses
#'   \item Number Of Periods
#'   \item Mean period
#'   \item Period SD
#'   \item Frac local unvoiced frames
#'   \item Voice breaks
#'   \item Degree voice breaks
#'   \item Jitter (local)
#'   \item Jitter (local, absolute)
#'   \item Jitter (rap)
#'   \item Jitter (ppq5)
#'   \item Jitter (ddp)
#'   \item Shimmer (local)
#'   \item Shimmer (local, absolute)
#'   \item Shimmer (apq3)
#'   \item Shimmer (apq5)
#'   \item Shimmer (apq11)
#'   \item Shimmer (dda)
#'   \item Mean Autocorrelation
#'   \item Mean noise-to-harmonics ratio
#'   \item Mean harmonics-to-noise ratio
#'   \item Mean intensity
#'   \item Median intensity
#'   \item Intensity standard deviation
#' }
#' @export
#'

praat_voice_report <- function(listOfFiles,beginTime=0,endTime=0,selectionOffset=0.0,selectionDuration=2.0,window="hanning",relativeWidth=1.0,returnWide=TRUE,praat_path=NULL){
  
  voice_report <- tjm.praat::wrap_praat_script(
    praat_location = get_praat(),
    script_code_to_run = readLines(
      file.path("inst","praat","praat_voice_report.praat")),
    return="info-window")
  
  for(currFile in listOfFiles){
    voice_report(currFile,
                 beginTime,
                 endTime,
                 selectionOffset,
                 selectionDuration,window,relativeWidth) -> info
  }
  values <- suppressWarnings(
    as.numeric(
      unlist(
        str_split(gsub("\\\\n","",
                       gsub("--undefined--","NA",info[2])),";"))))
  measures <- str_trim(unlist(str_split(gsub("\\\\n","",info[1]),";")),"both")
  
  data.frame(Measures=measures,Values=values) -> out
  if(returnWide){
    suppressWarnings({
      out <- out %>%
      tidyr::pivot_wider(names_from = "Measures",
                         values_from="Values",
                         values_fill=NA)
    }
    )
  }
  
     
  
  
  return(out)
}



# FOR INTERACTIVE TESTING
# library('testthat')
# test_file('tests/testthat/test_aaa_initDemoDatabase.R')
# test_file('tests/testthat/test_praat.R')




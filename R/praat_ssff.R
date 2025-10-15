
#' Make a Pitch object from an SSFF f0 track
#'
#' This function takes an SSFF object or a path to a file containing one, reads
#' the content of a specific `field` (and `channel`/column if required) and uses
#' Praat to construct a Pitch object and store it in a file. The function's primary purpose is to allow the user to re-use already computed, and possibly adjusted, pitch tracks in subsequent calculations in for instance [superassp::praat_sauce]. 
#' 
#' The user may optionally compute a Pitch object for only parts of the SSFF pitch track if the track is stored in a file. _An ability to subset an SSFF object directly is currently not implemented_.
#'
#'
#' @param inData An SSFF f0 track object, or the full path to one.
#' @param outputPath The directory where the Pitch file should be stored.
#' @param field The field / column in the SSFF object where the f0 track values
#'   are stored. The field may be indicated by the name of it, or its number.
#'   Often, the SSFF object will contain just one field, the one containing f0
#'   values, and `field=1` is therefore a good default.
#' @param channel The channel / column of the field which should be used to
#'   construct the pitch track. It is likely very uncommon that `channel` values
#'   higher than one will be used in real world applications.
#' @param start An optional start time (in s) for a part of the SSFF f0 track
#'   which should be converted to a Pitch object.
#' @param end An optional end time (in s) for a part of the SSFF f0 track which
#'   should be converted to a Pitch object.
#' @param zero.threshold The threshold below which f~0~ values will be considered to indicate an unvoiced frame. The default is zero, which is how the SSFF format and Praat pitch tracks encode lack of voicing. It is however also possible to trim away really low f$_0$ values before making the Pitch object by choosing a higher threshold.
#' @param dump.script If `TRUE`, the Praat script that is used to create a Pitch object will be dumped to a file for debuging or inspection. The file will be placed in the same directory as the output Pitch file, but with a '.praat'. extension.
#' @param soundFileDuration An explicit duration that will be set for the created
#'   Pitch object. If not explicitly supplied, the duration of the SSFF object track will be used instead.
#'
#' @return The path to the created Pitch file
#' @export
#' 
#' @seealso [ksvF0]
#' @seealso [superassp::praat_sauce]
#' @seealso [mhsF0]
#'

ssffToPitch <- function(inData,outputPath=NULL,field=1,channel=1,start=0.0,end=0.0,zero.threshold=0.0,dump.script=FALSE,soundFileDuration=NULL){

  if(is.character(inData) && file.exists(normalizePath(inData))){
    inData <- read.AsspDataObj(inData,begin=start,end=end)
  }
  #The only other option is an SSFF object
  if(! is.AsspDataObj(inData)){
    stop("The function toPitch requires either an SSFF or a path to an SSFF file")
  }
  
  if(is.character(field)){
    if(! field %in% names(inData) ){
      stop("Unable find the '",field,"' track in the SSFF file.")
    }  
  }else{
    if(is.integer(field) && field > length(inData)){
      stop("The SSFF file contain fewer than ",field," tracks.")
    }
  }
  if(channel > ncol(inData[[field]])){
    stop("There are no channel '",channel,"' in the SSFF track file.  Please give a valid value.")
  }
  pitchValues <- inData[[field]]
  sr <- attr(inData,"sampleRate")
  er <- attr(inData,"endRecord")
  
  
  pitchValues[pitchValues <= zero.threshold] <- NA
  maxF <- max(pitchValues,na.rm=TRUE)
  minF <- min(pitchValues,na.rm=TRUE)
  stepSize = 1 / sr
  times <- seq(attr(inData,"startRecord")-1,attr(inData,"endRecord")-1,1) * 1/ attr(inData,"sampleRate") + attr(inData,"startTime")

  soundFileDuration <- ifelse(is.null(soundFileDuration),
                              dur.AsspDataObj(inData), 
                              soundFileDuration)
  df <- data.frame(time=times,pitch=pitchValues)
  #Make a set of pitch values that should be placed on the pitch tier
  dfCompl <- na.omit(df)
  # These time frames will be unvoiced later once the Pitch object has been created
  dfNA <- df[is.na(df$pitch), ]

  #The devoicing formula uses a time window 
  devoiceStart <- dfNA$time - attr(inData,"startTime")
  devoiceEnd <- dfNA$time + attr(inData,"startTime")
  fp <- attr(s0,"filePath")
  outPath <- paste(tools::file_path_sans_ext(fp),"Pitch",sep=".")
  
  scriptPath <- paste(tools::file_path_sans_ext(fp),"praat",sep=".")
  
  #BEGIN Boilerplate
  formSetup <- paste("form Pitch creation procedure from a vector of times and values", " real SoundFileDuration 1.0", " real TimeStep 0.005"," real PitchFloor 60.0"," real PitchCeiling 400.0", "endform",sep="\n")
  ptSetup <- 'Create PitchTier: "newPitch", 0.0 , soundFileDuration'
  ptConvert <- 'To Pitch: timeStep, pitchFloor, pitchCeiling'
  writeFile <- paste0('Save as text file: "',outPath,'"',"\n")
  #END Boilerplate

  pitchPoints<- paste("Add point:",dfCompl[,"time"],",",dfCompl[,"pitch"])
  
  ptDevoice <- paste0("Formula: ~ if x>",devoiceStart," and x<",devoiceEnd," then 0 else self fi")

  script <- c(formSetup,ptSetup,pitchPoints,ptConvert,ptDevoice, writeFile)
  if(dump.script){
    readr::write_lines(script,scriptPath)
    script <- readr::read_lines(scriptPath)
  }
  makePitch <- tjm.praat::wrap_praat_script(praat_location = get_praat(),
                                               script_code_to_run = script
                                               ,return="last-argument")
  makePitch(soundFileDuration,stepSize,minF,maxF)
  
  return(outPath)
}


#' Make a Praat Formant object from an SSFF track track
#' 
#' This function takes an SSFF object or a path to a file containing one, reads
#' the content of a `fm.field` and `bw.field` and uses
#' Praat to construct a Formant object and store it in a file. The function's primary purpose is to allow the user to re-use already computed, and possibly adjusted, formant frequency and bandwidth tracks in subsequent calculations in for instance [superassp::praat_sauce]. 
#' 
#' The user may optionally compute a Formant object for only parts of the SSFF formant track if the track is stored in a file. _An ability to subset an SSFF object directly is currently not implemented_.
#' 
#' Please note that the process of converting formant tracks from SSFF is _very_ slow, as shown by the 
#' [microbenchmark] output below which indicates that constructing an Formant object using the current strategy takes about 
#' 353 times the time it takes to simply compute the formant tracks. 
#' 
#' ```
#' Unit: relative
#' expr         min          lq        mean      median          uq         max neval
#' read      1.0000      1.0000      1.0000      1.0000      1.0000      1.0000     1
#' forest    705.8529    705.8529    705.8529    705.8529    705.8529    705.8529     1
#' Formant 249267.5829 249267.5829 249267.5829 249267.5829 249267.5829 249267.5829     1
#' ```
#' 
#' @param inData An SSFF formant track object, or the full path to one.
#' @param outputPath The directory where the Formant file should be stored.
#' @param fm.field The field / column in the SSFF object where the formant frequency tracks
#'   are stored. The field may be indicated by the name of it, or its number.
#'   Often, the SSFF object will contain formant frequency values in the first field and bandwidths in the second,
#'    so the user could also give `field=1` as an argument.
#' @param bw.field The field / column in the SSFF object where the formant banwidth tracks
#'   are stored. The field may be indicated by the name of it, or its number.
#'   Often, the SSFF object will contain formant frequency values in the first field and bandwidths in the second,
#'    so the user could also give `field=2` as an argument.
#' @param start An optional start time (in s) for a part of the SSFF formant track
#'   which should be converted to a Formant object.
#' @param end An optional end time (in s) for a part of the SSFF formant track which
#'   should be converted to a Formant object.
#' @param windowShift The analysis window of the formant track (in ms).
#' @param nominalF1 An assumed F$_1$ value for the created Formant object.
#' @param zero.threshold The threshold below which formant frequencies or bandwidths should be considered as zero. 
#' @param dump.script If `TRUE`, the Praat script that is used to create a Pitch object will be dumped to a file for debuging or inspection. The file will be placed in the same directory as the output Formant file, but with a '.praat'. extension.
#' @param soundFileDuration An explicit duration that will be set for the created
#'   Formant object. If not explicitly supplied, the duration of the SSFF object track will be used instead.
#'   
#' @export
#' @seealso [superassp::praat_sauce]
#' @seealso [forest]
#' @seealso [superassp::praat_formant_burg]

ssffToFormant <- function(inData,outputPath=NULL,fm.field="fm",bw.field="bw",start=0.0,end=0.0,windowShift=5,nominalF1=500,zero.threshold=0.0,dump.script=FALSE,soundFileDuration=NULL){
  
  if(is.character(inData) && file.exists(normalizePath(inData))){
    inData <- read.AsspDataObj(inData,begin=start,end=end)
  }
  #The only other option is an SSFF object
  if(! is.AsspDataObj(inData)){
    stop("The function toPitch requires either an SSFF or a path to an SSFF file")
  }
  #Frequencies
  if(is.character(fm.field)){
    if(! fm.field %in% names(inData) ){
      stop("Unable find the '",fm.field,"' track in the SSFF file.")
    }  
  }else{
    if(is.integer(fm.field) && fm.field > length(inData)){
      stop("The SSFF file contain fewer than ",fm.field," tracks.")
    }
  }
  #Bandwidths
  if(is.character(bw.field)){
    if(! bw.field %in% names(inData) ){
      stop("Unable find the '",bw.field,"' track in the SSFF file.")
    }  
  }else{
    if(is.integer(bw.field) && bw.field > length(inData)){
      stop("The SSFF file contain fewer than ",bw.field," tracks.")
    }
  }
  
  freqValues <- inData[[fm.field]]
  bwValues <- inData[[bw.field]]
  nFormants <- ncol(freqValues)
  sr <- attr(inData,"sampleRate")
  er <- attr(inData,"endRecord")
  stepSize = 1 / sr
  times <- seq(attr(inData,"startRecord")-1,attr(inData,"endRecord")-1,1) * 1/ attr(inData,"sampleRate") + attr(inData,"startTime")
  
  soundFileDuration <- ifelse(is.null(soundFileDuration),
                              dur.AsspDataObj(inData), 
                              soundFileDuration)


  fp <- attr(s0,"filePath")
  outPath <- paste(tools::file_path_sans_ext(fp),"Formant",sep=".")
  scriptPath <- paste(tools::file_path_sans_ext(fp),"praat",sep=".")
  base <- basename(fp)
  
  
  #BEGIN Boilerplate
  formSetup <- paste("form Formant creation procedure from vectors of times and formant frequency and bandwidth values",  "sentence Name Hejsan","real Start_time_(s) 0.0","real End_time_(s) 1.0","real Time_step_(s) 0.01","integer Number_of_formants: 5","real Initial_first_formant_(Hz) 550.0","real Initial_formant_spacing_(Hz) 1100.0","real Initial_first_bandwidth_(Hz) 60.0","real Initial_badwidth_spacing_(Hz) 50.0","real Intensity_(Pa²) 0.1", "endform",sep="\n")
  ptSetup <- 'Create FormantGrid: name$, start_time, end_time, number_of_formants, initial_first_formant, initial_formant_spacing, initial_first_bandwidth, initial_badwidth_spacing'
  ptConvert <- 'To Formant: time_step, intensity'
  writeFile <- paste0('Save as text file: "',outPath,'"',"\n")
  #END Boilerplate
  
  formantPoints<- c()
  # This makes sense since contrary to a pitch track, for instance, you can have
  # an F1 cutback and then all other formants but F1 is defined in a frame. It makes no sense then to 
  # remove the entier frame (like what we do when creating a Pitch object).
  for(r in 1:nrow(freqValues)){
    for(c in 1:nFormants){

      if(freqValues[r,c] > zero.threshold){
        formantPoints <- c(formantPoints,
                           paste0("Add formant point: ",c,", ",times[r],", ", freqValues[r,c])
        )
      }
      if(! is.na(bwValues[r,c] ) && bwValues[r,c] > zero.threshold ){
        #The check above is likelly unecessary. It should not happen that a formant is defined as 
        # a frequency but not in terms of bandwidth.
        formantPoints <- c(formantPoints,
                           paste0("Add bandwidth point: ",c,", ",times[r],", ", bwValues[r,c])
                           )
      }
    
    }
  }

  removeFormantPoints<- c()
  # We need to remove formants explicitly once synthesized as a formant track
  for(r in 1:nrow(freqValues)){
    for(c in 1:nFormants){
      
      if(freqValues[r,c] <= zero.threshold){
        removeFormantPoints <- c(removeFormantPoints,
                           paste0("n = Get frame number from time: ",times[r],"\n",'Formula (frequencies): "if row = ',c,' and col = \'n\' then 0 else self fi"'),
                           paste0("\n",'Formula (bandwidths): "if row = ',c,' and col = \'n\' then 0 else self fi"')
                                  
        )
      }

      
    }
  }
  
  
  script <- c(formSetup,ptSetup,formantPoints,ptConvert,removeFormantPoints,writeFile)
  if(dump.script){
    readr::write_lines(script,scriptPath)
    script <- readr::read_lines(scriptPath)
  }
  makeFormants <- tjm.praat::wrap_praat_script(praat_location = get_praat(),
                                            script_code_to_run = script
                                            ,return="last-argument")
  #  formSetup <- paste("form Formant creation procedure from vectors of times and formant frequency and bandwidth values",  "sentence Name Hejsan","real Start_time_(s) 0.0","real End_time_(s) 1.0","real Time_step_(s) 0.01","integer Number_of_formants: 5","real Initial_first_formant_(Hz) 550.0","real Initial_formant_spacing_(Hz) 1100.0","real Initial_first_bandwidth_(Hz) 60.0","real Initial_badwidth_spacing_(Hz) 50.0","real Intensity_(Pa²) 0.1", "endform",sep="\n")
  makeFormants(base, start, soundFileDuration, windowShift/1000, nFormants, nominalF1,1100,60,50, 0.1 )
  
  return(outPath)
}

# forest("~/Desktop/kaa_yw_pb.wav")
# f <- read.AsspDataObj("~/Desktop/kaa_yw_pb.fms")
# ssffToFormant("~/Desktop/kaa_yw_pb.fms", dump.script = FALSE) -> out

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
















# FOR INTERACTIVE TESTING
#library('testthat')
# test_file('tests/testthat/test_aaa_initDemoDatabase.R')
# test_file('tests/testthat/test_praat.R')

#str(praat_pitch("~/Desktop/kaa_yw_pb.wav",toFile=FALSE))
#praat_pitch("~/Desktop/short_aaa.wav",toFile=FALSE,corr.only = FALSE) -> f0

# library(dplyr)
# library(logger)
# logger::log_threshold(logger::INFO)
# PRAAT_DEVEL <- TRUE
# 
# praat_formant_burg("~/Desktop/aaa.wav",toFile=FALSE) -> out
# praat_formantpath_burg("~/Desktop/aaa.wav",toFile=FALSE) -> out
# praat_sauce("~/Desktop/aaa.wav",toFile=FALSE) -> out
# praat_intensity("~/Desktop/aaa.wav",toFile=FALSE) -> out
# praat_pitch("~/Desktop/aaa.wav",toFile=FALSE) -> out
# praat_moments("~/Desktop/aaa.wav",toFile=FALSE) -> out




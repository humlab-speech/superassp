
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
#' @param zero.threshold The threshold below which f$_0$ values will be considered to indicate an unvoiced frame. The default is zero, which is how the SSFF format and Praat pitch tracks encode lack of voicing. It is however also possible to trim away really low f$_0$ values before making the Pitch object by choosing a higher threshold.
#' @param dump.script If `TRUE`, the Praat script that is used to create a Pitch object will be dumped to a file for debuging or inspection. The file will be placed in the same directory as the output Pitch file, but with a '.praat'. extension.
#' @param soundFileDuration An explicit duration that will be set for the created
#'   Pitch object. If not explicitly supplied, the duration of the SSFF object track will be used instead.
#'
#' @return The path to the created Pitch file
#' @export
#' 
#' @seealso [wrassp::ksvF0]
#' @seealso [superassp::praat_sauce]
#' @seealso [wrassp::mhsF0]
#'

ssffToPitch <- function(inData,outputPath=NULL,field=1,channel=1,start=0.0,end=0.0,zero.threshold=0.0,dump.script=FALSE,soundFileDuration=NULL){

  if(is.character(inData) && file.exists(normalizePath(inData))){
    inData <- wrassp::read.AsspDataObj(inData,begin=start,end=end)
  }
  #The only other option is an SSFF object
  if(! wrassp::is.AsspDataObj(inData)){
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
                              wrassp::dur.AsspDataObj(inData), 
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
#' @seealso [wrassp::forest]
#' @seealso [superassp::praat_formant_burg]

ssffToFormant <- function(inData,outputPath=NULL,fm.field="fm",bw.field="bw",start=0.0,end=0.0,windowShift=5,nominalF1=500,zero.threshold=0.0,dump.script=FALSE,soundFileDuration=NULL){
  
  if(is.character(inData) && file.exists(normalizePath(inData))){
    inData <- wrassp::read.AsspDataObj(inData,begin=start,end=end)
  }
  #The only other option is an SSFF object
  if(! wrassp::is.AsspDataObj(inData)){
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
                              wrassp::dur.AsspDataObj(inData), 
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

# wrassp::forest("~/Desktop/kaa_yw_pb.wav")
# f <- wrassp::read.AsspDataObj("~/Desktop/kaa_yw_pb.fms")
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



#' Formant estimation using the Praat burg algorithm implementation
#' 
#' Formants are estimated using Praat's built in function (burg algorithm). The function also computes the intensity (L) of the formant based on the power of the spectrum at the frequency of the formant. Naturally, if the algorithm failed to find a formant in a specified time frame, then the function will not return a formant frequency, bandwidth and intensity estimation.
#' 
#' If the user only want to estimate formant frequencies, computing them using the function [wrassp::forest] is _much_ quicker, and the user should therefore mainly consider using this function `praat_formant_burg` only if the use case specifically demands the use of the burg algorithm for computing formants, or if the user wants to also study the formant intensity levels (L_n) which  [wrassp::forest] does not do.
#'  
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
#' @param preemphasis the frequency from which a preemphasis will be applied.
#' @param trackFormants boolean; Should Praat attempt to gather short time formant frequency estimates into tracks?
#' @param numberOfTracks The number of tracks to follow (if `trackFormants` is `TRUE`), and the number of tracks in the output. Information on frequencies bandwidths of formants with numbers above  `numberOfTracks` will be discarded.
#' @param nominalF1 Described by the Praat manual as the preferred value near which the first track wants to be. For average (i.e. adult female) speakers, this value will be around the average F1 for vowels of female speakers, i.e. 550 Hz.
#' @param nominalF2 Described by the Praat manual as the preferred value near which the second track wants to be. A good value will be around the average F2 for vowels of female speakers, i.e. 1650 Hz.
#' @param nominalF3 Described by the Praat manual as the preferred value near which the third track wants to be. A good value will be around the average F3 for vowels of female speakers, i.e. 2750 Hz. This argument will be ignored if you choose to have fewer than three tracks, i.e., if you are only interested in F1 and F2.
#' @param frequencyCost Described by the Praat manual as the preferred value near which the five track wants to be. In the unlikely case that you want five tracks, a good value may be around 4950 Hz.Frequency cost (per kiloHertz)
#' @param bandwidthCost Described by the Praat manual as the local cost of having a bandwidth, relative to the formant frequency. For instance, if a candidate has a formant frequency of 400 Hz and a bandwidth of 80 Hz, and Bandwidth cost is 1.0, the cost of having this formant in any track is (80/400) · 1.0 = 0.200. So we see that the procedure locally favours the inclusion of candidates with low relative bandwidths.
#' @param transitionCost Described by the Praat manual as the cost of having two different consecutive formant values in a track. For instance, if a proposed track through the candidates has two consecutive formant values of 300 Hz and 424 Hz, and Transition cost is 1.0/octave, the cost of having this large frequency jump is (0.5 octave) · (1.0/octave) = 0.500.
#' @param windowShape the analysis window function used when extracting part of a sound file for analysis. De faults to "Hanning".
#' @param relativeWidth the relative width of the windowing function used.
#' @param spectWindowShape The shape of the windowing function used for constructing the spectrogram. 
#' @param spectResolution The frequency resolution of the spectrogram from which formant intensities will be collected.
#' @param toFile write the output to a file? The file will be written in  `outputDirectory`, if defined, or in the same directory as the soundfile. 
#' @param explicitExt the file extension that should be used.
#' @param outputDirectory set an explicit directory for where the signal file will be written. If not defined, the file will be written to the same directory as the sound file.
#' @param verbose Not implemented. Only included here for compatibility.  
#' @param praat_path give an explicit path for Praat. 
#'
#' @return Ar an SSFF track data object (if `toFile=FALSE`) containing three fields ("F", "B" and "L") containing formant frequencies, bandwidths and intensities.
#' 
#' @export
#' @seealso [wrassp::forest]
#' @seealso [superassp::praat_formantpath_burg]


praat_formant_burg <- function(listOfFiles,
                               beginTime=0,
                               endTime=0,
                               windowShift=5.0,
                               numFormants=5.0,
                               maxFormantHz=5500.0,
                               windowSize=30,
                               preemphasis=50.0,
                               trackFormants=TRUE,
                               numberOfTracks=3,
                               nominalF1=550,
                               nominalF2=1650,
                               nominalF3=2750,
                               frequencyCost=1.0,
                               bandwidthCost=1.0,
                               transitionCost=1.0,
                               windowShape="Gaussian1",
                               relativeWidth=1.0,
                               spectWindowShape="Gaussian",
                               spectResolution=40.0,
                               toFile=TRUE,
                               explicitExt="pfm",
                               outputDirectory=NULL,
                               verbose=FALSE,
                               praat_path=NULL){


  
  if(! have_praat(praat_path)){
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
  praat_dsp_directory <- make_dsp_environment()

  formant_burg <- cs_wrap_praat_script(praat_location = get_praat(),
                                    script_code_to_run = readLines(praat_script),
                                    directory=praat_dsp_directory
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
    
    formantTabFile <- tempfile(fileext = ".csv", tmpdir = praat_dsp_directory)

    #Required for preventing errors in the hand off of file names containing spaces and () characters
    # to Praat
    soundFile <- tempfile(fileext = ".wav", tmpdir = praat_dsp_directory)
    R.utils::createLink(soundFile,origSoundFile)
    #Alternative route - much slower
    #file.copy(origSoundFile,soundFile)

    # sentence SoundFile /Users/frkkan96/Desktop/kaa_yw_pb.wav
    # real BeginTime 0.0
    # real EndTime 0.0
    # real Time_step 0.005
    # real Number_of_formants 5.0
    # real MaxHzFormant 5500.0
    # real WindowLength 0.025
    # real Pre_emphasis 50.0
    # boolean Track_formants 1
    # natural Number_of_tracks 3
    # natural Reference_F1 550
    # natural Reference_F2 1650
    # natural Reference_F3 2750
    # natural Reference_F4 3850
    # natural Reference_F5 4950
    # real Frequency_cost 1.0
    # real Bandwidth_cost 1.0
    # real Transition_cost 1.0 
    # word WindowShape Gaussian1
    # real RelativeWidth 1.0
    # word Spectrogram_window_shape Gaussian
    # real Spectrogram_resolution 40.0
    # sentence TrackOut /Users/frkkan96/Desktop/kaa_yw_pb.FormantTab
    
    outFormantTabFile <- formant_burg(soundFile,
                            beginTime,
                            endTime,
                            windowShift/1000, #Praat takes seconds
                            numFormants,
                            maxFormantHz,
                            windowSize/1000, #Praat takes seconds
                            preemphasis,
                            ifelse(trackFormants,1,0),
                            numberOfTracks,
                            nominalF1,
                            nominalF2,
                            nominalF3,
                            nominalF1*7,
                            nominalF1*9,
                            frequencyCost,
                            bandwidthCost,
                            transitionCost,
                            windowShape,
                            relativeWidth,
                            spectWindowShape,
                            spectResolution,
                            formantTabFile)
    
    inTable <- read.csv(file=outFormantTabFile
                               ,header=TRUE
                               ,na.strings =c("--undefined--","NA"),
                        sep = ",")
      

    # We need the sound file to extract some information
    origSound <- wrassp::read.AsspDataObj(soundFile)

    starTime = inTable[1,"time.s."]
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16", "INT16", "INT16")
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
    
    outDataObj = wrassp::addTrack(outDataObj, "F", as.matrix(fmTable), "INT16")

    bwTable <- inTable %>%
      dplyr::select(tidyselect::starts_with("B",ignore.case = FALSE)) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))

    names(bwTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "B", as.matrix(bwTable), "INT16")

    intTable <- inTable %>%
      dplyr::select(tidyselect::starts_with("L",ignore.case = FALSE)) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    names(bwTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "L", as.matrix(intTable), "INT16")
    
    
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
                               ncol = ncol(outDataObj$F))


      missing_bw_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$B))
      missing_int_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$L))
      # prepend values
      outDataObj$F = rbind(missing_fm_vals, outDataObj$F)
      outDataObj$B = rbind(missing_fm_vals, outDataObj$B)
      outDataObj$L = rbind(missing_int_vals, outDataObj$L)
      
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
  #Make sure that we clear the environment so that we may run the script again without complaints
  clear_dsp_environment(  praat_dsp_directory )
  
  logger::log_trace("Computed a formant track for the input signal file '",listOfFiles,"'.")
  if(!toFile){
    return(outDataObj)
  }

}

attr(praat_formant_burg,"ext") <-  c("pfm") 
attr(praat_formant_burg,"tracks") <-  c("F", "B" , "L")
attr(praat_formant_burg,"outputType") <-  c("SSFF")

#' Formant estimation using the FormantPath functionality of Praat
#' 
#' This function exposes Praat's functionality for iteratively searching for the best fit formant track for a file by adjusting the maximum formant frequency (frequency ceiling). In each iteration, Praat's built in function (burg algorithm). See \insertCite{Escudero.2009.10.1121/1.3180321}{superassp} for example of how the procedure has been used and \insertCite{Weenink2015}{superassp} for a description of how the optimal formant track is identified.
#'  If `stepsUpDown` is zero, then this function and the `praat_formant_burg` function would produce the same result if identical settings are used. 
#' The function also computes the intensity (L) of the best fit formant tracks based on the power of the spectrum at the frequency of the formant. Naturally, if the algorithm failed to find a formant in a specified time frame, then the function will not return a formant frequency, bandwidth and intensity estimation.
#' 
#' If the user only want to estimate formant frequencies that should later be manually corrected, computing them using the function [wrassp::forest] or even `praat_formant_burg` is _much_ quicker. The user should consider this function only if the use case specifically demands an iterative serch for best fit formants.
#'  
#'
#' @param listOfFiles a vector of wav file paths to be processed by function.
#' @param beginTime the time where processing should end (in s) The default is 0 (zero) which means that the computation of formants will start at the start of the sound file.
#' @param endTime the time where processing should end (in s) The default is 0 (zero) which means that formants will be computed up to the end of the file.
#' @param windowShift the analysis window shift length (in ms).
#' @param numFormants the number of formants that the analysis should try to find 
#' @param maxFormantHz The maximum frequency under which the formants should be found
#' 
#' @param windowSize the analysis window length (in ms).
#' @param preemphasis the frequency from which a preemphasis will be applied.
#' @param ceilingStepSize The function multiple searches for formant tracks with the frequency ceiling set to `maxhzformant*exp(-ceilingStepSize*stepsUpDown)` to `maxhzformant*exp(ceilingStepSize*stepsUpDown)`.
#' @param stepsUpDown The number of iterations of increases and decreases of the frequency ceiling to use when trying to find the an optimal formant track.
#' @param trackFormants boolean; Should Praat attempt to gather short time formant frequency estimates into tracks?
#' @param numberOfTracks The number of tracks to follow (if `trackFormants` is `TRUE`), and the number of tracks in the output. Information on frequencies bandwidths of formants with numbers above  `numberOfTracks` will be discarded.
#' @param nominalF1 Described by the Praat manual as the preferred value near which the first track wants to be. For average (i.e. adult female) speakers, this value will be around the average F1 for vowels of female speakers, i.e. 550 Hz.
#' @param nominalF2 Described by the Praat manual as the preferred value near which the second track wants to be. A good value will be around the average F2 for vowels of female speakers, i.e. 1650 Hz.
#' @param nominalF3 Described by the Praat manual as the preferred value near which the third track wants to be. A good value will be around the average F3 for vowels of female speakers, i.e. 2750 Hz. This argument will be ignored if you choose to have fewer than three tracks, i.e., if you are only interested in F1 and F2.
#' @param frequencyCost Described by the Praat manual as the preferred value near which the five track wants to be. In the unlikely case that you want five tracks, a good value may be around 4950 Hz.Frequency cost (per kiloHertz)
#' @param bandwidthCost Described by the Praat manual as the local cost of having a bandwidth, relative to the formant frequency. For instance, if a candidate has a formant frequency of 400 Hz and a bandwidth of 80 Hz, and Bandwidth cost is 1.0, the cost of having this formant in any track is (80/400) · 1.0 = 0.200. So we see that the procedure locally favours the inclusion of candidates with low relative bandwidths.
#' @param transitionCost Described by the Praat manual as the cost of having two different consecutive formant values in a track. For instance, if a proposed track through the candidates has two consecutive formant values of 300 Hz and 424 Hz, and Transition cost is 1.0/octave, the cost of having this large frequency jump is (0.5 octave) · (1.0/octave) = 0.500.
#' @param windowShape the analysis window function used when extracting part of a sound file for analysis. De faults to "Hanning".
#' @param relativeWidth the relative width of the windowing function used.
#' @param spectWindowShape The shape of the windowing function used for constructing the spectrogram. 
#' @param spectResolution The frequency resolution of the spectrogram from which formant intensities will be collected.
#' @param toFile write the output to a file? The file will be written in  `outputDirectory`, if defined, or in the same directory as the sound file. 
#' @param explicitExt the file extension that should be used.
#' @param outputDirectory set an explicit directory for where the signal file will be written. If not defined, the file will be written to the same directory as the sound file.
#' @param verbose Not implemented. Only included here for compatibility.  
#' @param praat_path give an explicit path for Praat. If the praat 
#'
#' @return An SSFF track data object (if `toFile=FALSE`) containing three fields ("F", "B" and "L") containing formant frequencies, bandwidth and intensities.
#' 
#' @export
#'
#'@seealso praat_formant_burg
#'
#' @references 
#'   \insertAllCited{}
#'   
#' @seealso [wrassp::forest]
#' @seealso [superassp::praat_formant_burg]


praat_formantpath_burg <- function(listOfFiles,
                               beginTime=0,
                               endTime=0,
                               windowShift=5.0,
                               numFormants=5.0,
                               maxFormantHz=5500.0,
                               windowSize=30,
                               preemphasis=50.0,
                               ceilingStepSize=0.05,
                               stepsUpDown=4,
                               trackFormants=TRUE,
                               numberOfTracks=3,
                               nominalF1=550,
                               nominalF2=1650,
                               nominalF3=2750,
                               frequencyCost=1.0,
                               bandwidthCost=1.0,
                               transitionCost=1.0,
                               windowShape="Gaussian1",
                               relativeWidth=1.0,
                               spectWindowShape="Gaussian",
                               spectResolution=40.0,
                               toFile=TRUE,
                               explicitExt="pfp",
                               outputDirectory=NULL,
                               verbose=FALSE,
                               praat_path=NULL){
  

  
  if(! have_praat(praat_path)){
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
                         file.path("inst","praat","formantpath_burg.praat"),
                         file.path(system.file(package = "superassp",mustWork = TRUE),"praat","formantpath_burg.praat")
  )
  
  praat_dsp_directory <- make_dsp_environment()

  formantpath_burg <- cs_wrap_praat_script(praat_location = get_praat(),
                                          script_code_to_run = readLines(praat_script),
                                          directory=praat_dsp_directory,
                                          return="last-argument")
  
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
    
    formantTabFile <- tempfile(fileext = ".csv",tmpdir = praat_dsp_directory)
    
    #Required for preventing errors in the handoff of file names containing spaces and () characters
    # to Praat
    soundFile <- tempfile(fileext = ".wav",tmpdir = praat_dsp_directory)
    R.utils::createLink(soundFile,origSoundFile)

    
    # real BeginTime 0.0
    # real EndTime 0.0
    # real Time_step 0.005
    # real Number_of_formants 5.0
    # real MaxHzFormant 5500.0
    # real WindowLength 0.025
    # real Pre_emphasis 50.0
    # real Ceiling_step_size 0.05
    # natural Number_of_steps_each_direction 4
    # boolean Track_formants 1
    # natural Number_of_tracks 3
    # natural Reference_F1 550
    # natural Reference_F2 1650
    # natural Reference_F3 2750
    # natural Reference_F4 3850
    # natural Reference_F5 4950
    # real Frequency_cost 1.0
    # real Bandwidth_cost 1.0
    # real Transition_cost 1.0 
    # word WindowShape Gaussian1
    # real RelativeWidth 1.0
    # word Spectrogram_window_shape Gaussian
    # real Spectrogram_resolution 40.0
    # sentence TrackOut /Users/frkkan96/Desktop/kaa_yw_pb.FormantTab
    
    outFormantTabFile <- formantpath_burg(soundFile,
                                      beginTime,
                                      endTime,
                                      windowShift/1000, #Praat takes seconds
                                      numFormants,
                                      maxFormantHz,
                                      windowSize/1000, #Praat takes seconds
                                      preemphasis,
                                      ceilingStepSize,
                                      stepsUpDown,
                                      ifelse(trackFormants,1,0),
                                      numberOfTracks,
                                      nominalF1,
                                      nominalF2,
                                      nominalF3,
                                      nominalF1*7,
                                      nominalF1*9,
                                      frequencyCost,
                                      bandwidthCost,
                                      transitionCost,
                                      windowShape,
                                      relativeWidth,
                                      spectWindowShape,
                                      spectResolution,
                                      formantTabFile)
    
    inTable <- read.csv(file=outFormantTabFile
                        ,header=TRUE
                        ,na.strings =c("--undefined--","NA"),
                        sep = ",")
    
    
    
    # We need the sound file to extract some information
    origSound <- wrassp::read.AsspDataObj(soundFile)
    
    starTime = inTable[1,"time.s."]
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16", "INT16", "INT16")
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
    
    outDataObj = wrassp::addTrack(outDataObj, "F", as.matrix(fmTable), "INT16")
    
    bwTable <- inTable %>%
      dplyr::select(tidyselect::starts_with("B",ignore.case = FALSE)) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    names(bwTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "B", as.matrix(bwTable), "INT16")
    
    intTable <- inTable %>%
      dplyr::select(tidyselect::starts_with("L",ignore.case = FALSE)) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    names(bwTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "L", as.matrix(intTable), "INT16")
    
    
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
                               ncol = ncol(outDataObj$F))
      
      
      missing_bw_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$B))
      missing_int_vals = matrix(0,
                                nrow = nr_of_missing_samples,
                                ncol = ncol(outDataObj$L))
      # prepend values
      outDataObj$F = rbind(missing_fm_vals, outDataObj$F)
      outDataObj$B = rbind(missing_fm_vals, outDataObj$B)
      outDataObj$L = rbind(missing_int_vals, outDataObj$L)
      
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
  
  #Make sure that we clear the environment so that we may run the script again without complaints

  clear_dsp_environment(  praat_dsp_directory )

  logger::log_trace("Computed a formantpath track for the input signal file '",listOfFiles,"'.")
  
  if(!toFile){
    return(outDataObj)
  }
  
  
}

attr(praat_formantpath_burg,"ext") <-  c("pfp") 
attr(praat_formantpath_burg,"tracks") <-  c("F", "B" , "L")
attr(praat_formantpath_burg,"outputType") <-  c("SSFF")


#' Call the 'praat_sauce' analysis bundle to generate SSFF tracks 
#'
#' This function applies the \code{praat_sauce} bundle of Praat scripts
#' on a single file and puts the resulting data as tracks in an SSFF file.
#' By default, all analyses are applied to a windowed portion of the signal
#' every 5ms.  
#'
#' @param listOfFiles A vector of file names
#' @param beginTime The time in the sound file where analysis should start.
#' @param endTime The last time point to be included in the analysed sample. If zero (0), the sound file will be included until the end.
#' 
#' @param windowShift The time shift to the next analysis window. Defaults to every 5ms.
#' @param windowSize The analysis window length (in ms).
#' @param minF The minimal f0 to search for.
#' @param maxF The maximum f0 to search for.
#' @param formantTracking Boolean; Should the formant tracking abilities of Praat be used? Defaults to TRUE. If disabled, the raw window-by-window formant values will be used. 
#' @param numFormants The number of formants to be found within the frequency space.
#' @param maxFormantHz The cutoff frequency used when finding the `numFormants` formants.
#' @param nominalF1 The nominal F1 used in formant tracking.
#' @param nominalF2 The nominal F2 used in formant tracking.
#' @param nominalF3 The nominal F3 used in formant tracking.
#' @param preEmphFrom The frequency from which pre-emphasis will be applied.
#' @param useBandwidthFormula Should the bandwidth calculation metod of \insertCite{Hawks.1995.10.1121/1.412986;textual}{superassp} be used, instead of Praat's internal algorithm. Defaults to TRUE (use Hawks & Miller's method).
#' @param channel Which channel to use analyse. Defaults to the first channel.
#' @param resample_to_16k Resample the signal to 16000 Hz before processing? Defaults to TRUE.
#' @param toFile Should the SSFF track file be written to disk? Defaults to true, which by default means that it will be written as a file with a `explicitExt` extension  next to the sound file.
#' @param explicitExt The default file extension to use for the SSFF file.
#' @param outputDirectory A path to an alternative output directory.
#' @param verbose Verbose output. Not currently used.
#' @param praat_path=NULL
#'
#' @return
#' 
#' This function builds an SSFF track object and writes it to disk, or returns it (`toFile==FALSE`). The track object will
#' contain tracks with these fields:
#' \describe{
#'   \item{f0}{A track of the fundamental frequency computed by Praat from the signal as part of the PraatSauce analysis.}
#'   \item{fm}{Computed formant frequency measures, optionally tracked, one column per formant. }
#'   \item{bw}{Computed formant bandwiths, one column per formant.}
#'   \item{H}{The amplitudes of the first three harmonics, computed without accounting for the influence of adjecent formants (uncorrected).}
#'   \item{Hc}{The amplitudes of the first three harmonics, in which the influence of adjecent formants has been accounted for (corrected).}
#'   \item{A}{The amplitudes of the three harmonics that are closest to the center frequency of the first three formants, computed without accounting for the influence of the adjecent formants (uncorrected).}
#'   \item{Ac}{The amplitudes of the three harmonics that are closest to the center frequency of the first three formants, in which the influence of the adjecent formants has been accounted for (corrected).}
#'   \item{HH}{The  differences between the (uncorrected) amplitudes of the first and second (column 1) and second and fourth (column 2).}
#'   \item{HHc}{The  differences between the amplitudes of the first and second (column 1) and second and fourth (column 2), computed with correction with regards to neighbouring formants.}
#'   \item{HA}{The  differences between the (uncorrected) amplitudes of the harmonics closest to the center frequencies of the first and second (column 1) and second and fourth (column 2) formants.}
#'   \item{HAc}{The  differences between the amplitudes of the harmonics closest to the center frequencies of the first and second (column 1) and second and fourth (column 2) formants, corrected for the influence of the adjecent formants themselves.}
#'   \item{cpp}{A track of containing the smoothed Cepstral Peak Prominence \insertCite{Fraile.2014.10.1016/j.bspc.2014.07.001,Hillenbrand.1994.10.1044/jshr.3704.769}{superassp} across the acoustic signal.}
#'   \item{hnr}{Harmonic-to-noise ratios as defined by \insertCite{Krom.1993.10.1044/jshr.3602.254}{superassp} for frequencies 0-500 Hz, 0-1500 Hz, 0-2500 Hz, and 0-3500 Hz, respectively. (four columns)}
#' }
#' 
#' @export
#' 
#' @references 
#'   \insertAllCited{}
#'
#' @importFrom dplyr %>%
#' @examples
#' \dontrun{
#' 
#' }

praat_sauce <- function(listOfFiles,
                        beginTime=NULL,
                        endTime=NULL,
                        windowShift=5.0,
                        windowSize=25,
                        minF=50,
                        maxF=300,
                        formantTracking=TRUE,
                        numFormants=5,
                        maxFormantHz=5000,
                        nominalF1=500,
                        nominalF2=1500,
                        nominalF3=2500,
                        preEmphFrom=50,
                        useBandwidthFormula=FALSE,
                        channel=1,
                        resample_to_16k=TRUE,
                        toFile=TRUE,
                        explicitExt="psa",
                        outputDirectory=NULL,
                        verbose=FALSE,
                        praat_path=NULL){


  if(! have_praat(praat_path)){
    stop("Could not find praat. Please specify a full path.")
  }
  
  if(length(listOfFiles) > 1 & ! toFile){
    stop("length(listOfFiles) is > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.")
  }
  
  tryCatch({
    fileBeginEnd <- data.frame(
      listOfFiles = listOfFiles, 
      beginTime = ifelse(is.null(beginTime),0, beginTime),
      endTime=ifelse(is.null(endTime),0, endTime)
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
  praat_dsp_directory <- make_dsp_environment()

  praatsauce <- cs_wrap_praat_script(praat_location = get_praat(),
                                           script_code_to_run = readLines(praat_script),
                                           directory=praat_dsp_directory,
                                           return="last-argument")
  
  #Copy additional files
  copied <- file.copy(additional_scripts,praat_dsp_directory,overwrite = TRUE)

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
    
    outputfile <- tempfile(fileext = ".csv",tmpdir = praat_dsp_directory)
    
    #Required for preventing errors in the handoff of file names containing spaces and () characters
    # to Praat
    soundFile <- tempfile(fileext = ".wav",tmpdir = praat_dsp_directory)
    R.utils::createLink(soundFile,origSoundFile)

# 
#     form File and measures
#     sentence inputfile /Users/frkkan96/Desktop/a2.wav
#     real beginTime 0
#     real endTime 0
#     integer channel 1
#     comment 1= n equidistant points
#     comment 2=   option every n milliseconds
#     natural measure 2
#     natural Points 5
#     boolean resample_to_16k 1    
#     boolean pitchTracking 1
#     boolean formantMeasures 1
#     boolean spectralMeasures 1
#     positive windowLength 0.025
#     positive windowPosition 0.5
#     positive maxFormantHz 5000
#     positive spectrogramWindow 0.005
#     boolean useExistingPitch 0 
#     positive f0min 50
#     positive f0max 300
#     real timeStep 0
#     integer maxNumFormants 5
#     positive preEmphFrom 50
#     boolean formantTracking 1
#     positive F1ref 500
#     positive F2ref 1500
#     positive F3ref 2500
#     boolean useExistingFormants 0
#     boolean useBandwidthFormula 0
#     sentence outputfile /Users/frkkan96/Desktop/spectral_measures.txt
#     endform
        
    
    outputfile <- praatsauce(soundFile,
                            ifelse(is.null(beginTime),0,beginTime),
                            ifelse(is.null(endTime),0,endTime),
                            channel,
                            2,
                            windowShift,
                            ifelse(resample_to_16k,1,0),
                            1,
                            1,
                            1,
                            windowSize/1000, #Praat takes seconds.
                            0.5,
                            maxFormantHz,
                            0.005,
                            0, #ifelse(useExistingPitch,1,0), TO IMPLEMENT LATER
                            minF,
                            maxF,
                            windowShift/1000, #Praat takes seconds.
                            numFormants,
                            preEmphFrom,
                            ifelse(formantTracking,1,0),
                            nominalF1,
                            nominalF2,
                            nominalF3,
                            0, #ifelse(useExistingFormants,1,0), TO IMPLEMENT LATER
                            ifelse(useBandwidthFormula,1,0),
                            outputfile)
    
    inTable <- read.csv(file=outputfile
                        ,header=TRUE
                        ,na.strings =c("--undefined--","NA"),
                        sep = ",")
    
    

    #### Create the SSFF object     ####

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
    
    

    #### f0 values are placed in the track "f0"     #### 

    # Too elaborate for a single column manipulation, but it is consistent with subsequent processing of 
    # multiple columns
    f0Table <- inTable %>%
      dplyr::select(f0) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(f0=as.integer(f0))
    
    names(f0Table) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "f0", as.matrix(f0Table), "INT16")
    
      
  
    
    ### 
    # This properties are only extracted if the user has specified that formants should be computed,
    # OR spectral measures are to be computed
    ### 
    #if(formantMeasures){
        

    #### Formant frequencies are placed in the track "fm"    #### 
    
    fmTable <- inTable %>%
      dplyr::select(tidyselect::starts_with("F",ignore.case = FALSE)) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
  
    
    names(fmTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "fm", as.matrix(fmTable), "INT16")
    

    #### Formant bandwidths are placed in the track "bw"    #### 
    
    bwTable <- inTable %>%
      dplyr::select(tidyselect::starts_with("B",ignore.case = FALSE)) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    names(bwTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "bw", as.matrix(bwTable), "INT16")
    
    #}
    
    
    ### 
    # Extraction of the more special voice characteristics
    ###    
    
    
    #if(spectralMeasures){
      

    #### Amplitudes (uncorrected) of harmonics  are placed in the track "H"    #### 
    
    harmTable <- inTable %>%
      dplyr::select(H1u:H4u) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    names(harmTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "H", as.matrix(harmTable), "INT16")


    #### Corrected amplitudes  of harmonics  are placed in the track "Hc"    #### 

    
    harmTable <- inTable %>%
      dplyr::select(H1c:H4c) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    names(harmTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "Hc", as.matrix(harmTable), "INT16")
    
    
    #### The (uncorrected) amplitudes of harmonics  closest to F1-F3 are placed in the track "A"    #### 
    
    harmTable <- inTable %>%
      dplyr::select(A1u:A3u) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    names(harmTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "A", as.matrix(harmTable), "INT16")
    
    #### The corrected amplitudes of harmonics closest to F1-F3 are placed in the track "Ac"    #### 
    
    harmTable <- inTable %>%
      dplyr::select(A1c:A3c) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    names(harmTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "Ac", as.matrix(harmTable), "INT16")
    
    #### The (uncorrected) first and second harmonics closest to 2 and 5k Hz respectively are placed in the columns in "H25K"    #### 

    harmTable <- inTable %>%
      dplyr::select(H2Ku,H5Ku) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    names(harmTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "H25K", as.matrix(harmTable), "INT16")

    

    #### The  differences between the (uncorrected) amplitudes of the first and second harmonics and the second and fourth are stored in the "HH" field    #### 

    harmTable <- inTable %>%
      dplyr::select(H1H2u,H2H4u) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    names(harmTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "HH", as.matrix(harmTable), "INT16")
    
    #### The differences between the corrected amplitudes of the first and second harmonics and the second and fourth are stored in the "HHc" field    #### 
    
    harmTable <- inTable %>%
      dplyr::select(H1H2c,H2H4c) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    names(harmTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "HHc", as.matrix(harmTable), "INT16")

    #### The differences between the (uncorrected) amplitudes of the first harmonic and the harmonics closest to F1,F2, and F3 are placed as columns in the "HA" field    #### 
    
    harmTable <- inTable %>%
      dplyr::select(H1A1u,H1A2u,H1A3u) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    names(harmTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "HA", as.matrix(harmTable), "INT16")
    
    #### The  differences between the corrected amplitudes of the first harmonic and the harmonics closest to F1,F2, and F3 are placed as columns in the "HAc" field    #### 
    
    harmTable <- inTable %>%
      dplyr::select(H1A1c,H1A2c,H1A3c) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    names(harmTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "HAc", as.matrix(harmTable), "INT16")

    #### The cepstral peak prominence is inserted into the "cpp" field    #### 
    
    cppTable <- inTable %>%
      dplyr::select(CPP) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    names(cppTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "cpp", as.matrix(cppTable), "INT16")
    
    #### The harmonic-to-noise ratio measured from 0 to 500, 1500, 2500 and 3500 Hz respectively is inserted into columns of the field "hnr"    #### 
    
    hnrTable <- inTable %>%
      dplyr::select(HNR05,HNR15,HNR25,HNR35) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    names(hnrTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "hnr", as.matrix(hnrTable), "INT16")     

    #}
    
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
  
  #Make sure that we clear the environment so that we may run the script again without complaints
  clear_dsp_environment(  praat_dsp_directory )
  
  logger::log_trace("Computed a praat_sauce track for the input signal file '",listOfFiles,"'.")
  if(!toFile){
    return(outDataObj)
  }
  
}
attr(praat_sauce,"ext") <-  c("psa") 
attr(praat_sauce,"tracks") <-  c("f0","fm", "bw","H","Hc","A","Ac","H25K","HH","HHc","HA","HAc","cpp","hnr")
attr(praat_sauce,"outputType") <-  c("SSFF")

#' Compute a sound signal intensity track using Praat
#'
#' This function calls Praat to compute an intensity contour for a sound file,
#' and return the result as a Assp Data Object that is stored in an SSFF file
#' (the default), or returned as an object. The analysis will attempt to
#' compensate for the effect of the periodicity introduced by pitch of the
#' produced speech, as well as the constant pressure in the recording
#' environment if supplied with appropriate arguments.
#'
#' @param listOfFiles A list of sound files that should be analyzed.
#' @param beginTime The time (in s) where the analysis should begin.
#' @param endTime The time (in s) where the analysis should end
#' @param windowShift The time step (in ms) of the resulting intensity contour.
#' @param minF The minimum periodicity of the sound signal. If set too
#'   high, the intensity variations within a pitch period will influence the
#'   computed intensity contour. If set too low, the smearing of the intensity
#'   contour may hide rapid intensity variations.
#' @param subtractMean Should the average intensity be subtracted in order to
#'   compensate for the constant pressure of the recording environment?
#' @param window Which windowing function should be applied when extracting part
#'   of the recording for analysis? Allowed windowing functions are
#'   "rectangular", "triangular", "parabolic", "Hanning", "Hamming",
#'   "Gaussian1", "Gaussian2", "Gaussian3", "Gaussian4", "Gaussian5", "Kaiser1",
#'   and "Kaiser2". Consult the Praat manual for more details.
#' @param relativeWidth The relative with of the window used when extracting
#'   part of a sound file for analysis.
#' @param toFile Should the SSFF track be written to file (TRUE) or returned as
#'   an object (FALSE)
#' @param explicitExt The SSFF file written to disk will have the same name as
#'   the original sound file, but with this file extension.
#' @param outputDirectory The directory where the SSFF track will be stored.
#'   Defaults to the same directory as the sound file.
#' @param verbose For comparability with wrassp functions, and expected by EmuR. Nothing happens if you set it to FALSE:
#' @param praat_path An explicit path to the Praat executable.
#'
#' @return An SSFF object containing the intensity track (if `toFile==FALSE`).
#'   
praat_intensity <- function(listOfFiles,beginTime=0,endTime=0,windowShift=5.0,minF=80,subtractMean=TRUE,window="Gaussian1",relativeWidth=1.0,toFile=TRUE,explicitExt="int",outputDirectory=NULL,verbose=TRUE,praat_path=NULL){
  
  # real BeginTime 0.0
  # real EndTime 0.0
  # real Time_step 0.0
  # real Minimal_f0_Frequency 50.0
  # boolean Subtract_mean 1
  # word WindowShape Gaussian1
  # real RelativeWidth 1.0
  
  if(! have_praat(praat_path)){
    stop("Could not find praat. Please specify a full path.")
  }
  
  if(length(listOfFiles) > 1 & ! toFile){
    stop("length(listOfFiles) is > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.")
  }
  
  wt <- c("rectangular", "triangular", "parabolic", "Hanning", "Hamming","Gaussian1", "Gaussian2", "Gaussian3", "Gaussian4", "Gaussian5", "Kaiser1","Kaiser2")
  
  if (! tolower(window) %in% tolower(wt)) {
    stop("The supplied window function is not recognised. \n",
    "The window needs to be one of: ", paste(wt, collapse = ", "),
    "\nPlease consult the Praat manual (Intensity) for more information.")
  }
  
  tryCatch({
    fileBeginEnd <- data.frame(
      listOfFiles = listOfFiles, 
      beginTime = beginTime,
      endTime=endTime
    )
  },error=function(e){stop("The beginTime and endTime must either be a single value or the same length as listOfFiles")})
  
  
  
  praat_script <- ifelse(PRAAT_DEVEL== TRUE,
                         file.path("inst","praat","intensity.praat"),
                         file.path(system.file(package = "superassp",mustWork = TRUE),"praat","intensity.praat")
  )
  
  praat_dsp_directory <- make_dsp_environment()
  
  praat_intensity <- cs_wrap_praat_script(praat_location = get_praat(),
                                          script_code_to_run = readLines(praat_script),
                                          directory=praat_dsp_directory,
                                          return="last-argument")
  
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
    
    intensityTabFile <- tempfile(fileext = ".csv",tmpdir = praat_dsp_directory)
    
    #Required for preventing errors in the handoff of file names containing spaces and () characters
    # to Praat
    soundFile <- tempfile(fileext = ".wav",tmpdir = praat_dsp_directory)
    R.utils::createLink(soundFile,origSoundFile)

    
    outIntensityTabFile <- praat_intensity(soundFile,
                                      beginTime,
                                      endTime,
                                      windowShift/1000, #Praat expects shift in s
                                      minF,
                                      ifelse(subtractMean,1,0),
                                      window,
                                      relativeWidth,
                                      intensityTabFile)
    
    inTable <- read.csv(file=outIntensityTabFile
                        ,header=TRUE
                        ,na.strings =c("--undefined--","NA"),
                        sep = ",")
    

    
    # We need the sound file to extract some information
    origSound <- wrassp::read.AsspDataObj(soundFile)
    
    starTime = inTable[1,"Time..s."]
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16")
    #Use the time separation between second and first formant measurement time stamps to compute a sample frequency.
    sampleRate <-  as.numeric(1 / (inTable[2,"Time..s."] - inTable[1,"Time..s."]))
    attr(outDataObj, "sampleRate") <- sampleRate
    
    attr(outDataObj, "origFreq") <-  as.numeric(attr(origSound, "sampleRate"))
    startTime <- as.numeric(inTable[1,"Time..s."])
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
    class(outDataObj) = "AsspDataObj"
    
    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    intensityTable <- inTable %>%
      dplyr::select(`Intensity..dB.`) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
    
    noValues <- nrow(intensityTable)
    noCols <- 1
    
    names(intensityTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "int", as.matrix(intensityTable), "INT16")
    
    
    ## Apply fix from Emu-SDMS manual
    ##https://raw.githubusercontent.com/IPS-LMU/The-EMU-SDMS-Manual/master/R/praatToFormants2AsspDataObj.R
    
    # add missing values at the start as Praat sometimes
    # has very late start values which causes issues
    # in the SSFF file format as this sets the startRecord
    # depending on the start time of the first sample
    if( startTime > (1/sampleRate) ){
      
      nr_of_missing_samples = as.integer(floor(startTime / (1/sampleRate)))
      
      missing_int_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$int))
      
      
      
      # prepend values
      outDataObj$int = rbind(missing_int_vals, outDataObj$int)

      
      # fix start time
      attr(outDataObj, "startTime") = startTime - nr_of_missing_samples * (1/sampleRate)
    }
    
    assertthat::assert_that(wrassp::is.AsspDataObj(outDataObj),
                            msg = paste("The AsspDataObj created by the praat_intensity function is invalid.\nPlease check the table file '",tabfile,"' for errors.",sep=""))
    
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
  
  #Make sure that we clear the environment so that we may run the script again without complaints
  clear_dsp_environment(  praat_dsp_directory )
  logger::log_trace("Computed an intensity track for the input signal file '",listOfFiles,"'.")
  if(!toFile){
    return(outDataObj)
  }
  
}

attr(praat_intensity,"ext") <-  c("int") 
attr(praat_intensity,"tracks") <-  c("int")
attr(praat_intensity,"outputType") <-  c("SSFF")




#' Compute spectral moments using Praat
#'
#' This function takes a sound file, or an indicated part of a sound file, and
#' computes the first (Center of gravity), second (Standard deviation), third
#' (skewness), and fourth (kurtosis) spectral moments at regular intervals. The
#' moments are stored as separate tracks in an SSFF signal object on disk, or
#' optionally returned as an object (toFile=FALSE).
#' 
#' By default, the moments are based on the entire spectrum of the sound file, which means that they will depend on the sampling frequency of the file. The user may however choose to include just a portion of the spectrum up to a cutoff frequency if required. This way, the user may opt to limit the analysis to half of the smallest sampling rate for a collection of sound files and thereby produce comparable results.
#'
#'
#' @param listOfFiles The full paths of the files to be processed.
#' @param beginTime The start time of the portion of a wave that should be included. This argument needs to be either be a single value or a vector of the same length as `listOfFiles`.
#' @param endTime The end time of the portion of a wave that should be included. Like `beginTime` this argument needs to be either be a single value or a vector of the same length as `listOfFiles`.
#' @param windowShift The time step between (time) analysis windows (in ms). 
#' @param windowSize The size of the time aligned analysis window (in ms). 
#' @param freqBinSize The spectral resolution.
#' @param power The power to be used when computing the spectral moments. If `power=1` the spectral moments will be computed from the absolute spectrum, and if `power=2` they will be computed based on the power spectrum. 
#' @param maximumFrequency The cutoff frequency (in Hz) used when computing the spectrum. Frequencies above this cutoff will not be included when computing spectral moments.
#' @param windowShape The window type used for extracting a section of the wave file for analysis. Permitted values are "rectangular", "triangular", "parabolic", "Hanning", "Hamming", "Gaussian1", "Gaussian2", "Gaussian3", "Gaussian4", "Gaussian5", "Kaiser1", and "Kaiser2". See the Praat manual for a descriptio of these window shapes.
#' @param relativeWidth The relative width of the windowing function used for extracting part of the sound file.
#' @param toFile Should the SSFF signal tracks be store on disk or returned?
#' @param explicitExt The signal file extension.
#' @param outputDirectory Where should the signal file be stored. If `NULL`, the signal file will be stored in the same directory as the wave file.
#' @param verbose Produce verbose output?
#' @param praat_path The location where the Praat executable is stored.
#'
#' @return An SSFF file (if toFile=FALSE), or nothing.
#' @export
#'
praat_moments <- function(listOfFiles,
                        beginTime=NULL,
                        endTime=NULL,
                        windowShift=5.0,
                        windowSize=25,
                        freqBinSize=20,
                        power=2,
                        maximumFrequency=NULL,
                        windowShape="Gaussian1",
                        relativeWidth=1.0,
                        toFile=TRUE,
                        explicitExt="pmo",
                        outputDirectory=NULL,
                        verbose=FALSE,
                        praat_path=NULL){
  
  
  if(! have_praat(praat_path)){
    stop("Could not find praat. Please specify a full path.")
  }
  
  if(length(listOfFiles) > 1 & ! toFile){
    stop("length(listOfFiles) is > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.")
  }
  
  tryCatch({
    fileBeginEnd <- data.frame(
      listOfFiles = listOfFiles, 
      beginTime = ifelse(is.null(beginTime),0, beginTime),
      endTime=ifelse(is.null(endTime),0, endTime)
    )
  },error=function(e){stop("The beginTime and endTime must either be a single value or the same length as listOfFiles")})
  
  
  
  praat_script <- ifelse(PRAAT_DEVEL== TRUE,
                         file.path("inst","praat","praat_spectral_moments.praat"),
                         file.path(system.file(package = "superassp",mustWork = TRUE),"praat","praat_spectral_moments.praat")
  )
  praat_dsp_directory <- make_dsp_environment()

  pmoments <- cs_wrap_praat_script(praat_location = get_praat(),
                                  script_code_to_run = readLines(praat_script),
                                  directory = praat_dsp_directory,
                                  return="last-argument")

  #Check that all files exists before we begin
  filesEx <- file.exists(listOfFiles)
  if(!all(filesEx) ){
    filedNotExists <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ",paste(filedNotExists, collapse = ", "))
  }

  
  #The empty vector of file names that should be returned
  outListOfFiles <- c()
  
  for(i in 1:nrow(fileBeginEnd)){ 
    origSoundFile <- fileBeginEnd[i, "listOfFiles"]
    
    beginTime <- fileBeginEnd[i, "beginTime"]
    endTime <- fileBeginEnd[i, "endTime"]
    
    outputfile <- tempfile(fileext = ".csv",tmpdir = praat_dsp_directory)
    
    #Required for preventing errors in the handoff of file names containing spaces and () characters
    # to Praat
    soundFile <- tempfile(fileext = ".wav",tmpdir = praat_dsp_directory)
    R.utils::createLink(soundFile,origSoundFile)
    
    # form Compute the spectral moments 1-4 
    # sentence SoundFile /Users/frkkan96/Desktop/a1.wav
    # real BeginTime 0.0
    # real EndTime 0.0
    # real WindowLength 0.005
    # real Maximum_frequency_(Hz) 0.0
    # real Time_step 0.005
    # real Frequency_step 20.0
    # real Power 2
    # word WindowShape Gaussian1
    # real RelativeWidth 1.0
    # sentence TrackOut /Users/frkkan96/Desktop/a1.mom
    # endform
    
    outputfile <- pmoments(soundFile,
                            ifelse(is.null(beginTime),0,beginTime),
                            ifelse(is.null(endTime),0,endTime),
                           windowSize/1000, #Praat takes seconds.
                           ifelse(is.null(maximumFrequency),0.0,maximumFrequency),
                           windowShift/1000, #Praat takes seconds.
                           freqBinSize/1000,
                           power,
                           windowShape,
                           1.0,
                           outputfile)
    
    inTable <- read.csv(file=outputfile
                        ,header=TRUE
                        ,na.strings =c("--undefined--","NA"),
                        sep = ",")
    
    
    #return(inTable)
    #### Create the SSFF object     ####
    
    # We need the sound file to extract some information
    origSound <- wrassp::read.AsspDataObj(soundFile)
    
    startTime = inTable[1,"Time"]
    
    outDataObj = list()
    attr(outDataObj, "trackFormats") <- rep("INT16",ncol(inTable)-1) #All but the "t" column
    #Use the time separation between second and first formant measurement time stamps to compute a sample frequency.
    sampleRate <-  as.numeric(1 / (inTable[2,"Time"] - inTable[1,"Time"]))
    attr(outDataObj, "sampleRate") <- sampleRate
    
    attr(outDataObj, "origFreq") <-  as.numeric(attr(origSound, "sampleRate"))
    startTime <- as.numeric(inTable[1,"Time"])
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
    class(outDataObj) = "AsspDataObj"
    
    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    
    
    #### First spectral moment : Center of gravity     #### 
    
    cogTable <- inTable %>%
      dplyr::select(CenterOfGravity) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(CenterOfGravity=as.integer(CenterOfGravity))
    
    names(cogTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "cog", as.matrix(cogTable), "INT16")

    #### Second spectral moment : Standard deviation     #### 
    
    sdTable <- inTable %>%
      dplyr::select(SD) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(SD=as.integer(SD))
    
    names(sdTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "sd", as.matrix(sdTable), "INT16")

    #### Third spectral moment : Skewness     #### 
    
    skewTable <- inTable %>%
      dplyr::select(Skewness) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(Skewness=as.integer(Skewness))
    
    names(sdTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "skew", as.matrix(skewTable), "INT16")
    
    #### Fourth spectral moment : Kurtosis     #### 
    
    kurtTable <- inTable %>%
      dplyr::select(Kurtosis) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(Kurtosis=as.integer(Kurtosis))
    
    names(sdTable) <- NULL
    
    outDataObj = wrassp::addTrack(outDataObj, "kurt", as.matrix(kurtTable), "INT16")
    
    
    
    ## Apply fix from Emu-SDMS manual
    ##https://raw.githubusercontent.com/IPS-LMU/The-EMU-SDMS-Manual/master/R/praatToFormants2AsspDataObj.R
    
    # add missing values at the start as Praat sometimes
    # has very late start values which causes issues
    # in the SSFF file format as this sets the startRecord
    # depending on the start time of the first sample
    if( startTime > (1/sampleRate) ){
      
      nr_of_missing_samples = as.integer(floor(startTime / (1/sampleRate)))
      
      missing_cog_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$cog))
      
      missing_sd_vals = matrix(0,
                                nrow = nr_of_missing_samples,
                                ncol = ncol(outDataObj$sd))
      missing_skew_vals = matrix(0,
                                nrow = nr_of_missing_samples,
                                ncol = ncol(outDataObj$skew))
      missing_kurt_vals = matrix(0,
                                nrow = nr_of_missing_samples,
                                ncol = ncol(outDataObj$kurt))
      # prepend values
      outDataObj$cog = rbind(missing_cog_vals, outDataObj$cog)
      outDataObj$sd = rbind(missing_sd_vals, outDataObj$sd)
      outDataObj$skew = rbind(missing_skew_vals, outDataObj$skew)
      outDataObj$kurt = rbind(missing_kurt_vals, outDataObj$kurt)

      
      # fix start time
      attr(outDataObj, "startTime") = startTime - nr_of_missing_samples * (1/sampleRate)
    }
    
    assertthat::assert_that(wrassp::is.AsspDataObj(outDataObj),
                            msg = paste("The AsspDataObj created by the praat_moments function is invalid.\nPlease check the table file '",tabfile,"' for errors.",sep=""))
    
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
  
  #Make sure that we clear the environment so that we may run the script again without complaints
  clear_dsp_environment(  praat_dsp_directory )
  logger::log_trace("Computed spectral moments tracks for the input signal file '",listOfFiles,"'.")
  if(!toFile){
    return(outDataObj)
  }
  
  
}
attr(praat_moments,"ext") <-  c("pmo") 
attr(praat_moments,"tracks") <-  c("cog","sd","skew","kurt")


#' Compute f0 tracks using Praat 
#' 
#' This function calls Praat to compute f0 tracks. Both the auto-correlation and cross-correlation methods are used, and the results are stored in separate fields in the returned SSFF track object. 
#' Most arguments to the function map directly to formal arguments to the underlying Praat procedure, and the description of these are therefore replicated here. See the Praat manual for more information.
#'
#'
#' @param listOfFiles A vector of file paths to wav files.
#' @param beginTime The start time of the section of the sound file that should be processed.
#' @param endTime The end time of the section of the sound file that should be processed.
#' @param windowShift  The measurement interval (frame duration), in seconds. If you supply 0, Praat will use a time step of 0.75 / (pitch floor), e.g. 0.01 seconds if the pitch floor is 75 Hz; in this example, Praat computes 100 pitch values per second.
#' @param minF candidates below this frequency will not be recruited. This parameter determines the effective length of the analysis window: it will be 3 longest periods long, i.e., if the pitch floor is 75 Hz, the window will be effectively 3/75 = 0.04 seconds long. Note that if you set the time step to zero, the analysis windows for consecutive measurements will overlap appreciably: Praat will always compute 4 pitch values within one window length, i.e., the degree of oversampling is 4.
#' @param maxF Candidates above this frequency will be ignored.
#' @param very.accurate If FALSE, the window is a Hanning window with a physical length of 3 / (pitch floor). If TRUE, the window is a Gaussian window with a physical length of 6 / (pitch floor), i.e. twice the effective length.
#' @param max.f0.candidates The maximum numbrf of f0 candidates to consider
#' @param silence.threshold Frames that do not contain amplitudes above this threshold (relative to the global maximum amplitude), are probably silent.
#' @param voicing.threshold The strength of the unvoiced candidate, relative to the maximum possible autocorrelation. To increase the number of unvoiced decisions, increase this value.
#' @param octave.cost The degree of favoring of high-frequency candidates, relative to the maximum possible autocorrelation. This is necessary because even (or: especially) in the case of a perfectly periodic signal, all undertones of f0 are equally strong candidates as f0 itself. To more strongly favour recruitment of high-frequency candidates, increase this value.
#' @param octave.jump.cost Degree of disfavoring of pitch changes, relative to the maximum possible autocorrelation. To decrease the number of large frequency jumps, increase this value. In contrast with what is described by \insertCite{Boersma1993}{superassp}, this value will be corrected for the time step: multiply by 10ms / windowShift to get the value in the way it is used in the formulas in the article.
#' @param voiced.voiceless.cost Degree of disfavoring of voiced/unvoiced transitions, relative to the maximum possible autocorrelation. To decrease the number of voiced/unvoiced transitions, increase this value. In contrast with what is described in the article, this value will be corrected for the time step: multiply by 10 ms / windowShift to get the value in the way it is used in the formulas in \insertCite{Boersma1993}{superassp}.
#' @param corr.only boolean; Compute autocorrelation (AC) and cross-correlation (CC) estimates of f0 only. If `FALSE` (the default) the function will additionally estimate f0 using a Spatial Pitch Network (SPINET) model \insertCite{Cohen.1995.10.1121/1.413512}{superassp} as well as as using a spectral compression (SHS) model \insertCite{Hermes.10.1121/1.396427}{superassp}. The computational load is increased considerably by these f0 estimates, and should be avoided if not explicitly needed by setting this parameter to TRUE.
#' @param windowSize the window size used for computing the SPINET model. 
#' @param min.filter.freq the minimum filter frequency used when computing the SPINET model.
#' @param max.filter.freq the maximum filter frequency used when computing the SPINET model.
#' @param filters the number of filters used when computing the SPINET model.
#' @param max.freq.components higher frequencies will not be considered when computing SHS.
#' @param subharmonics the maximum number of harmonics that add up to the pitch in SHS.
#' @param compression the factor by which successive compressed spectra are multiplied before the summation in SHS.
#' @param points.per.octave determines the sampling of the logarithmic frequency scale in SHS.
#' @inheritParams praat_formant_burg
#'
#' @references 
#'   \insertAllCited{}
#' @return An SSFF object containing the f0 tracks (if `toFile==FALSE`).
#' 
#' @export
#'
#'
#'


praat_pitch <- function(listOfFiles,
                        beginTime=0,
                        endTime=0,
                        windowShift=5.0,
                        minF=75,
                        maxF=600,
                        max.f0.candidates=15,
                        very.accurate=TRUE,
                        silence.threshold=0.03,
                        voicing.threshold=0.45,
                        octave.cost=0.01,
                        octave.jump.cost=0.35,
                        voiced.voiceless.cost=0.14,
                        corr.only=FALSE,
                        windowSize=40,
                        min.filter.freq=70.0,
                        max.filter.freq=5000.0,
                        filters=250,
                        max.freq.components=1250.0,
                        subharmonics=15,
                        compression=0.84,
                        points.per.octave=48,
                        windowShape="Gaussian1",
                        relativeWidth=1.0,
                        toFile=TRUE,
                        explicitExt="pf0",
                        outputDirectory=NULL,
                        verbose=FALSE,
                        praat_path=NULL){
  
  
  
  if(! have_praat(praat_path)){
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
                         file.path("inst","praat","praat_pitch.praat"),
                         file.path(system.file(package = "superassp",mustWork = TRUE),"praat","praat_pitch.praat")
  )
  
  praat_dsp_directory <- make_dsp_environment()

  pitch <- cs_wrap_praat_script(praat_location = get_praat(),
                                script_code_to_run = readLines(praat_script),
                                directory=praat_dsp_directory,
                                return="last-argument")
  
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
    
    pitchTabFile <- tempfile(fileext = ".csv",tmpdir = praat_dsp_directory)
    
    #Required for preventing errors in the handoff of file names containing spaces and () characters
    # to Praat
    soundFile <- tempfile(fileext = ".wav",tmpdir = praat_dsp_directory)
    R.utils::createLink(soundFile,origSoundFile)
    #Alternative route - much slower
    #file.copy(origSoundFile,soundFile)
    
    # sentence SoundFile /Users/frkkan96/Desktop/kaa_yw_pb.wav
    # real BeginTime 0.0
    # real EndTime 0.0
    # real Time_step 0.005
    # real Window_length_(s) 0.040
    # real Minimum_f0 75.0
    # real Maximum_f0 600
    # boolean Very_accurate 1
    # natural Number_of_cancidates 15
    # real Silence_threshold 0.03
    # real Voicing_threshold 0.45
    # real Octave_cost 0.01
    # real Octave_jump_cost 0.35
    # real Voiced/unvoiced_cost 0.14
    # boolean Only_correlation_methods 1
    # real Minimum_filter_frequency_(Hz) 70.0
    # real Maximum_filter_frequency_(Hz) 5000.0
    # natural Number_of_filters 250
    # real Maximum_frequency_components 1250.0
    # natural Maximum_number_of_subharmonics 15
    # real Compression_factor 0.84
    # natural Number_of_points_per_octave 48
    # word WindowType Gaussian1
    # real RelativeWidth 1.0
    # sentence TrackOut /Users/frkkan96/Desktop/kaa_yw_pb.f0Tab

    outPitchTabFile <- pitch(soundFile,
                             beginTime,
                             endTime,
                             windowShift/1000, #Praat takes seconds
                             windowSize/1000,
                             minF,
                             maxF,
                             ifelse(very.accurate,1,0),
                             max.f0.candidates,
                             silence.threshold,
                             voicing.threshold,
                             octave.cost,
                             octave.jump.cost,
                             voiced.voiceless.cost,
                             ifelse(corr.only,1,0),
                             min.filter.freq,
                             max.filter.freq,
                             filters,
                             max.freq.components,
                             subharmonics,
                             compression,
                             points.per.octave,
                             windowShape="Gaussian1",
                             relativeWidth,
                             pitchTabFile)
    
    inTable <- read.csv(file=outPitchTabFile
                        ,header=TRUE
                        ,na.strings =c("--undefined--","NA","?"),
                        sep = ",")
    
    
    
    # We need the sound file to extract some information
    origSound <- wrassp::read.AsspDataObj(soundFile)
    
    sampleRate <-  1/ windowShift * 1000
    
    startTime = windowShift
    #inTable$time <- 1:nrow(inTable) * windowShift / 1000

    outDataObj = list()
    attr(outDataObj, "trackFormats") <- c("INT16", "INT16")
    #Use the time separation between second and pitch measurement time stamps to compute a sample frequency.
   
    attr(outDataObj, "sampleRate") <- sampleRate
    
    attr(outDataObj, "origFreq") <-  as.numeric(attr(origSound, "sampleRate"))
    startTime <- 1/sampleRate
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(nrow(inTable))
    class(outDataObj) = "AsspDataObj"
    
    wrassp::AsspFileFormat(outDataObj) <- "SSFF"
    wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
    
    # Cross-correlation track
    ccTable <- inTable %>%
      dplyr::select(cc) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
 
    
    noCCValues <- nrow(ccTable)
    names(ccTable) <- NULL
    outDataObj = wrassp::addTrack(outDataObj, "cc", as.matrix(ccTable[,1]), "INT16")
  
    # Auto-correlation track
    acTable <- inTable %>%
      dplyr::select(ac) %>%
      replace(is.na(.), 0) %>%
      dplyr::mutate(
        dplyr::across(
          tidyselect::everything(),as.integer))
  
    noACValues <- nrow(acTable)
    names(acTable) <- NULL
    outDataObj = wrassp::addTrack(outDataObj, "ac", as.matrix(acTable[,1]), "INT16")
    
    if(! corr.only){
      attr(outDataObj, "trackFormats") <- c("INT16", "INT16","INT16", "INT16")
      # SPINET track
      spinetTable <- inTable %>%
        dplyr::select(spinet) %>%
        replace(is.na(.), 0) %>%
        dplyr::mutate(
          dplyr::across(
            tidyselect::everything(),as.integer))
      
      names(spinetTable) <- NULL
      outDataObj = wrassp::addTrack(outDataObj, "spinet", as.matrix(spinetTable[,1]), "INT16")
      
      # SHS track
      shsTable <- inTable %>%
        dplyr::select(shs) %>%
        replace(is.na(.), 0) %>%
        dplyr::mutate(
          dplyr::across(
            tidyselect::everything(),as.integer))
      
      names(shsTable) <- NULL
      outDataObj = wrassp::addTrack(outDataObj, "shs", as.matrix(shsTable[,1]), "INT16")
    }
    
    #return(outDataObj)
    ## Apply fix from Emu-SDMS manual
    ##https://raw.githubusercontent.com/IPS-LMU/The-EMU-SDMS-Manual/master/R/praatToFormants2AsspDataObj.R
    
    # add missing values at the start as Praat sometimes
    # has very late start values which causes issues
    # in the SSFF file format as this sets the startRecord
    # depending on the start time of the first sample
    if( startTime > (1/sampleRate) ){
      
      nr_of_missing_samples = as.integer(floor(startTime / (1/sampleRate)))
      
      missing_cc_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$cc))
      missing_ac_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$ac))
      
      # prepend values
      outDataObj$cc = rbind(missing_cc_vals, outDataObj$cc)
      outDataObj$ac = rbind(missing_ac_vals, outDataObj$ac)
      
      if(! corr.only){
        missing_spinet_vals = matrix(0,
                                 nrow = nr_of_missing_samples,
                                 ncol = ncol(outDataObj$spinet))
        missing_shs_vals = matrix(0,
                                 nrow = nr_of_missing_samples,
                                 ncol = ncol(outDataObj$shs))
        
        outDataObj$spinet = rbind(missing_spinet_vals, outDataObj$spinet)
        outDataObj$shs = rbind(missing_shs_vals, outDataObj$shs)
      }
      
      # fix start time
      attr(outDataObj, "startTime") = startTime - nr_of_missing_samples * (1/sampleRate)
    }
    
    assertthat::assert_that(wrassp::is.AsspDataObj(outDataObj),
                            msg = paste("The AsspDataObj created by the praat_pitch function is invalid.\nPlease check the table file '",tabfile,"' for errors.",sep=""))
    
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
  
  #Make sure that we clear the environment so that we may run the script again without complaints
  clear_dsp_environment(  praat_dsp_directory )
  logger::log_trace("Computed an f0 track for the input signal file '",listOfFiles,"'.")
  if(!toFile){
    return(outDataObj)
  }
  
}

attr(praat_pitch,"ext") <-  c("pf0") 
attr(praat_pitch,"tracks") <-  c("cc","ac","spinet","shs")
attr(praat_pitch,"outputType") <-  c("SSFF")



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




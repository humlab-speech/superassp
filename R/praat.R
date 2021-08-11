

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
#' If the user only want to estimate formant frequencies, computing them using the function [wrassp::forest] is _much_ quicker, and the user should therefore mainly consider using this function `praat_formant_burg` only if the use case specicifally demands the use of the burg algorithm for computing formants, or if the user wants to also study the formant intensity levels (L_n) which  [wrassp::forest] does not do.
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
#' @param preemphasis the frequency from which a preemphasis will be applied..
#' @param windowShape the analysis window function used when extracting part of a sound file for analysis. De faults to "Hanning".
#' @param relativeWidth the relative width of the windowing function used.
#' @param spectWindowShape The shape of the windowing function used for constructing the spectrogram. 
#' @param spectResolution The frequency resolution of the spectrogram from which formant intensities will be collected.
#' @param toFile write the output to a file? The file will be written in  `outputDirectory`, if defined, or in the same directory as the soundfile. 
#' @param explicitExt the file extension that should be used.
#' @param outputDirectory set an explicit directory for where the signal file will be written. If not defined, the file will be written to the same directory as the sound file.
#' @param verbose Not implemented. Only included here for compatibility.  
#' @param praat_path give an explicit path for Praat. If the praat 
#'
#' @return The number of processed files, or an SSFF track data object (if `toFile=FALSE`) containing three fields ("F", "B" and "L") containing formant frequencies, bandwidth and intensities.
#' 
#' @export
#'



praat_formant_burg <- function(listOfFiles,
                               beginTime=0,
                               endTime=0,
                               windowShift=5.0,
                               numFormants=5.0,
                               maxhzformant=5500.0,
                               windowSize=30,
                               preemphasis=50.0,
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

    # form Compute a formant track
    # sentence SoundFile /Users/frkkan96/Desktop/kaa_yw_pb.wav
    # real BeginTime 0.0
    # real EndTime 0.0
    # real Time_step 0.005
    # real Number_of_formants 5.0
    # real MaxHzFormant 5500.0
    # real WindowLength 0.025
    # real Pre_emphasis 50.0
    # word WindowShape Gaussian1
    # real RelativeWidth 1.0
    # word Spectrogram_window_shape Gaussian
    # real Spectrogram_resolution 40.0
    # sentence TrackOut /Users/frkkan96/Desktop/kaa_yw_pb.FormantTab
    # endform    
    
    outFormantTabFile <- formant_burg(soundFile,
                            beginTime,
                            endTime,
                            windowShift/1000, #Praat takes seconds
                            numFormants,
                            maxhzformant,
                            windowSize/1000, #Praat takes seconds
                            preemphasis,
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

  if(toFile){
    return(length(outListOfFiles))
  }else{
    return(outDataObj)
  }
    
}

attr(praat_formant_burg,"ext") <-  c("pfm") 
attr(praat_formant_burg,"tracks") <-  c("F", "B" , "L")
attr(praat_formant_burg,"outputType") <-  c("SSFF")


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
  
  if(toFile){
    return(length(outListOfFiles))
  }else{
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
#' @return An Assp Data Object with a field \code{intensity} containing the
#'   intensity values (in dB) obtained for each analysis window (which will be
#'   \code{windowShift} ms apart)
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
  
  praat_intensity <- tjm.praat::wrap_praat_script(praat_location = get_praat(),
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
    
    intensityTabFile <- tempfile(fileext = ".csv")
    
    #Required for preventing errors in the handoff of file names containing spaces and () characters
    # to Praat
    soundFile <- tempfile(fileext = ".wav")
    R.utils::createLink(soundFile,origSoundFile)
    #Alternative route - much slower
    
    #function(listOfFiles,beginTime=0,endTime=0,windowShift=5.0,minF=80,subtractMean=TRUE,window="Gaussian1",relativeWidth=1.0,toFile=TRUE,explicitExt="fms",outputDirectory=NULL,verbose=FALSE,praat_path=NULL){
    
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
    
    outDataObj = wrassp::addTrack(outDataObj, "intensity", as.matrix(intensityTable), "INT16")
    
    
    ## Apply fix from Emu-SDMS manual
    ##https://raw.githubusercontent.com/IPS-LMU/The-EMU-SDMS-Manual/master/R/praatToFormants2AsspDataObj.R
    
    # add missing values at the start as Praat sometimes
    # has very late start values which causes issues
    # in the SSFF file format as this sets the startRecord
    # depending on the start time of the first sample
    if( startTime > (1/sampleRate) ){
      
      nr_of_missing_samples = as.integer(floor(startTime / (1/sampleRate)))
      
      missing_intensity_vals = matrix(0,
                               nrow = nr_of_missing_samples,
                               ncol = ncol(outDataObj$intensity))
      
      
      
      # prepend values
      outDataObj$intensity = rbind(missing_intensity_vals, outDataObj$intensity)

      
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
  
  if(toFile){
    return(length(outListOfFiles))
  }else{
    return(outDataObj)
  }
  
}

attr(praat_intensity,"ext") <-  c("int") 
attr(praat_intensity,"tracks") <-  c("int")
attr(praat_intensity,"outputType") <-  c("SSFF")




#' Compute spectral moments using praat
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
  

  pmoments <- tjm.praat::wrap_praat_script(praat_location = get_praat(),
                                             script_code_to_run = readLines(praat_script)
                                             ,return="last-argument")

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
    
    outputfile <- tempfile(fileext = ".csv")
    
    #Required for preventing errors in the handoff of file names containing spaces and () characters
    # to Praat
    soundFile <- tempfile(fileext = ".wav")
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
    
    #### Fourtn spectral moment : Skewness     #### 
    
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
  
  if(toFile){
    return(length(outListOfFiles))
  }else{
    return(outDataObj)
  }
  
}
attr(praat_moments,"ext") <-  c("pmo") 
attr(praat_moments,"tracks") <-  c("cog","sd","skew","kurt")
attr(praat_moments,"outputType") <-  c("SSFF")


# FOR INTERACTIVE TESTING
#praat_moments("~/Desktop/aaa.wav",toFile=FALSE) -> out1
#praat_moments("~/Desktop/aaa.wav",maximumFrequency = 15000,toFile=FALSE) -> out2
#Use this to mark that a Praat script is being developed (and the version in the installed copy of the library can not be used)
# library(emuR)
# # 
# library(rstudioapi)
# PRAAT_DEVEL = TRUE
# sv <- data.frame("absolute_file_path"=c("/Users/frkkan96/Documents/src/superassp/tests/signalfiles/msajc003.wav","/Users/frkkan96/Documents/src/superassp/tests/signalfiles/msajc003.wav"),
#                  "start"= c(2.1839846102785367,2.1839846102785367),
#                  "end"=c(2.2545007519283704,2.2545007519283704))
# cs <- data.frame("absolute_file_path"=c("/Users/frkkan96/Documents/src/superassp/tests/signalfiles/msajc003.wav","/Users/frkkan96/Documents/src/superassp/tests/signalfiles/msajc003.wav","/Users/frkkan96/Documents/src/superassp/tests/signalfiles/msajc003.wav"),
#                  "start"= c(0.8949000559078937, 0.19097576470165953,1.0297467127470492),
#                  "end"=c(2.0392409511025624,2.2545007519283704, 2.4524407986647456))
# 
# praat_avqi(svDF=sv,csDF=cs)
# list.files(file.path(tempdir(check=TRUE),"avqtemp"))

# FOR INTERACTIVE TESTING
#Use this to mark that a Praat script is being developed (and the version in the installed copy of the library can not be used)
#library(emuR)
# 
#library(rstudioapi)
#PRAAT_DEVEL = TRUE
#
# path2demoData = file.path(tempdir(),"emuR_demoData")
# unlink(path2demoData, recursive = TRUE)
# 
# emuR::create_emuRdemoData()
# 
# ae <- emuR::load_emuDB(file.path(path2demoData,"ae_emuDB"))
# 
# attr(praat_sauce,"tracks")
# 
# add_trackDefinition(ae,
#                      name="cpp",
#                      fileExtension ="psa",columnName = "cpp",
#                      onTheFlyFunctionName = "praat_sauce",onTheFlyParams = list(pitchTracking=TRUE,formantMeasures=TRUE,spectralMeasures=TRUE))
# 
# add_ssffTrackDefinition(ae,name="FORMANTS",columnName = "fm",fileExtension = "psa")
# add_ssffTrackDefinition(ae,name="Bandwidths",columnName = "bw",fileExtension = "psa")
# add_ssffTrackDefinition(ae,name="Harmonics",columnName = "H",fileExtension = "psa")
# add_ssffTrackDefinition(ae,name="HNR",columnName = "hnr",fileExtension = "psa")
# add_ssffTrackDefinition(ae,name="H25K",columnName = "H25K",fileExtension = "psa")
# add_ssffTrackDefinition(ae,name="f0",columnName = "f0",fileExtension = "psa")
# 
# add_perspective(ae,"PSA")
# 
# set_levelCanvasesOrder(ae,"PSA",c("Phonetic","Tone"))
# set_signalCanvasesOrder(ae,"PSA",c("OSCI","SPEC","FORMANTS","Bandwidths","H25K"))
# 
# set_specOverlay(ae,"PSA","FORMANTS")
# 
# add_perspective(ae,"Harmonics")
# set_levelCanvasesOrder(ae,"Harmonics",c("Phonetic"))
# set_signalCanvasesOrder(ae,"Harmonics",c("OSCI","SPEC","Harmonics","HNR","H25K"))
# 


#rm(dbConfig)
#rstudioapi::navigateToFile(file.path(ae$basePath,"ae_DBconfig.json"),line=395)
#httpuv::stopAllServers()
#serve(ae,autoOpenURL = NULL)


# library('testthat')
# test_file('tests/testthat/test_aaa_initDemoDatabase.R')
# test_file('tests/testthat/test_praat.R')




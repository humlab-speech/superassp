

has_praat <- function(praat_path=NULL){
  return(ifelse(is.na(get_praat(praat_path)) ,FALSE,TRUE))
}

get_praat <- function(praat_path=NULL){
  sysname <- as.vector(Sys.info()["sysname"])
  
  if(is.null(praat_path)){
    praat_path <- switch(sysname,
                         Darwin = "/Applications/Praat.app/Contents/MacOS/Praat",
                         Sys.which("praat")
    )    
  }
  if(! file.exists(praat_path)){
    praat_path <- NA
    warning("The function get_praat could not find your praat binary. Please provide a full path directly.")
  }
  
  return(praat_path)
  
}

praat_formant_burg <- function(soundFile,ssff_file=NULL,time_step=0.0,formants=5.0,maxhzformant=5500.0,windowlength=0.025,preemph=50.0,praat_path=NULL){
  if(! has_praat(praat_path)){
    stop("Could not find praat. Please specify a full path.")
  }
  
  if(!file.exists(soundFile)) {stop("Unable to find the sound file ",soundFile)}
  
  formant_burg <- wrap_praat_script(praat_location = get_praat(),
                                      script_code_to_run = readLines(file.path("src","formant_burg.praat"))
                                      ,return="last-argument")
  outFile <- ifelse(is.null(ssff_file),tempfile(fileext = ".Table"),ssff_file)
  
  tabfile <- formant_burg(soundFile,time_step,formants,maxhzformant,windowlength,preemph,outFile)
  
  # We need the sound file to extract some information
  
  origSound <- wrassp::read.AsspDataObj(soundFile)
  
  inTable  <- readr::read_csv(tabfile,na=c("--undefined--",""))
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
    select(starts_with("F",ignore.case = FALSE)) %>%
    replace(is.na(.), 0) %>%
    mutate(across(everything(),as.integer)) 
  
  noFormantsValues <- nrow(fmTable)
  
  outDataObj = wrassp::addTrack(outDataObj, "fm", as.matrix(fmTable), "INT16")
  
  bwTable <- inTable %>%
    dplyr::select(starts_with("F",ignore.case = FALSE)) %>%
    replace(is.na(.), 0) %>%
    mutate(across(everything(),as.integer)) 
  
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
  
  assertthat::assert_that(is.AsspDataObj(outDataObj),
                          msg = paste("The AsspDataObj created by the praat_formant_burg function is invalid.\nPlease check the table file '",tabfile,"' for errors.",sep=""))
  
  
  if(! is.null(ssff_file)){
    attr(outDataObj,"filePath") <- as.character(ssff_file)
    wrassp::write.AsspDataObj(dobj=outDataObj,file=ssff_file)
  }
  
  return(outDataObj)
}


# > str(fms)
# List of 2
# $ fm: int [1:581, 1:4] 0 1239 0 0 0 0 0 0 0 0 ...
# $ bw: int [1:581, 1:4] 0 644 0 0 0 0 0 0 0 0 ...
# - attr(*, "trackFormats")= chr [1:2] "INT16" "INT16"
# - attr(*, "sampleRate")= num 200
# - attr(*, "filePath")= chr "/Users/frkkan96/Documents/src/superassp/tests/msajc003.fms"
# - attr(*, "origFreq")= num 20000
# - attr(*, "startTime")= num 0.0025
# - attr(*, "startRecord")= int 1
# - attr(*, "endRecord")= int 581
# - attr(*, "class")= chr "AsspDataObj"
# - attr(*, "fileInfo")= int [1:2] 20 2


# > str(outDataObj)
# List of 2
# $ fm: num [1:463, 1:5] 0 0 0 0 981 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:5] "F1(Hz)" "F2(Hz)" "F3(Hz)" "F4(Hz)" ...
# $ bw: num [1:463, 1:5] 0 0 0 0 981 ...
# ..- attr(*, "dimnames")=List of 2
# .. ..$ : NULL
# .. ..$ : chr [1:5] "F1(Hz)" "F2(Hz)" "F3(Hz)" "F4(Hz)" ...
# - attr(*, "trackFormats")= chr [1:2] "INT16" "INT16"
# - attr(*, "sampleRate")= num 160
# - attr(*, "origFreq")= num 50000
# - attr(*, "startTime")='data.frame':	1 obs. of  1 variable:
#   ..$ time(s): num 0.00285
# - attr(*, "endRecord")= int 459
# - attr(*, "class")= chr "AsspDataObj"
# - attr(*, "fileInfo")= int [1:2] 20 2

# Error in wrassp::write.AsspDataObj(dobj = outDataObj, file = ssff_file) : 
#   REAL() can only be applied to a 'numeric', not a 'integer'

# FOR INTERACTIVE TESTING
# library(emuR)
# unlink(file.path(tempdir(),"emuR_demoData"),recursive=TRUE)
# emuR::create_emuRdemoData()
# load_emuDB(file.path(tempdir(),"emuR_demoData","ae_emuDB")) -> ae_test


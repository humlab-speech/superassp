

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
  sampleRate <-  as.integer(1 / (inTable[2,"time(s)"] - inTable[1,"time(s)"]))
  attr(outDataObj, "sampleRate") <- as.integer(sampleRate)
    
  attr(outDataObj, "origFreq") <-  as.integer(attr(origSound, "sampleRate"))
  startTime <- inTable[1,"time(s)"]
  attr(outDataObj, "startTime") <- as.integer(startTime)
  attr(outDataObj, "endRecord") = nrow(inTable)
  class(outDataObj) = "AsspDataObj"
  
  wrassp::AsspFileFormat(outDataObj) <- "SSFF"
  wrassp::AsspDataFormat(outDataObj) <- as.integer(2) # == binary
  
  fmTable <- inTable %>%
    select(starts_with("F",ignore.case = FALSE)) %>%
    mutate(across(everything(),as.integer))
  
  outDataObj = wrassp::addTrack(outDataObj, "fm", as.matrix(fmTable), "INT16")
  
  bwTable <- inTable %>%
    dplyr::select(starts_with("F",ignore.case = FALSE)) %>%
    dplyr::mutate(across(everything(),as.integer))
  
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
  
  if(is.null(ssff_file)){
    return(outDataObj) #Return the ssff track object if a filename was not provided
  }else{
    attr(outDataObj,"filePath") <- as.character(ssff_file)
    wrassp::write.AsspDataObj(dobj=outDataObj,file=ssff_file)
  }
  
}

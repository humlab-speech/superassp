##' Prepare a file path for signal processing functions
##' 
##' Normalise a list of filenames so that they can be passed to a signal processing function
##' 
##' @param listOfFiles The list of file names to process
##' @return A normalised list of filenames
##' @author Raphael Winkelmann

prepareFiles <- function(listOfFiles) {
    
	listOfFiles = gsub("^file://","", listOfFiles)
	listOfFiles = normalizePath(path.expand(listOfFiles))
    
	return(listOfFiles)
}


makeOutputDirectory <- function(outputDirectory,logToFile, funName){
  if (!is.null(outputDirectory) && is.character(outputDirectory)) {
    outputDirectory = normalizePath(path.expand(outputDirectory),mustWork = FALSE)
    
    if( file.exists(outputDirectory) ){
      finfo  <- file.info(outputDirectory)
      if (! finfo$isdir)
        cli::cli_abort("The path {.path outputDirectory} exists but is not a directory.")
    }else{
      if (!dir.create(outputDirectory, recursive=TRUE))
        cli::cli_abort("Unable to create the output directory {.path outputDirectory}.")
    }
    
    if(logToFile){
      logger::log_appender()
      cli::cli_inform("Storing the processing log in {.path {outputDirectory}}.")
      logger::log_threshold(logger::TRACE,namespace=funName)
      logger::log_layout(logger::layout_glue_generator(
        format = "{level} [{format(time, \"%Y-%m-%d %H:%M:%S\")}] {msg}"
      ),namespace=funName)
      logger::log_appender(logger::appender_file(
        file=file.path(outputDirectory,paste(funName,"log",sep="."))
      ),namespace=funName)
      logger::log_trace("{currCall}",namespace=funName)
    }
    
  }
}


convertInputMediaFiles <- function(listOfFiles,nativeFiletypes,preferedFiletype,knownLossless,funName,convertOverwrites,keepConverted,verbose){
  
  if(verbose){
    convert_pb <- list(name="Converting media files",
                       format="Preparation: converting to {.field {cli::pb_extra$preferedFT}} format {cli::pb_bar} {cli::pb_current}/{cli::pb_total}",
                       show_after=1,
                       clear=FALSE,
                       extra=list(preferedFT=preferedFiletype)
    )
  }else{
    convert_pb <- NULL
  }
  
  ## Check and fix input file paths
  listOfFiles = gsub("^file://","", listOfFiles)
  listOfFiles = normalizePath(path.expand(listOfFiles))
  
  if (is.null(listOfFiles) || length(listOfFiles) == 0 || ! all(file.exists(listOfFiles)) ) {
    cli::cli_abort(c("!"="The {.arg listOfFiles} has to contain a vector of working full paths to speech recordings."))
  }
  #Verify that the av library can read the input file, which is assumed to mean that it can be converted
  readcheck <- purrr::map(listOfFiles,mediacheck)
  cant_read <- purrr::map_lgl(readcheck,purrr::is_null) && tools::file_ext(listOfFiles) %in% nativeFiletypes
  if(any(cant_read)){
    cli::cli_abort("The input file {.path {listOfFiles[cant_read]} cannot be converted or read by {.fun {funName}}")
  }
  
  
  notLossless <- listOfFiles[! tools::file_ext(listOfFiles) %in% knownLossless]
  
  if(length(notLossless) > 0){
    cli::cli_warn(c("w"="Found {.val {length(notLossless)}} recording{?s} stored in lossy compression formats",
                    "i"="If lossy compression was used when storing the signal, the result {.fun {funName}} may not be accurate",
                    "x"="Please use a lossless file format (file extensions {.or {.val { knownLossless}}}) for acoustic analysis of speech recordings if possible."))
    if(verbose)
      cli::cli_inform(c("i"="Lossy files: {.file {basename(notLossless)}}"))
  }
  
  # Convertion code for non-"native" files
  
  listOfFilesDF <- data.frame(audio=listOfFiles,
                              output = paste(tools::file_path_sans_ext(listOfFiles),preferedFiletype,sep=".")
  )
  
  
  toConvert <- subset(listOfFilesDF, ! (tools::file_ext(listOfFiles) %in% nativeFiletypes) )
  
  if(!convertOverwrites){ 
    toConvert <- subset( toConvert,!file.exists(output) && audio != output)
  }
  
  if(nrow(toConvert) > 0 ){
    cli::cli_inform(c("Found {.val {nrow(toConvert)}} recording{?s} that require conversion",
                      "x"="If a file format is not nativelly suppored by {.fun {funName}} it will have to be converted before application of the function",
                      "i"="Please use {.or {.val {nativeFiletypes}}} which are nativelly supported by {.fun {funName}} to eliminate this conversion."))
    if(verbose){
      
      cli::cli_inform(c("i"="Converted files: {.file {basename(toConvert$audio)}}"))
    }
    
    
    purrr::pwalk(.l=toConvert,.f=av::av_audio_convert,verbose = FALSE,channels=1,format=NULL,.progress = convert_pb)
    
    # Here we explicitly make a new listOfFiles that points to wav files
    # so that upcoming legacy code can be used
    listOfFiles <- listOfFilesDF$output
    
  }
  return(list(listOfFiles,toConvert))
}

cleanupConvertedInputMediaFiles <- function(toConvert, keepConverted,verbose){
  #Clear the wavs created in the conversion step
  # It is assumed that files that had not been created in the conversion should also not be cleared
  if(! keepConverted && nrow(toConvert) > 0){
    if(verbose) cli::cli_inform("Cleaning up temporary (converted) versions of media files")
    purrr::pwalk(.l=toConvert,.f= \(audio, output) unlink(output,recursive = FALSE, force = FALSE, expand = FALSE))
  }
}



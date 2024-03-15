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


makeOutputDirectory <- function(outputDirectory,logToFile=FALSE, funName){
  
  if (!is.null(outputDirectory) && is.character(outputDirectory)) {
    outputDirectory = normalizePath(path.expand(outputDirectory),mustWork = FALSE)

    if( file.exists(outputDirectory) ){
      finfo  <- file.info(outputDirectory)
      if (! finfo$isdir) cli::cli_abort("The path {.path outputDirectory} exists but is not a directory.")
    }else{
      if (!dir.create(outputDirectory, recursive=TRUE)) cli::cli_abort("Unable to create the output directory {.path outputDirectory}.")
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


convertInputMediaFiles <- function(listOfFiles,beginTime, endTime, windowShift=5.0,nativeFiletypes,preferedFiletype,knownLossless,funName,keepConverted,verbose){
  
  if(length(listOfFiles) < 1) return(list(c(),c()))
  #Fix begin and endtime argument for C code
  if(is.null(beginTime)) beginTime <- 0
  if(is.null(endTime)) endTime <- 0  
  if(any(endTime <= beginTime && endTime != 0 )) cli::cli_abort(c("The {.par endTime} needs to come later than {.par beginTime}",
                                                  "x"="Times for {.file {basename(listOfFiles[endTime <= beginTime && endTime != 0])}} are difficult to interpret as a time window",
                                                  "i"="Start and end times {.val {paste(paste(beginTime[endTime <= beginTime && endTime != 0], endTime[endTime <= startTime && endTime != 0],sep=\"->\"))}}")
                                                  )
  
  
  # A function that is safe to use for checking the possibility to convert media

  mediacheck <- purrr::possibly(av::av_media_info, otherwise = NULL,quiet = TRUE)
  
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
  #Check begin and end times
  assertthat::assert_that( length(beginTime) == 1 || (length(beginTime) > 1 && length(beginTime) != length(listOfFiles)))
  assertthat::assert_that( length(endTime) == 1 || (length(endTime) > 1 && length(endTime) != length(endTime)))
 
  # Make a summary of input files and their output files (if conversion is needed)

  
  
  listOfFilesDF <- data.frame(audio=listOfFiles) |>
    dplyr::mutate(audio_ext =tools::file_ext(audio) ) |>
    dplyr::mutate(output_ext =ifelse(audio_ext %in% nativeFiletypes, audio_ext, preferedFiletype)) |>
    dplyr::mutate(lossless = audio_ext %in% knownLossless) |>
    dplyr::mutate(convert_file = (! audio_ext %in% nativeFiletypes) ) |>
    dplyr::mutate(convert_timewindow = ( convert_file &  (endTime != beginTime) )) |>
    dplyr::mutate(dsp_input = ifelse(convert_timewindow , 
                                  tempfile(pattern=basename(audio),fileext = paste0(".",output_ext)), ## Process a time window to temp file
                                  paste(tools::file_path_sans_ext(audio),output_ext,sep=".") #Place the convertered whole file next to the original
    ))
  
  notLossless <- listOfFilesDF[! listOfFilesDF$lossless,"audio"]
  
  if(length(notLossless) > 0){
    
    if(all(listOfFilesDF[["beginTime"]] == 0) && all(listOfFilesDF[["endTime"]] == 0)){
      losslessMessage <- c("w"="Found {.val {length(notLossless)}} recording{?s} stored in not optimal formats (e.g. using lossy compression)",
        "i"="If lossy compression was used when storing the signal, the result {.fun {funName}} may not be accurate",
        "x"="Please use a nativelly supported lossless file format ({.or {.val {intersect(nativeFiletypes,knownLossless)}}} files) for acoustic analysis of speech recordings whenever possible.")
    }else{
      losslessMessage <- c("w"="Found {.val {length(notLossless)}} recording{?s} stored in lossy compression formats",
                           "i"="If lossy compression was used when storing the signal, there may be issues in deducing the time window you wanted, and the result {.fun {funName}} may not be accurate",
                           "x"="Please use a nativelly supported lossless file format ({.or {.val {intersect(nativeFiletypes,knownLossless)}}} files) for acoustic analysis of speech recordings whenever possible.")
    }
      
    cli::cli_warn(losslessMessage)
    if(verbose)
      cli::cli_inform(c("i"="Files that are in lossy formats : {.file {basename(notLossless)}}"))
  }
  
  toConvert <- listOfFilesDF |>
    dplyr::filter(convert_file | convert_timewindow ) |> 
    dplyr::rename(output= dsp_input) |>
    dplyr::mutate(start_time=max(beginTime-(windowShift/1000),0)) |> # Negative values are assumed to  be an error and mean 0
    dplyr::mutate(duration= purrr::map_dbl(audio, ~ purrr::pluck(av::av_media_info(.x),"duration"))) |>
    dplyr::mutate(endTime = ifelse(endTime == 0 || is.null(endTime),duration, endTime)) |>
    dplyr::mutate(total_time=(endTime-beginTime)+(windowShift/1000) ) |>
    dplyr::select(audio,output,start_time, total_time) 
  
  if(nrow(toConvert) > 0 ){
    
    # TODO: Activate this code again when we know how to check that a file can be read
    # #Verify that the av library can read the input file, which is assumed to mean that it can be converted
    # cant_read <- toConvert  |>
    #   dplyr::select(audio) |>
    #   purrr::map(mediacheck) |>
    #   purrr::map_lgl(purrr::is_null) 
    # 
    # #return(listOfFilesDF[cant_read])
    # 
    # if(any(cant_read)){
    #   listOfFilesDF[cant_read]  |>
    #     dplyr::mutate(dirname= dirname(audio)) |>
    #     dplyr::group_by(dirname) |>
    #     dplyr::group_walk(~ cli::cli_alert_danger("The input {cli::qty(basename(.x$audio))} file{?s} {.file {basename(.x$audio)}} in {.path {(.y$dirname)}} cannot be read by {.fun {funName}} or converted."))
    #   cli::cli_alert_warning("Excluding {listOfFilesDF[cant_read]} file{'s} from conversion due to read error.")
    #   #Exlude from conversion
    #   toConvert <- toConvert |>
    #     dplyr::filter(!cant_read)
    # }
    
    
    cli::cli_inform(c("Found {.val {nrow(toConvert)}} recording{?s} that require conversion",
                      "x"="If a file format is not nativelly suppored by {.fun {funName}} it will have to be converted before application of the function",
                      "i"="Please use {.or {.val {nativeFiletypes}}} which are nativelly supported by {.fun {funName}} to eliminate this conversion."))
    if(verbose){
      
      cli::cli_inform(c("i"="Converted files: {.file {basename(toConvert$audio)}}"))
    }
    
    assertthat::assert_that(all(names(toConvert) %in% formalArgs(av::av_audio_convert) ))
    purrr::pwalk(.l=toConvert,.f=av::av_audio_convert,verbose = FALSE,channels=1,format=NULL,.progress = convert_pb)
    
    assertthat::assert_that(all(file.exists(toConvert$output) ))

  }

  # Here we explicitly make a new listOfFiles that points to wav files
  # so that upcoming legacy code can be used
  toClear <-  toConvert$output

                     

  return(list(listOfFilesDF,toConvert, toClear))
}

writeSSFFOutputFile <- function(ssffobj,filename,ext, outputDirectory=NULL,verbose=FALSE){
  if(!is.null(outputDirectory) && !is.character(outputDirectory) ) cli::cli_abort("Invalid output directory")
  
  if(!is.null(outputDirectory) && is.character(outputDirectory) && !dir.exists(outputDirectory)){
    dir.create(outputDirectory,recursive = TRUE,verbose=TRUE)
    if(verbose) cli::cli_inform("Creating output directory {.path {outputDirectory}}")
  }
  #here we create the output file name
  if(is.null(outputDirectory)){
    outputDirectory <- dirname(filename)
  }
  outputfile <- file.path(outputDirectory,paste(basename(tools::file_path_sans_ext(filename)),ext,sep="."))

  if(verbose) cli::cli_inform("Writing SSFF object with tracks {.field {names(ssffobj)}} to output file {.file {outputfile}}")
  write.AsspDataObj(ssffobj,outputfile)
  return(file.exists(outputfile))
}


  
cleanupConvertedInputMediaFiles <- function(toClear, keepConverted,verbose){

  #Clear the wavs created in the conversion step
  # It is assumed that files that had not been created in the conversion should also not be cleared
  if(! keepConverted && length(toClear) > 0){
    if(verbose) cli::cli_inform("Cleaning up temporary (converted) versions of media files")
    unlink(toClear,recursive = FALSE, force = FALSE, expand = FALSE)
  }

}


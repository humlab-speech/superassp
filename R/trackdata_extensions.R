#' Get available tracks for a function or an SSFF file
#'
#'
#' @param x The name of function defined to output trackdata, or the full path to a file that can be read using [wrassp::read.AsspDataObj].
#'
#' @return A vector of tracks that the function is defined to return, or are contained within the file.
#' @export
#'
#' @examples
#' get_definedtracks("forest")
#' get_definedtracks("praat_formant_burg")
#'
get_definedtracks <- function(x){
  
  if(is.character(x) && file.exists(x)){
    tryCatch({
      read.AsspDataObj(x) -> r
      tr <- tracks.AsspDataObj(r)
    },error=function(e) {stop("Could not open the file ",x," for reading. Please ensure that it is an SSFF file!")})
    return(tr)
  }
  
  if(is.character(x)){
    fun <- get0(x)
  }else{
    fun <- x
  }
  
  #Check that the function has been prepared for use with this function by
  # giving it the the required additional attributes "ext" and "tracks"
  if(is.null(attr(fun,"tracks")) ){
    
    stop("The function ",x," is not defined correctly.\nPlease provide attributes \"tracks\" for the function.\n See ?attr for details, as well as attributes(praat_formant_burg) for an example." )
    
  }
  return(attr(fun,"tracks"))
}

#' Get the (default) extension for an SSFF producing function or a signal file 
#'
#'
#' @param x The name of function defined to output trackdata, or the full path to a file that can be read using [wrassp::read.AsspDataObj].
#'
#' @return A string indicating the default file extension of the SSFF generating function, or the file extension of the signal file.
#' @export
#'
#' @examples
#' get_extension("forest")
#' get_extension("praat_formant_burg")
#'

get_extension <- function(x){
  
  if(is.character(x) && file.exists(x)){
    tryCatch({
      read.AsspDataObj(x,begin=0, end=1, samples = TRUE) -> r
    },error=function(e) {stop("Could not open the file ",x," for reading. Please ensure that it is an SSFF file!")})
    return(tools::file_ext(x))
  }
  
  if(is.character(x)){
    fun <- get0(x)
  }else{
    fun <- x
  }
  
  #Check that the function has been prepared for use with this function by
  # giving it the the required additional attributes "ext"
  if(is.null(attr(fun,"ext")) ){
    
    stop("The function ",x," is not defined correctly.\nPlease provide attributes \"ext\" for the function.\n See ?attr for details, as well as attributes(praat_formant_burg) for an example." )
    
  }
  return(attr(fun,"ext"))
  
}


#' Get the return format of a wrassp/superassp speech signal processing function 
#'
#' Speech signal processing functions in superassp may return an SSFF track with one or more columns, or a simple, not nested list. This function may be used to learn the output type of a specific function. 
#' @param x The name of a speech signal processing function that are defined in the superassp or wrassp packages.
#' @param package The name of the package where the function is defined.
#'
#' @return Either "SSFF" or "list".
#' @export
#'
#' @examples
#' get_outputType("forest")
#' get_outputType("praat_avqi")
#' get_outputType("praat_formant_burg")

get_outputType <- function(x,package="superassp"){
  
  
  if(is.character(x)){
    fun <- utils::getFromNamespace(x,package)
  }else{
    fun <- x
  }
 
  #Check that the function has been prepared for use with this function by
  # giving it the the required additional attributes "ext" and "tracks"
  if(is.null(attr(fun,"outputType")) ){
    
    stop("The function ",x," is not defined correctly. See ?attr for details, as well as attributes(superassp::praat_formant_burg) for an example." )
    
  }
  return(attr(fun,"outputType"))
}
  
#' Summary table of superassp DSP function output
#' 
#' The summary table produced by this function lists the DSP function names,
#' default file extension, and a summary of produced SSFF tracks in the output file, or 
#' the number of fields in slice producing functions. 
#' 
#' @details 
#' If the number of tracks or fields produced by the function is very large, then the output 
#' is truncated to a summary of the number of tracks.  
#' 
#' @return 
#' A data.frame with function names as row labels, and with "extension" and "tracks" columns.
#' The output is ordered by file extension in alphabetical order by default to make it iseasier 
#' to make sure DSP data are not overwritten when multiple functions are applied to the same recordings.
#'
#' @examples
#' superassp_summary()
superassp_summary <- function(){
  funs <- stringr::str_sort(ls("package:superassp"))
  
  summaryTable <- data.frame(row.names = funs)
  
  for(f in funs){
    
    tr <- attr(get(f),"tracks")
    ext <- attr(get(f),"ext")
    type <- attr(get(f),"outputType")
    if(!is.null(ext)) summaryTable[f,"extension"] <-paste( ext,collapse = ";")
    if(!is.null(ext)) summaryTable[f,"extension"] <-paste( ext,collapse = ";")
    
    if(! is.null(type) ){
      ntracks <- length(tr)
      if (type == "SSFF"){
        trlabel <- " tracks)"
      }
      if (type == "list"){
        trlabel <- " fields)"
      }
      
      if(ntracks > 3){
        trs <- paste0(paste(c(tr[1:2],"..."),collapse = ","), " (",ntracks,trlabel)
        if(stringr::str_length(trs) > 10){
          trs <- paste0("(",ntracks,trlabel)
        }
      }else{
        
        trs <- paste(tr,collapse = ",")
        
      }
      
      summaryTable[f,"tracks"] <- paste(trs,collapse = ",")
    }else{
      
    }
  }
  summaryTable <- summaryTable[complete.cases(summaryTable),] %>% 
    dplyr::arrange(extension)
  
  return(summaryTable)
}



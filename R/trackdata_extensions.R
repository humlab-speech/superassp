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
get_definedtracks <- function(x,package="superassp"){
  
  if(is.character(x) && file.exists(x)){
    tryCatch({
      read.AsspDataObj(x) -> r
      tr <- tracks.AsspDataObj(r)
    },error=function(e) {stop("Could not open the file ",x," for reading. Please ensure that it is an SSFF file!")})
    return(tr)
  }
  
  if(is.character(x)){
    fun <- utils::getFromNamespace(x,package)
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

get_extension <- function(x,package="superassp"){
  
  if(is.character(x) && file.exists(x)){
    tryCatch({
      read.AsspDataObj(x,begin=0, end=1, samples = TRUE) -> r
    },error=function(e) {stop("Could not open the file ",x," for reading. Please ensure that it is an SSFF file!")})
    return(tools::file_ext(x))
  }
  
  if(is.character(x)){
    fun <- utils::getFromNamespace(x,package)
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
  


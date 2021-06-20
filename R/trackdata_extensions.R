#' Get available tracks for a function or an SSFF file
#'
#'
#' @param x The name of function defined to output trackdata, or the full path to a file that can be read into a [wrassp::AsspDataObj].
#'
#' @return A vector of tracks that the function is defined to return, or are contained within the file.
#' @export
#'
#' @examples
#' get_definedtracks("forest")
#' get_definedtracks("praat_formant_burg")
#'
get_definedtracks <- function(x){
  
  if(file.exists(x)){
    tryCatch({
      read.AsspDataObj(x) -> r
      tr <- tracks.AsspDataObj(r)
    },error=function(e) {stop("Could not open the file ",x," for reading. Please ensure that it is an SSFF file!")})
    return(tr)
  }
  
  if( is.null(x) || !is.null(wrassp::wrasspOutputInfos[[x]])){
    #Wrassp function
    return(wrassp::wrasspOutputInfos[[x]][["tracks"]])
  } else {
    if(exists(x) && is.function(get(x))){
      
      fun <- get(x)
      #Check that the function has been prepared for use with this function by
      # giving it the the required additional attributes "ext" and "tracks"
      if(is.null(attr(fun,"tracks")) ){
        
        warning("The function ",onTheFlyFunctionName," is not defined correctly. NULL is returned. \nPlease provide it with the attributes \"ext\" and \"tracks\".\n See ?attr for details, as well as attributes(praat_formant_burg) for an example." )
        
      }
      return(attr(fun,"tracks"))
    }
  }

  

  
}

#' Get the (default) extension for an SSFF producing function or a signal file 
#'
#'
#' @param x The name of function defined to output trackdata, or the full path to a file that can be read into a [wrassp::AsspDataObj].
#'
#' @return A string indicating the default file extension of the SSFF generating function, or the file extension of the signal file.
#' @export
#'
#' @examples
#' get_extension("forest")
#' get_extension("praat_formant_burg")
#'

get_extension <- function(x){
  
  if(file.exists(x)){
    tryCatch({
      read.AsspDataObj(x) -> r
    },error=function(e) {stop("Could not open the file ",x," for reading. Please ensure that it is an SSFF file!")})
    return(tools::file_ext(r))
  }
  
  if(exists(x) & is.function(get(x))){
    
    fun <- get(x)
    #Check that the function has been prepared for use with this function by
    # giving it the the required additional attributes "ext" and "tracks"
    if(!is.null(attr(fun,"ext")) ){
      return(attr(fun,"ext"))
    }
  }
  
  if( is.null(x) || !is.null(wrassp::wrasspOutputInfos[[x]])){
    #Wrassp function
    return(wrassp::wrasspOutputInfos[[x]][["ext"]])
  }
  
}
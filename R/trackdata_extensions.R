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
  
  if(exists(x) & is.function(get(x))){
    
    fun <- get(x)
    #Check that the function has been prepared for use with this function by
    # giving it the the required additional attributes "ext" and "tracks"
    if(!is.null(attr(fun,"tracks")) ){
      return(attr(fun,"tracks"))
    }
  }
  
  if( is.null(x) || !is.null(wrassp::wrasspOutputInfos[[x]])){
    #Wrassp function
    return(wrassp::wrasspOutputInfos[[x]][["tracks"]])
  }
  
}

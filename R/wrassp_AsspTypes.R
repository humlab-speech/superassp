##' returns all valid AsspWindowTypes according to the assp library
##'
##' wrapper function for AsspWindowTypes of wrassp
##' @title AsspWindowTypes
##' @return vector containing window types
##' @author Raphael Winkelmann
##' @useDynLib superassp, .registration = TRUE
##' @export
'AsspWindowTypes' <- function(){
	
	return(.Call("AsspWindowTypes_", PACKAGE = "superassp"))
  
}

##' returns all valid AsspLpTypes according to the assp library
##'
##' wrapper function for AsspLpTypes of wrassp
##' @title AsspLpTypes
##' @return vector containing lp types
##' @author Raphael Winkelmann
##' @useDynLib superassp, .registration = TRUE
##' @export
'AsspLpTypes' <- function(){
	
	return(.Call("AsspLpTypes_", PACKAGE = "superassp"))

}

##' returns all valid AsspSpectTypes according to the assp library
##'
##' wrapper function for AsspSpectTypes of wrassp
##' @title AsspSpectTypes
##' @return vector containing spectrogram types
##' @author Raphael Winkelmann
##' @useDynLib superassp, .registration = TRUE
##' @export
'AsspSpectTypes' <- function(){
	
	return(.Call("AsspSpectTypes_", PACKAGE = "superassp"))

}

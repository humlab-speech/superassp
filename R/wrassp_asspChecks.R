##' checks if given string is a valid AsspWindowType according to the assp library
##' 
##' @title isAsspWindowType
##' @param windowName name of window
##' @return (BOOL) true if windowName is valid; false otherwise
##' @author Raphael Winkelmann
##' @useDynLib superassp, .registration = TRUE
##' @export
isAsspWindowType <- function(windowName) {
	if (missing(windowName)) {
		stop("No windowName given!")
	}

	toupper(windowName) %in% AsspWindowTypes()

}

##' checks if given string is a valid AsspLpType according to the assp library
##' 
##' @title isAsspLpType
##' @param lpName name of lp type
##' @return (BOOL) true if lpName is valid; false otherwise
##' @author Raphael Winkelmann
##' @useDynLib superassp, .registration = TRUE
##' @export
"isAsspLpType" <- function(lpName = NULL) {
	if (is.null(lpName)) {
		stop("No lpName given!")
	}

	lpTypes = AsspLpTypes()

	isValidLp = FALSE

	for (type in lpTypes) {
		if (lpName == type) {
			isValidLp = TRUE
			break
		}
	}
	return(isValidLp)
}

##' checks if given string is a valid AsspSpectType according to the assp library
##'
##' @title isAsspSpectType
##' @param spectName name of lp type
##' @return (BOOL) true if spectName is valid; false otherwise
##' @author Raphael Winkelmann
##' @useDynLib superassp, .registration = TRUE
##' @export
"isAsspSpectType" <- function(spectName = NULL) {
  if (is.null(spectName)) {
    stop("No lpName given!")
  }
  
  spectTypes = AsspSpectTypes()
  
  isValidSpect = FALSE
  
  for (type in spectTypes) {
    if (spectName == type) {
      isValidSpect = TRUE
      break
    }
  }
  return(isValidSpect)
}

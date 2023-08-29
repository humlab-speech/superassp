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

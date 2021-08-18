

#' Derivation of SSFF track objects 
#' 
#' This function takes an SSFF object or file and computes the `order`th derivative of the tracks in it. The user may also specify a lag of differentiation. If `lag=1`, then ordinary differences between consecutive values are computed. If `lag=2`, then the difference between the 1st and 3rd value will be returned, and so on. The user may specify an order of differentiation too, and thereby cause the differentiation to be conducted in multiple iterations.
#' 
#' Differentiation always results in loss of data, and the user may specify how to align the differentiation output. Initial zero padding values will be inserted so that the vector length of the input and the output will always be the same. If `padLeft=TRUE` (the default) the initial zero values will be inserted into the tracks so that the differentiation result aligns in time with the occurrence of a value change. That is, if `lag=2` a value in the output vector indicates that at that place in the input vector a change has happened of the indicated size compared to the value two positions back in the vector. If `padLeft=FALSE` and `lag=2` then the value indicates the change that will have occurred when looking two positions forward in the vector. This is likely an unusual use case, and therefore not the default behavior.
#' 
#' Padding the signal with zeros is performed after all iterations of differentiation have been performed completely, and the padding zeros will therefore never be differentiated themselves.
#' 
#'
#' @param inSSFF The SSFF oject, or a full path to a file that contains an SSFF object and may be read as such by [wrassp::read.AsspDataObj].
#' @param order The number of iterations in which the vector will be differentiated. The first order differentiation gives the size of changes in consecutive values (with an indicated lag). The second order differentiation gives the rate of change, and so on.
#' @param onlyTracks Only differentiate certain tracks, and leave the others as is. Defaults to process all tracks.
#' @param padLeft Should initial zeros be inserted into the vector from the left?
#' @param toFile Should the resulting SSFF object be written to file, or returned?
#' @param explicitExt By default, a character "d" will be prepended to the file name suffix when writing the output to file. The user can also specify an explicit extension which will be used instead.
#' @param overwrite Should an existing file be overwritten when writing the output?
#'
#' @return
#'   The function will return an SSFF object if `toFile` is `TRUE`. Otherwise, nothing is returned.
#' @export
#'
difftrack <- function(inSSFF, order=1,onlyTracks=NULL,padLeft=TRUE,toFile=TRUE,explicitExt=NULL,overwrite=FALSE){
  
  if(! class(inSSFF) %in% c("character","AsspDataObj")){
    stop("The 'difftrack' function can only be applies to SSFF objects or files containing such objects.")
  }
  if(is.character(inSSFF) && ! base::file.exists(base::normalizePath(inSSFF))){
    stop("The SSFF file does not exists.")
  }
  
  
  if(toFile ){
    # Make sure we could write the output file
    fp <- normalizePath(attr(i1,"filePath"))
    d <- dirname(fp)
    #This works because a directory "" does not exist. 
    if(!dir.exists(d)) stop("The directory '",d,"' does not exist so an output file cannot be created there.")
    #Construct an output file name
    if(is.null(explicitExt)){
      fp <- paste(tools::file_path_as_absolute(fp),explicitExt,sep=".")
    }else{
      fp <- paste0(tools::file_path_as_absolute(fp),".d",tools::file_ext(fp))
    }
   
    if(!overwrite && file.exists(fp)) stop("The file '",fp,"' already exists. Please set overwrite to TRUE if you wanted to overwrite it.")
  }

  tracks <- names(inSSFF)
  if(!is.null(onlyTracks)){
    
    if(class(onlyTracks) == "character"){
      tracks <- tracks[which(tracks %in% c("H","A"))]
    }
    if(class(onlyTracks) == "numeric"){
      tracks <- tracks[onlyTracks]
    }
  }
  
  outSSFF <- inSSFF
  
  if(padLeft){
    # In this case, vector of diffs will be pushed to later in the output vectors so that 
    # the the change will happen at that point in the vector
    lPad <- rep(0,order )
    rPad <- c()
  }else{
    lPad <- c()
    rPad <- rep(0,order)
  }
  
  for(tr in tracks){
    for(c in 1:ncol(outSSFF[[tr]])){
      newVec <- c(lPad,base::diff(inSSFF[[tr]][,c],lag=1,differences = order),rPad)
      #return(newVec)
      outSSFF[[tr]][,c] <- newVec
    }
  }
  #Rename all tracks to reflect the derivation
  dlab <- paste(rep("d",order),collapse = "")
  newNames <- paste0(dlab, names(outSSFF) )
  names(outSSFF) <- newNames
  if(!toFile){
    return(outSSFF)
  }else{
    wrassp::write.AsspDataObj(outSSFF,file = fp)
  }
}

## INTERACTIVE TESTING
# testFile <- file.path("tests","signalfiles","msajc003.wav")
# wrassp::forest(testFile,toFile=FALSE) -> inSSFF
# difftrack(inSSFF,toFile=FALSE,padLeft = TRUE) -> outSSFF

# tail(cbind(as.data.frame(outSSFF$bw),as.data.frame(inSSFF$bw)))
# 
# lag <- 2; order=1
# print(c(rep(0,order),diff(c(1,2,44,2,1),lag=lag,differences = order),rep(0,lag-1)))
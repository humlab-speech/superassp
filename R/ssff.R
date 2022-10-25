as_tibble.AsspDataObj <- function(x,track=1){
  df <- x[[track]]
  colnames(df) <- paste0("T",seq(1,ncol(df)))
  
  times <- seq(from=attr(x,"startTime"),
               by=1/attr(x,"sampleRate"),
               length.out=nrow(df))
  out <- cbind(data.frame(times=times), df)
  out <- as_tibble(out) %>%
    na_if(tidyselect::matches("T[0-9]+"))
  return(out)
  
}


#' Derivation of SSFF track objects 
#' 
#' This function takes an SSFF object or file and computes the `order`th derivative of the tracks in it. The user may also specify a lag of differentiation. If `lag=1`, then ordinary differences between consecutive values are computed. If `lag=2`, then the difference between the 1st and 3rd value will be returned, and so on. The user may specify an order of differentiation too, and thereby cause the differentiation to be conducted in multiple iterations.
#' 
#' Differentiation always results in loss of data, and the user may specify how to align the differentiation output. Initial zero padding values will be inserted so that the vector length of the input and the output will always be the same. If `padLeft=TRUE` (the default) the initial zero values will be inserted into the tracks so that the differentiation result aligns in time with the occurrence of a value change. That is, if `lag=2` a value in the output vector indicates that at that place in the input vector a change has happened of the indicated size compared to the value two positions back in the vector. If `padLeft=FALSE` and `lag=2` then the value indicates the change that will have occurred when looking two positions forward in the vector. This is likely an unusual use case, and therefore not the default behavior.
#' 
#' Padding the signal with zeros is performed after all iterations of differentiation have been performed completely, and the padding zeros will therefore never be differentiated themselves.
#' 
#'
#' @param inSSFF The SSFF object, or a full path to a file that contains an SSFF object and may be read as such by [wrassp::read.AsspDataObj].
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
differentiate <- function(inSSFF, order=1,onlyTracks=NULL,padLeft=TRUE,toFile=TRUE,explicitExt=NULL,overwrite=FALSE){
  
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



#' Compute the harmonic structure from f0 measurements
#'
#' This function takes a pre-computed f0 track and derive `n` harmonic tracks
#' from it so that the vector of f0 values are now a matrix with `n` columns.
#' Each column then encode the `n`th harmonic values.
#'
#' @param f0 An f0 track, either as an SSFF object or as the name of an SSFF formatted file.
#' @param n The number of harmonics to compute.
#' @param explicitExt The output file extension.
#' @param toFile boolean;Should the SSFF track be returned or stored on disc?
#'
#' @return The SSFF track object, if required. `NULL` otherwise.
#' @export
#' 
#' 
#' 
harmonics <- function(f0, n=10, explicitExt="har",toFile=TRUE){
  
  if(! wrassp::is.AsspDataObj(f0) && file.exists(f0)){
    #We have a name of a file and need to read it in
    wrassp::read.AsspDataObj(f0) -> f0
  }
  
  
  
  # Now we are sure to have an SSFF object
  #Deduce an output file extension
  ext <- ifelse(is.null(explicitExt),
                paste0("m",tools::file_ext(attr(f0,"filePath"))),
                explicitExt)
  
  out <- f0
  for(i in 1:length(f0)){
    if(ncol(f0[[i]]) > 1){
      stop("The harmonic funcion does not work for multidimensional tracks.")
    }
    out[[i]] <- f0[[i]] %*% t(seq(1,n,1))  
    
  }
  outPath <- paste(tools::file_path_sans_ext(attr(f0,"filePath")),ext,sep=".")
  attr(out,"filePath") <- outPath
  
  if(toFile){
    wrassp::write.AsspDataObj(out,file = outPath)
  }else{
    return(out)
  }
  
}
attr(harmonics,"ext") <-  c("har") 
attr(harmonics,"outputType") <-  c("SSFF")

F_boundaries <- function(x, columnName = "fm",explicitExt="fbo",toFile=TRUE){
  
  if(! wrassp::is.AsspDataObj(x) && file.exists(x)){
    #We have a name of a file and need to read it in
    wrassp::read.AsspDataObj(x) -> x
  }
  
  formants <- x[[columnName]]
  F2 <- formants[,2]
  F2 <- ifelse(F2 > 10, F2,NA )
  F1 <- formants[,1]
  F1 <- ifelse(F1 > 10, F1,NA )
  F1med <- median(F1,na.rm=TRUE)
  F1min <- min(F1, na.rm=TRUE)
  F1max <- max(F1, na.rm=TRUE)
  F2med <- median(F2,na.rm=TRUE)
  F2min <- min(F2, na.rm=TRUE)
  F2max <- max(F2, na.rm=TRUE)
  
  ch <- geometry::convhulln(F2,F1)
    #articulated::VSD(F2,F1)
  
  outf2norm <- ch$p[ ch$hull[,1] ]
  outf1norm <- ch$p[ ch$hull[,2] ]
  
  outf1 <- (F1max-F1min) * outf1norm + F1med
  
  outf2 <- (F2max-F2min) * outf2norm + F2med
  
  
  return(list(hull=ch, x=outf2, y=outf1))
}

## INTERACTIVE TESTING
# testFile <- file.path("tests","signalfiles","msajc003.wav")
# wrassp::forest(testFile,toFile=FALSE) -> inSSFF
# difftrack(inSSFF,toFile=FALSE,padLeft = TRUE) -> outSSFF

# tail(cbind(as.data.frame(outSSFF$bw),as.data.frame(inSSFF$bw)))
# 
# lag <- 2; order=1
# print(c(rep(0,order),diff(c(1,2,44,2,1),lag=lag,differences = order),rep(0,lag-1)))
#"/Users/frkkan96/Desktop/a1.wav" -> fi
#"/Users/frkkan96/Desktop/aaa.f0" -> f0
#"/Users/frkkan96/Desktop/aaa.fms" -> fm
#read.AsspDataObj(fi) -> a
#read.AsspDataObj(f0) -> af0
#wrassp::read.AsspDataObj(fm) -> afm

#harmonics(wrassp::ksvF0(fi,toFile=FALSE),toFile=FALSE) -> mult

#F_boundaries(afm) -> out


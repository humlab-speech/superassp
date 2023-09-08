# #' Convert an AsspDataObj to a tibble
# #'
# #' This function converts an `AsspDataObj` to [tibble::as_tibble] form.
# #' @details
# #' This function redefines the [tibble::as_tibble] method for `AsspDataObj` so that 
# #' the output columns works well with `reindeer::quantify` and to replicate the output of [emuR::get_trackdata].
# #' 
# #' Contrary to [emuR::get_trackdata] this function assumes that track data that are exactly zero (0) are actually missing measurements, 
# #' acknowleging that this is how missing measurements are stored in the SSFF file format. The user should supply an argument `na.zeros=FALSE` if 
# #' that assumption is risky in the SSFF file that is processed.
# #' 
# #' 
# #' 
# #' @param x An object of class `AsspDataObj` (usually created by calling [read.AsspDataObj])
# #' @param field An optional argument indicating either the field name or field index number to extract. If not given (NULL), all fields will be extracted.
# #' @param start The start time of the portition of the SSFF track that was converted to a tibble. Defaults to zero (0) which means that the extracted portion is expected to start at the beginning of the signal.
# #' @param na.zeros Replace all zero (0) values in the track data columns with `NA` value? Defaults to `TRUE` so that subsequent applications of summary statistics functions do not risk confusing the zero values as actual measurements.   
# #'
# #' @return A [tibble::tibble] containing columns `timed_orig`, `times_rel`, `times_norm`, followed by one column for each track and track field (separated by '_')
# #' that is, if the user has chosen to convert the output of [forest] to a tibble, then the tibble will have columns `times_orig times_rel times_norm  fm_1  fm_2  fm_3  fm_4  bw_1  bw_2  bw_3  bw_4`.
# #' If the user only wanted the first field (or the "fm" field), then the output tibble will have columns `times_orig times_rel times_norm  fm_1  fm_2  fm_3  fm_4`.
# #' 
# #' 
# #' @export
# #'
# #' 
# as_tibble.AsspDataObj <- function(x,field=NULL,start=0,na.zeros=TRUE){
#   
#   if(!is.null(field)){
#     if(is.numeric(field) && field <= length(tracks.AsspDataObj(x))){
#       field <- tracks.AsspDataObj(x)[field]
#     }
#     # Now remove all other fields
#     for(todel in setdiff(tracks.AsspDataObj(x),field)){
#       x <- delTrack(x,todel)
#     }
#   }
#   
#   baseDF <- :as_tibble.AsspDataObj(x)
#   fixColName <- function(x) {stringr::str_replace(x,"(.+)([0-9]+)$","\\1_\\2")}
#   
#   out <- baseDF %>%
#     #dplyr::rename(time=frame_time) %>%
#     dplyr::mutate(times_orig=(frame_time /1000) + start ,
#                   times_rel=seq(from=0,to=(attr(x,"endRecord")-1)* 1000/attr(x,"sampleRate") ,by=1000/attr(x,"sampleRate")),
#                   times_norm=times_rel / (max(times_rel) - min(times_rel))
#     ) %>%
#     dplyr::select(-frame_time) %>%
#     dplyr::relocate(times_orig, times_rel,times_norm) %>%
#     dplyr::rename_with(.fn=fixColName)
#   
#   if(na.zeros){
#     out <- out %>%
#       dplyr::mutate(dplyr::across(!times_orig & !times_rel & !times_norm, ~ dplyr::na_if(.,0)))
#   }
#   
#   return(out)
#   
# }



#' Derivation of SSFF track objects 
#' 
#' This function takes an SSFF object or file and computes the `order`th derivative of the tracks in it. The user may also specify a lag of differentiation. If `lag=1`, then ordinary differences between consecutive values are computed. If `lag=2`, then the difference between the 1st and 3rd value will be returned, and so on. The user may specify an order of differentiation too, and thereby cause the differentiation to be conducted in multiple iterations.
#' 
#' Differentiation always results in loss of data, and the user may specify how to align the differentiation output. Initial zero padding values will be inserted so that the vector length of the input and the output will always be the same. If `padLeft=TRUE` (the default) the initial zero values will be inserted into the tracks so that the differentiation result aligns in time with the occurrence of a value change. That is, if `lag=2` a value in the output vector indicates that at that place in the input vector a change has happened of the indicated size compared to the value two positions back in the vector. If `padLeft=FALSE` and `lag=2` then the value indicates the change that will have occurred when looking two positions forward in the vector. This is likely an unusual use case, and therefore not the default behavior.
#' 
#' Padding the signal with zeros is performed after all iterations of differentiation have been performed completely, and the padding zeros will therefore never be differentiated themselves.
#' 
#'
#' @param inSSFF The SSFF object, or a full path to a file that contains an SSFF object and may be read as such by [read.AsspDataObj].
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
    write.AsspDataObj(outSSFF,file = fp)
  }
}



#' Compute the harmonic frequency structure from f0 measurements
#'
#' This function takes a pre-computed f0 track and derive `n` harmonic tracks
#' from it so that the vector of f0 values are now a matrix with `n` columns.
#' Each column then encode the `n`th harmonic values.
#' 
#' @details 
#' The stored harmonic frequencies are simply multiples of the fundamental frequency (f0) track, and not derived independently from the speech signal. Therefore, errors in the frequency tracking of the f0 signal will be carried over to these tracks.
#' The primary use case for the track is to have have estimates of the harmonic frequencies to visualize harmonic frequency (n*f~0* ) against harmonic amplitude ( *L~{1-n}* ) .
#' 
#'
#' @param track An f0 track, either as an SSFF object or as the name of an SSFF formatted file. It is recommended that the 
#' 
#' @param column The name of the column to use as an f0 track.
#' @param n The number of harmonics to compute.
#' @param explicitExt The output file extension.
#' @param toFile boolean;Should the SSFF track be returned or stored on disc?
#'
#' @return The SSFF track object, if required. `NULL` otherwise.
#' @export
#' 
#' 
#' 
harmonics <- function(track, column="f0",n=5, explicitExt="har",toFile=TRUE){
  
  if(! is.AsspDataObj(track) && file.exists(track)){
    #We have a name of a file and need to read it in
    read.AsspDataObj(track) -> track
  }
  
  outDataObj = list()

  harTable <- as.data.frame(track[[column]] %*% t(seq(1,n,1)) ) 
    
  outPath <- paste(tools::file_path_sans_ext(attr(track,"filePath")),explicitExt,sep=".")
  
  #Copy attributes over
  attr(outDataObj, "trackFormats") <- c("INT16")
  attr(outDataObj, "sampleRate") <- attr(track, "sampleRate")
  attr(outDataObj, "origFreq") <-  attr(track, "origFreq")
  attr(outDataObj, "startTime") <- attr(track, "startTime")
  attr(outDataObj, "startRecord") <- attr(track, "startRecord")
  attr(outDataObj, "endRecord") <- attr(track, "endRecord")
  class(outDataObj) = "AsspDataObj"
  
  AsspFileFormat(outDataObj) <- "SSFF"
  AsspDataFormat(outDataObj) <- as.integer(2) # == binary
  
 
  
  noHARValues <- nrow(harTable)
  names(harTable) <- NULL
  outDataObj = addTrack(outDataObj, "har", as.matrix(harTable), "INT16")
  
  
  if(toFile){
    write.AsspDataObj(outDataObj,file = outPath)
  }else{
    return(outDataObj)
  }
  
}

attr(harmonics,"ext") <-  c("har") 
attr(harmonics,"outputType") <-  c("SSFF")

F_boundaries <- function(x, columnName = "fm",explicitExt="fbo",toFile=TRUE){
  
  if(! is.AsspDataObj(x) && file.exists(x)){
    #We have a name of a file and need to read it in
    read.AsspDataObj(x) -> x
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
  formantsRev <- as.matrix(na.omit(data.frame("F1"=F1,"F2"=F2)))

  ch <- geometry::convhulln(formantsRev,
                            output.options="FA")
    #articulated::VSD(F2,F1)
  
  outf2norm <- ch$p[ ch$hull[,1] ]
  outf1norm <- ch$p[ ch$hull[,2] ]
  
  outf1 <- (F1max-F1min) * outf1norm + F1med
  outf1 <- outf1norm
  
  outf2 <- (F2max-F2min) * outf2norm + F2med
  outf2 <- outf2norm
  
  return(list(hull=ch, x=outf2, y=outf1))
}

checkRWLossless <- function(x,knownLossless = c("wav","flac","aiff","wv","tta","caf")){
  outputextensions <- setdiff(knownLossless,tools::file_ext(x)) # Do not overwrite the original input file
  bn <- tools::file_path_sans_ext(x)
  
  outnames <- normalizePath(paste(bn,outputextensions,sep="."),mustWork =FALSE)
  combinations <- expand.grid(audio=c(normalizePath(x,mustWork =FALSE),outnames),
                                      output=outnames,
                              stringsAsFactors=FALSE) |>  
    dplyr::filter(output != audio)
  
  init <- subset(combinations, grepl(".*wav",audio))
  rest <- subset(combinations, ! grepl(".*wav",audio))
  init |> 
    dplyr::bind_rows(rest) |> 
    purrr::pwalk(av::av_audio_convert,verbose=FALSE)

}



#' Cut out a portion of an SSFF object
#'
#' @param obj An in-memory object of class `AsspDataObj`.
#' @param where The time point (relative to the signal's total duration) where the cutout will be centered (0.0-1.0)
#' @param n_preceeding The number of data samples to include in addition to, and *before*, the sample closest to the `where` time point. 
#' @param n_following  The number of data samples to include in addition to, and *after*, the sample closest to the `where` time point.
#'
#' @return an `AsspDataObj` object that contains just the data samples at, and possibly surrounding, the `where` relative time reference in the input obj, but with time references relative to the timeline of the original object.
#' @export

cut.AsspDataObj <- function(obj,where,n_preceeding,n_following){
  if(where > 1 || where < 0 ) cli::cli_abort("A {.arg where} in the 0.0-1.0 range is required.")
  records <- attr(obj,"endRecord") - attr(obj,"startRecord")  + 1
  at <- round(records * where,0)
  startTime <- attr(obj, "startTime")
  sr <- attr(obj, "sampleRate")
  start <- at - n_preceeding
  end <- at + n_following
  cutout <- seq(start,end,1)
  for(n in names(obj)){
    obj[[n]] <- obj[[n]][cutout, ]
  }
  attr(obj,"startRecord") <- as.integer(start)
  attr(obj,"endRecord") <- as.integer(end)
  attr(obj,"startTime") <- startTime + ( (start - 1)  / sr ) 
  return(obj)
}

## INTERACTIVE TESTING
# testFile <- file.path("tests","signalfiles","msajc003.wav")
# forest(testFile,toFile=FALSE) -> inSSFF
# difftrack(inSSFF,toFile=FALSE,padLeft = TRUE) -> outSSFF

# tail(cbind(as.data.frame(outSSFF$bw),as.data.frame(inSSFF$bw)))
# 
# lag <- 2; order=1
# print(c(rep(0,order),diff(c(1,2,44,2,1),lag=lag,differences = order),rep(0,lag-1)))
#"/Users/frkkan96/Desktop/a1.wav" -> fi
#"/Users/frkkan96/Desktop/aaa.f0" -> f0
#"/Users/frkkan96/Desktop/a1.fms" -> fm
#forest(fi)
#read.AsspDataObj(fi) -> a
#read.AsspDataObj(f0) -> af0
#read.AsspDataObj(fm) -> afm

#harmonics(ksvF0(fi,toFile=FALSE),toFile=FALSE) -> mult

#F_boundaries(afm) -> out


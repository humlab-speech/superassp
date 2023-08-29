
make_sliceFileName <- function(mediaFileName,fileExtention,outputDirectory=NULL){
  if(is.null(outputDirectory)){
    sliceFileName <- gsub("[.][a-zA-Z0-9]+$",paste0(".",fileExtention),mediaFileName)
  }else{
    sliceFileBasename <- gsub("[.][a-zA-Z0-9]+$",paste0(".",fileExtention),basename(mediaFileName))
    sliceFileName <- file.path(outputDirectory,sliceFileBasename)
  }
  return(sliceFileName)
}

ensure_sparseSliceFile <- function(mediaFileName,measures,fileExtention="sli",outputDirectory=NULL){
  

  sliceFileName <- make_sliceFileName(mediaFileName=mediaFileName,fileExtention=fileExtention,outputDirectory=outputDirectory)
  measurenames <- names(measures)
  
  dbHandle <- DBI::dbConnect(sliceFileName,drv = RSQLite::SQLite()) 
  
  mf <- read.AsspDataObj(mediaFile)
  sr <- rate.AsspDataObj(mf)
  hash <- digest::sha1(mediaFile)

  
  if("slices" %in% DBI::dbListTables(dbHandle)){
    required_fields <- c("start_sample","end_sample","samplerate","sha",measurenames)
    if(! all(required_fields  %in% DBI::dbListFields(dbHandle,"slices"))){
      stop("The 'slices' table of the slice file '",sliceFileName,"' is not correctly formanted. Please make sure that it has the fields ",paste(required_fields,collapse=","))
    }
  } else {
    #Create the table definition
    cslices <- "CREATE TABLE slices ( \
      `start_sample` INTEGER NOT NULL, \
      `end_sample` INTEGER NOT NULL,  \
      `samplerate` INTEGER NOT NULL, \
      `sha` TEXT NOT NULL, \
      PRIMARY KEY (start_sample, end_sample, sha) \
      );"
    dbExecute(dbHandle,cslices)
    
    for(m in measurenames){
      cat(m)
      type <- switch(class(measures[[m]]),
                     character = "TEXT",
                     logical = "INTEGER",
                     numeric = "REAL") 
      dbExecute(dbHandle,paste("ALTER TABLE slices ADD COLUMN ",m,type," DEFAULT NULL;"))
    }
    
  }
  DBI::dbDisconnect(dbHandle)
  return(file.exists(sliceFileName))
}




#' Provides the ability to store a multidimensional feature set related to a
#' part of a signal.
#'
#' This function allow you associate a feature vector of high dimensionality to
#' the part of the signal from which is computed. The part can be a time slice
#' (with defined start and end times), or the entire signal. The collection of
#' time windows from which data has been obtained is viewed as sparse collection
#' of slices. A slice defined by its start and end times. The key features of
#' the idea of sparse collection of slices is that even though the set of
#' features extracted for the slice collection is consistent, a consistent
#' spacing of slices against the time line of the signal file is not required.
#' This design choice helps with the storage of output of computationally
#' intensive procedures that summarize a whole recording or just the part of it
#' where it makes sense to compute it.
#'
#' The primary example of a feature set belonging to a sparse collection of
#' slices is a large collection of acoustic measurements of voice obtained from
#' just the portion of the recording where the participant produces a maximally
#' prolonged vowel. Computing the output may be computationaly intensive, and
#' applying it across the entire signal just makes no sense for the application.
#' As an example of this situation, please refer to the \code{\link{vat}}
#' function which takes approximately 60x the duration of the sample to
#' complete, and is only valid for vowel productions. The participant may also
#' produce more than one prolonged vowel in a recording, and there is therefore
#' no unified analysis window to apply to the entire speech file (such as what
#' is the case when constructing a spectrogram from a set of equally spaced
#' spectral slices). The specification of `start_sample` and `start_sample` for
#' each slice additionaly allows for partially overlapping slices in flexible
#' way, and can be added to the sparse collection iteratively. The user may, for
#' instance, apply the analysis to the entire vowel, and then just to a 2s
#' portion starting 1s into the vowel, and store the results in the same
#' datafile for later use.
#'
#' An example of a whole file feature set would be a high-dimensional acoustic
#' description of an entire speaker or recording. This set is not well suited to
#' be aligned to a particularportion of the speech signal, and in this case
#' `start_time` and `start_time` is not set (`NULL`). Even though the the
#' feature set is alined to a maximally large portion of the speech signal, it
#' is still considered a slice in the sparse collection of slices.
#'
#' Sparse collection of slices are flexible, but the implementation of them
#' assures that two important features are upheld at all times:
#'
#' 1. The results stored in a sparse collection of slices must be identially
#' shaped. That is, they may only hold the features defined when you first
#' insert data into them. It is of course possible to have multiple sparse
#' collections (with unique file extensions) where the user may store outputs of
#' different analyses associated with a single speech signal.
#'
#' 2. Two slices may not be associated with exacly the same portion of the
#' speech signal. If the user supplies a new result to the same time slice
#' (including a slice including the entier file), the old values will be
#' overwritten. This is to ensure that the user will always be able to deternine
#' which set of output values is obtained from the just the start and end times
#' supplied to the \code{\link{get_slicedata}}, and also assures the most
#' efficient retrieval of the data. If the user needs to, for instance, apply
#' the same analysis multiple times to the same portion of the signal file, then
#' the user should instead use multiple sparse slice collection files.
#'
#'
#' @param mediaFileName The signal file from which the set of values has been
#'   obtained. Primarily, this will be the name of the speech signal file, but
#'   all signal files handled well by the \code{\link[wrassp]{read.AsspDataObj}}
#'   and can return a sample rate and number of samples will likely work.
#' @param values A list of values to be stored. The identity of each measurement
#'   value should be indicates as informative list names, as it will help
#'   retrieval of values. If the list values are not named, they will be
#'   identified by an integer value instead. The first time the user stores
#'   values from an analyse in a particular sparse slice collection file will
#'   determine what features may be stored and what names may be used for them.
#' @param start_sample The first sample of the signal file that was submitted
#'   for analysis, and for which the results should now be stored. If
#'   \code{NULL}, all samples from the first sample to the \code{end_sample}
#'   will be included.
#' @param end_sample The last sample of the signal file that was submitted for
#'   analysis, and for which the results should now be stored. If \code{NULL},
#'   all samples until the end of the media file will be used.
#' @param fileExtention This file extension will be used when making the sparse
#'   slice name from the \code{mediaFileName}. Defaults to ".sli".
#' @param outputDirectory The directory where the slice file should be stored.
#'   If not defiled (NULL), the sparse slice file will placed in the same folder
#'   as the media file.
#'
#' @export
#'
#' 

store_slice <- function(mediaFileName,values,measureNames,start_sample=NULL,end_sample=NULL,fileExtention="sli",outputDirectory=NULL){
  
  if(!ensure_sparseSliceFile(mediaFileName=mediaFileName,fileExtention=fileExtention,outputDirectory=outputDirectory)){
    sliceFileName <- make_sliceFileName(mediaFileName=mediaFileName,fileExtention=fileExtention,outputDirectory=outputDirectory)
    stop("Could not find the slice file '",sliceFileName,"'")
  }
  
  #The data we want to store must be in a row format, so that the measure is 
    
}  

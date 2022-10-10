

get_VoiceAnalysisToolkitCommand <- function(){
  candidates <- c("/Applications/voice_analysis_directory/application/run_voice_analysis_directory.sh","/matlab_install/voice_analysis_directorystandaloneApplication/run_voice_analysis_directory.sh")
  
  ex <- file.exists(candidates)
  return(head(candidates[ex],n=1))
}

get_MatlabRuntime <- function(){
  candidates <- c("/Applications/MATLAB/MATLAB_Runtime/R2022b","/usr/local/MATLAB/MATLAB_Runtime/v912/")
  
  ex <- file.exists(candidates)
  return(head(candidates[ex],n=1))
}



#' Applies the Voice Analysis Toolbox on all vowel samples in a directory
#'
#' The Voice Analysis Toolbox \insertCite{Tsanas:2011cb,Tsanas:2012un,tsanas2013automatic}{superassp} applies a wide range of Matlab™-implemented voice
#' analysis procedures on a single sustained vowel (usually \[a::\]) and computes
#' 339 acoustic quantities. This function calls a compiled application to
#' compute the outcome measures from each sustained vowel recording in a
#' directory, and collect the results into a [base::data.frame].
#'
#' @details
#' The user should be aware that applying this procedure to a directory of sound files may take 9-40 times the total duration of the sound files to perform, depending on the machine and how the application was compiled.
#' The user should therefore make sure to capture the tibble once returned (or retrieve it immidiately from the `.Last.value` variable once the command completes).
#' 
#' Under the hood, the Voice Analysis Toolbox utilizes several other Matlab™ toolboxes which are also compiled into the runtime binary that performs the procedure. 
#' The permissive open source licences of these published external toolboxes \insertCite{little2006nonlinear,little2007exploiting,brookes2011voicebox}{superassp} are gratefully acknowledged.
#' 
#'
#' @param directory The directory where sustained vowel samples are stored. The
#'   directory needs to be writable by the user, since the compiled code will
#'   also store the results of computations there.
#'
#' @return A [tibble::tibble] with one row for each sustained vowel sample,
#'   and with 340 columns. The first column (`listOfFiles`) contains the file
#'   names of recordings. The following columns (2 to 340) then contain the
#'   measurement values for the recording.
#' @export
#' 
#' @references 
#'   \insertAllCited{}
#' 

voice_analysis <- function(directory){
  
  directory <- normalizePath(directory)
  currwd <- getwd()
  
  if(!dir.exists(directory)) stop("The directory '",directory,"' does not extist!")
  
  setwd(directory)
  
  results <- base::system2(get_VoiceAnalysisToolkitCommand(),
                c(get_MatlabRuntime(),directory), stderr = TRUE, stdout=TRUE)
  
  #voice_analysis_directory_output.mat
  logger::log_trace(paste(results,sep="\n"))
  matout <- R.matlab::readMat(file.path(directory,"voice_analysis_directory_output.mat"))
  
  #Restore the working directory
  setwd(currwd)
  
  nfiles <- length(purrr::pluck(matout,"va.results"))
  
  measurenames <- as.vector(unlist(purrr::pluck(matout,"va.results",1,1,3,1)))
  
  d <- matrix(NA, nrow=nfiles,ncol=(1+length(measurenames)))
  d <- as.data.frame(d)
  names(d) <- c("listOfFiles",measurenames)
  
  for(file in 1:nfiles ){
    filename <- as.vector(unlist(purrr::pluck(matout,"va.results",file,1,1,1)))
    values <-  as.vector(unlist(purrr::pluck(matout,"va.results",file,1,2,1)))
    d[file,1] <- filename
    d[file,2:(length(measurenames)+1)] <- values
  }
  
  return(tibble::as_tibble(d))
}





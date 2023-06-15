

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
#' compute the outcome measures from each sustained vowel recording , and collect the results into a [base::list].
#'
#' @details
#' The user should be aware that applying this procedure to a directory of sound files may take 9-40 times the total duration of the sound files to perform, depending on the machine and how the application was compiled.
#' The user should therefore make sure to capture the tibble once returned (or retrieve it immidiately from the `.Last.value` variable once the command completes).
#' 
#' Under the hood, the Voice Analysis Toolbox utilizes several other Matlab™ toolboxes which are also compiled into the runtime binary that performs the procedure. 
#' The permissive open source licences of these published external toolboxes \insertCite{little2006nonlinear,little2007exploiting,brookes2011voicebox}{superassp} are gratefully acknowledged.
#' 
#'
#' @param listOfFiles A list of input sound files. Currently, the function takes only .wav files.
#' @param knownLossless A vector of file extensions known to be associated with losslessly encoded speech data. Please note that this parameter is currently not used by the underlying code since only .wav files are processed.
#'
#' @return A named list of 339 voice measurements measurement values for the recording.
#' @export
#' 
#' @references 
#'   \insertAllCited{}
#' 

voice_analysis_toolkit <- function(listOfFiles,
                           knownLossless=c("wav","flac","aiff","aif")){
  
  listOfFiles <- normalizePath(listOfFiles)
  
  currwd <- getwd()
  

  innerFunction <- function(file){
    
    directory <- superassp:::make_dsp_environment()
    
    file.copy(from = file,to=directory)
    setwd(directory)
    results <- base::system2(get_VoiceAnalysisToolkitCommand(),
                             c(get_MatlabRuntime(),directory), stderr = TRUE, stdout=TRUE)
    

    matout <- R.matlab::readMat(file.path(directory,"voice_analysis_directory_output.mat"))
    #Restore the working directory
    setwd(currwd)
    
    nfiles <- length(purrr::pluck(matout,"va.results"))
    
    measurenames <- as.vector(unlist(purrr::pluck(matout,"va.results",1,1,3,1)))
    values <- as.vector(unlist(purrr::pluck(matout,"va.results",1,1,2,1)))
    outList <- as.list(values)
    names(outList) <- measurenames
    superassp:::clear_dsp_environment(directory)
    
    return(outList)
    
  }

  purrr::map(.x=listOfFiles,.f=innerFunction)

}
attr(voice_analysis_toolkit,"ext") <-  "vat" 
attr(voice_analysis_toolkit,"outputType") <-  "list"
attr(voice_analysis_toolkit,"nativeFiletypes") <-  c("wav")
#,"flac","aiff","aif","aifc","au","ogg","opus","mp3","mp4","m4a")



## Interactive testing
#voice_analysis(c("~/Desktop/a1.wav","~/Desktop/a2.wav")) -> out




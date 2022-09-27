
#' A utility function to make a concurrency-safe DSP environment
#' 
#' The function will make a unique directory associated with a
#' particular speech signal file, which will then be unique ans safe for calling 
#' Praat DSP functions concurrently.
#' 
#' Teardown of the environment should be done using the [clear_dsp_environment] function,
#' which should then be called using the same signal file name. The file path is used "as is" when 
#' creating the concurrency safe environment, so the user should make sure to do any file path
#'normalization prior to using this function and [clear_dsp_environment] so that the unique 
#'identifier of the environment being constructed and tore down are ensured to be the same.
#'
#' @param signalFileFullPath The full path of the signal file that should be processed. 
#' @param salt an additional string to "salt" the identifier of the environment by make it unique. 
#' Usually the name of the function setting the environment up.
#'
#' @return The full path of the constructed directory (string).
#' @export
#' @seealso clear_dsp_environment
#'
make_dsp_environment <- function(signalFileFullPath, salt){
  hash <- digest::sha1(paste(signalFileFullPath, salt,sep="_"))
  tdir <- file.path(tempdir(check = TRUE),hash)
  if(dir.exists(tdir)){
    logger::log_warn("There is an existing Praat DSP environment found for the file '",signalFileFullPath,"'.")
    
  }
  logger::log_debug("Creating a Praat DSP environment found for the file '",signalFileFullPath,"'.")
  dir.create(tdir,recursive=FALSE,showWarnings = FALSE)
  return(tdir)
}

#' Tear down a DSP environment
#'
#' This function will remove the environment that has been created to process a
#' specific speech signal file `signalFileFullPath`. The same signal file path
#' should be given to this function and the [make_dsp_environment] function to
#' ensure that the correct environment is removed.
#'
#'
#' @param signalFileFullPath The full path of the signal file that the
#'   environment was originally set up for.
#' @param salt an additional string to "salt" the identifier of the environment by make it unique. 
#' Usually the name of the function setting the environment up.
#'
#' @return The function returns `FALSE` if the function was *unable* to remove
#'   the environment. If the function returns `TRUE`, the environment associated
#'   with the speech signal file was successfully removed, or did not exist in
#'   the first place.This means that the user can be reasonably sure in this
#'   case to not have clashing environments for Praat to run in if this function returns 
#'   `TRUE`.
#'   
clear_dsp_environment <- function(signalFileFullPath, salt){
  hash <- digest::sha1(paste(signalFileFullPath, salt,sep="_"))
  tdir <- file.path(tempdir(check = TRUE),hash)
  if(file.exists(tdir)){
    logger::log_debug("A Praat DSP environment exists for the file '",signalFileFullPath,"'. Deleting it.")
    success <- unlink(tdir, recursive = TRUE)
  }else{
    success <- 0
  }
  return(ifelse(success == 0,TRUE,FALSE))
}


#' A concurrency safe way to wrap a Praat script into an R function
#'
#' This function is an extension of the [tjm.praat::wrap_praat_script] function
#' that wraps a Praat script into an R function, but with an extra demand to
#' specify a directory where the Praat script should be stored before execution.
#' This means that rather than having multiple  Praat script executing from the
#' same temporary directory (which is the case when using
#' [tjm.praat::wrap_praat_script]), an environment is expected to have been
#' craeted using the [make_dsp_environment] already, and supplied to this
#' function via the `directory` argument. Any additional sound files or batches of sound files 
#' will then only be available to the particular instance of Praat and DSP Praat functions can be called concurrently.
#'
#' @param praat_location path to the Praat executable
#' @param script_code_to_run Praat script to run
#' @param directory the full path of a directory set up by [make_dsp_environment]
#' @param return value to return. "last-argument" returns the last argument to the Praat script. "info-window" returns the contents of the Praat Info Window.
#'
#' @return see `return` argument
#' @export

#' 
cs_wrap_praat_script <- function (praat_location,
                                  script_code_to_run,
                                  directory,
                                  return = c("last-argument",
                                             "info-window")
){
  return <- match.arg(return)
  script_file_to_run <- tempfile(fileext = ".praat", tmpdir=directory)
  writeLines(script_code_to_run, con = script_file_to_run)
  function(...) {
    if (return == "info-window") {
      results <- system2(praat_location,
                         c(
                           "--utf8",
                           "--run",
                           shQuote(script_file_to_run),
                           vapply(list(...),
                                  shQuote, "")
                         ),
                         stdout = TRUE)
      return(results)
    }
    else if (return == "last-argument") {
      results <- system2(praat_location, c(
        "--utf8",
        "--run",
        shQuote(script_file_to_run),
        vapply(list(...),
               shQuote, "")
      ))
      return(...elt(...length()))
    }
  }
}



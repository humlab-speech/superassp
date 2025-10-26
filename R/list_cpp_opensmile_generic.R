#' Generic OpenSMILE Feature Extraction using C++
#'
#' Internal function that extracts OpenSMILE features using the C++ implementation
#'
#' @param file Path to audio file
#' @param config_name Name of the config file (without extension)
#' @param config_dir Directory path relative to inst/opensmile/config
#' @param feature_set_name Name for verbose output
#' @param beginTime Start time in seconds
#' @param endTime End time in seconds  
#' @param verbose Print processing information
#' @return Named list with acoustic features
#' @keywords internal
opensmile_extract_generic <- function(file, config_name, config_dir, 
                                     feature_set_name = "features",
                                     beginTime = 0, endTime = 0,
                                     verbose = FALSE) {
  
  # Convert time parameters
  bt <- if (beginTime == 0) 0 else beginTime
  et <- if (endTime == 0) NULL else endTime
  
  # Get config file path
  config_file <- system.file("opensmile", "config", config_dir, config_name,
                            package = "superassp")
  
  if (config_file == "" || !file.exists(config_file)) {
    stop("OpenSMILE config file not found: ", config_dir, "/", config_name)
  }
  
  # Load audio with av package (universal format support)
  # Resample to 16kHz as expected by OpenSMILE configs
  audio_obj <- av_to_asspDataObj(
    file,
    start_time = bt,
    end_time = et,
    target_sample_rate = 16000
  )
  
  # Call C++ extraction function
  result <- opensmile_extract_cpp(audio_obj, config_file, feature_set_name, verbose)
  
  return(result)
}


#' Compute eGeMAPS Features (C++ Implementation)
#'
#' Extracts Extended Geneva Minimalistic Acoustic Parameter Set (eGeMAPS) features
#' using OpenSMILE C++ library.
#'
#' @param file Path to audio file
#' @param beginTime Start time in seconds (default: 0)
#' @param endTime End time in seconds (default: 0 = end of file)
#' @param verbose Print processing information (default: FALSE)
#' @return Named list with 88 eGeMAPS features
#' @keywords internal
lst_eGeMAPS_cpp <- function(file, beginTime = 0, endTime = 0, verbose = FALSE) {
  opensmile_extract_generic(
    file,
    config_name = "eGeMAPSv02_external.conf",
    config_dir = "egemaps/v02",
    feature_set_name = "eGeMAPS",
    beginTime = beginTime,
    endTime = endTime,
    verbose = verbose
  )
}


#' Compute emobase Features (C++ Implementation) 
#'
#' Extracts emobase features using OpenSMILE C++ library.
#'
#' @param file Path to audio file
#' @param beginTime Start time in seconds (default: 0)
#' @param endTime End time in seconds (default: 0 = end of file)
#' @param verbose Print processing information (default: FALSE)
#' @return Named list with 988 emobase features
#' @keywords internal
lst_emobase_cpp <- function(file, beginTime = 0, endTime = 0, verbose = FALSE) {
  opensmile_extract_generic(
    file,
    config_name = "emobase_external.conf",
    config_dir = "emobase",
    feature_set_name = "emobase",
    beginTime = beginTime,
    endTime = endTime,
    verbose = verbose
  )
}


#' Compute ComParE 2016 Features (C++ Implementation)
#'
#' Extracts ComParE 2016 features using OpenSMILE C++ library.
#'
#' @param file Path to audio file
#' @param beginTime Start time in seconds (default: 0)
#' @param endTime End time in seconds (default: 0 = end of file)
#' @param verbose Print processing information (default: FALSE)
#' @return Named list with 6373 ComParE features
#' @keywords internal
lst_ComParE_2016_cpp <- function(file, beginTime = 0, endTime = 0, verbose = FALSE) {
  opensmile_extract_generic(
    file,
    config_name = "ComParE_2016_external.conf",
    config_dir = "compare16",
    feature_set_name = "ComParE 2016",
    beginTime = beginTime,
    endTime = endTime,
    verbose = verbose
  )
}

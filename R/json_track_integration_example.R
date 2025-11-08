#' Example: Integrating JSTF with lst_* Functions
#' 
#' This file demonstrates how to add toFile support to lst_* functions
#' using the JSON Track Format (JSTF).
#'
#' @name json_track_integration
#' @keywords internal
NULL

#' Example lst_* function with JSTF support
#'
#' Template showing how to integrate JSTF writing into existing lst_* functions.
#'
#' @param listOfFiles Character vector of audio file paths
#' @param beginTime Start time in seconds (default: 0.0)
#' @param endTime End time in seconds (default: 0.0 = full duration)
#' @param toFile Logical, write results to JSTF file (default: FALSE)
#' @param explicitExt File extension (default: auto from function name)
#' @param outputDirectory Output directory (default: same as input)
#' @param verbose Logical, print progress (default: TRUE)
#' @param ... Additional function-specific parameters
#'
#' @return If toFile=FALSE: list of results. If toFile=TRUE: invisibly returns file paths
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' # In-memory mode (existing behavior)
#' results <- lst_example("audio.wav")
#' # Returns: list(measure1 = 123, measure2 = 456, ...)
#' 
#' # File output mode (new behavior)
#' lst_example("audio.wav", toFile = TRUE)
#' # Creates: audio.exm (JSON track file)
#' 
#' # Read back later
#' track <- read_track("audio.exm")
#' df <- as.data.frame(track)
#' }
lst_example_with_jstf <- function(listOfFiles,
                                  beginTime = 0.0,
                                  endTime = 0.0,
                                  toFile = FALSE,
                                  explicitExt = "exm",
                                  outputDirectory = NULL,
                                  verbose = TRUE,
                                  ...) {
  
  # Normalize time parameters
  beginTime <- if(is.null(beginTime)) 0.0 else beginTime
  endTime <- if(is.null(endTime)) 0.0 else endTime
  n_files <- length(listOfFiles)
  beginTime <- fast_recycle_times(beginTime, n_files)
  endTime <- fast_recycle_times(endTime, n_files)
  
  # Process each file
  results_list <- vector("list", n_files)
  
  for (i in seq_along(listOfFiles)) {
    file_path <- listOfFiles[i]
    
    if (verbose) {
      message("Processing: ", basename(file_path))
    }
    
    # Load audio
    audio_data <- av::read_audio_bin(
      file_path,
      start_time = if (beginTime[i] > 0) beginTime[i] else NULL,
      end_time = if (endTime[i] > 0) endTime[i] else NULL
    )
    
    sample_rate <- attr(audio_data, "sample_rate")
    duration <- length(audio_data) / sample_rate
    
    # ========================================
    # YOUR DSP PROCESSING HERE
    # ========================================
    results <- list(
      measure1 = rnorm(1, 100, 10),
      measure2 = rnorm(1, 50, 5),
      measure3 = rnorm(1, 200, 20)
    )
    # ========================================
    
    if (toFile) {
      # Create JsonTrackObj
      json_obj <- create_json_track_obj(
        results = results,
        function_name = "lst_example",
        file_path = file_path,
        sample_rate = sample_rate,
        audio_duration = duration,
        beginTime = beginTime[i],
        endTime = if (endTime[i] > 0) endTime[i] else duration,
        parameters = list(...)
      )
      
      # Construct output path
      base_name <- tools::file_path_sans_ext(basename(file_path))
      out_dir <- if (is.null(outputDirectory)) dirname(file_path) else outputDirectory
      output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))
      
      # Write to file
      write_json_track(json_obj, output_path, pretty = FALSE)
      
      results_list[[i]] <- output_path
      
    } else {
      # Return results directly
      results_list[[i]] <- results
    }
  }
  
  # Simplify output for single file
  if (n_files == 1) {
    results_list <- results_list[[1]]
  }
  
  if (toFile) {
    return(invisible(results_list))
  } else {
    return(results_list)
  }
}

# Set function attributes
attr(lst_example_with_jstf, "ext") <- "exm"
attr(lst_example_with_jstf, "outputType") <- "JSTF"
attr(lst_example_with_jstf, "format") <- "JSON"
attr(lst_example_with_jstf, "nativeFiletypes") <- c("wav")

#' Update existing lst_* function to support JSTF
#'
#' Helper function to add toFile support to existing lst_* functions.
#' This is a minimal patch that can be applied to existing functions.
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' # Minimal changes to add JSTF support:
#' 
#' lst_vat <- function(listOfFiles, ..., 
#'                     toFile = FALSE,              # ADD
#'                     explicitExt = "vat",        # ADD
#'                     outputDirectory = NULL) {   # ADD
#'   
#'   # ... existing processing code ...
#'   results <- voice_analysis(audio)
#'   
#'   # ADD this block at the end:
#'   if (toFile) {
#'     json_obj <- create_json_track_obj(
#'       results, "lst_vat", file_path,
#'       sample_rate, duration, beginTime, endTime
#'     )
#'     output_path <- file.path(
#'       outputDirectory %||% dirname(file_path),
#'       paste0(tools::file_path_sans_ext(basename(file_path)), ".", explicitExt)
#'     )
#'     write_json_track(json_obj, output_path)
#'     return(invisible(output_path))
#'   }
#'   
#'   return(results)  # existing return
#' }
#' 
#' # Set attributes
#' attr(lst_vat, "ext") <- "vat"
#' attr(lst_vat, "outputType") <- "JSTF"
#' }
add_jstf_support_to_lst <- function() {
  message("See json_track_integration_example.R for implementation guide")
}

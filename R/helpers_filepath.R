##' Prepare a file path for signal processing functions (Rcpp-optimized)
##' 
##' Normalise a list of filenames so that they can be passed to a signal processing function
##' 
##' @param listOfFiles The list of file names to process
##' @return A normalised list of filenames
##' @author Raphael Winkelmann
##' @keywords internal
prepareFiles <- function(listOfFiles) {
  # Use Rcpp to strip file:// protocol
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles))
  return(listOfFiles)
}


##' Create and validate output directory
##' 
##' @param outputDirectory Character path to output directory
##' @param logToFile Logical, whether to create log file
##' @param funName Character, name of calling function
##' @keywords internal
makeOutputDirectory <- function(outputDirectory, logToFile = FALSE, funName) {
  
  if(is.null(outputDirectory) || !is.character(outputDirectory)) {
    return(invisible(NULL))
  }
  
  outputDirectory <- normalizePath(path.expand(outputDirectory), mustWork = FALSE)
  
  if(file.exists(outputDirectory)) {
    finfo <- file.info(outputDirectory)
    if(!finfo$isdir) {
      cli::cli_abort("The path {.path {outputDirectory}} exists but is not a directory.")
    }
  } else {
    if(!dir.create(outputDirectory, recursive = TRUE)) {
      cli::cli_abort("Unable to create the output directory {.path {outputDirectory}}.")
    }
  }
  
  if(logToFile) {
    cli::cli_inform("Storing the processing log in {.path {outputDirectory}}.")
    logger::log_threshold(logger::TRACE, namespace = funName)
    logger::log_layout(
      logger::layout_glue_generator(
        format = "{level} [{format(time, \"%Y-%m-%d %H:%M:%S\")}] {msg}"
      ),
      namespace = funName
    )
    logger::log_appender(
      logger::appender_file(
        file = file.path(outputDirectory, paste0(funName, ".log"))
      ),
      namespace = funName
    )
  }
  
  invisible(outputDirectory)
}


##' Convert input media files (Rcpp-optimized)
##' 
##' @param listOfFiles Character vector of input file paths
##' @param beginTime Numeric vector of begin times
##' @param endTime Numeric vector of end times
##' @param windowShift Numeric window shift in milliseconds
##' @param nativeFiletypes Character vector of natively supported formats
##' @param preferedFiletype Character string of preferred conversion format
##' @param knownLossless Character vector of known lossless formats
##' @param funName Character name of calling function
##' @param keepConverted Logical whether to keep converted files
##' @param verbose Logical whether to show progress messages
##' @return List with processed file information and files to clean
##' @keywords internal
convertInputMediaFiles <- function(listOfFiles, beginTime, endTime, 
                                   windowShift = 5.0, nativeFiletypes, 
                                   preferedFiletype, knownLossless, 
                                   funName, keepConverted, verbose) {
  
  if(length(listOfFiles) < 1) {
    return(list(
      data.frame(
        audio = character(), dsp_input = character(), 
        beginTime = numeric(), endTime = numeric(),
        stringsAsFactors = FALSE
      ),
      data.frame(),
      character()
    ))
  }
  
  # Fix begin and endtime argument for C code
  beginTime <- if(is.null(beginTime)) 0.0 else beginTime
  endTime <- if(is.null(endTime)) 0.0 else endTime
  
  # Use Rcpp for fast file preparation
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles))

  # Check file existence
  files_exist <- file.exists(listOfFiles)
  if(!all(files_exist)) {
    missing_files <- listOfFiles[!files_exist]
    cli::cli_abort(c(
      "!" = "Some files do not exist:",
      "x" = "{.file {fast_basename(missing_files)}}"
    ))
  }

  # Validate time windows using Rcpp
  fast_validate_times(beginTime, endTime, fast_basename(listOfFiles))

  # Ensure time vectors match file count
  n_files <- length(listOfFiles)
  if(length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if(length(endTime) == 1) endTime <- rep(endTime, n_files)

  # Use Rcpp for combined format checking (single pass, faster)
  format_check <- fast_check_file_formats(listOfFiles, nativeFiletypes, knownLossless)
  audio_ext <- format_check$extensions
  is_native <- format_check$is_native
  is_lossless <- format_check$is_lossless
  
  # Determine which files need conversion
  needs_conversion <- !is_native
  needs_timewindow <- needs_conversion & (endTime != beginTime)
  
  # Check for lossless formats
  not_lossless <- listOfFiles[!is_lossless]
  
  if(length(not_lossless) > 0) {
    lossless_msg <- if(all(beginTime == 0) && all(endTime == 0)) {
      c(
        "w" = "Found {length(not_lossless)} recording{?s} in non-optimal format{?s}",
        "i" = "Lossy compression may affect {.fun {funName}} accuracy",
        "x" = "Use native lossless formats: {.val {intersect(nativeFiletypes, knownLossless)}}"
      )
    } else {
      c(
        "w" = "Found {length(not_lossless)} recording{?s} with lossy compression",
        "i" = "May affect time window extraction and {.fun {funName}} accuracy",
        "x" = "Use native lossless formats: {.val {intersect(nativeFiletypes, knownLossless)}}"
      )
    }
    
    cli::cli_warn(lossless_msg)
    
    if(verbose) {
      cli::cli_inform(c(
        "i" = "Lossy format files: {.file {fast_basename(not_lossless)}}"
      ))
    }
  }
  
  # Use Rcpp to generate output paths efficiently
  output_ext <- ifelse(is_native, audio_ext, preferedFiletype)
  dsp_input <- fast_generate_output_paths(
    listOfFiles,
    preferedFiletype,
    is_native,
    needs_timewindow
  )
  
  # Build main dataframe
  listOfFilesDF <- data.frame(
    audio = listOfFiles,
    audio_ext = audio_ext,
    output_ext = output_ext,
    lossless = is_lossless,
    convert_file = needs_conversion,
    convert_timewindow = needs_timewindow,
    dsp_input = dsp_input,
    beginTime = ifelse(needs_timewindow, 0.0, beginTime),
    endTime = ifelse(needs_timewindow, 0.0, endTime),
    stringsAsFactors = FALSE
  )
  
  # Handle conversions if needed
  toConvert <- data.frame()
  toClear <- character()
  
  n_conversions <- sum(needs_conversion | needs_timewindow)
  
  if(n_conversions > 0) {
    cli::cli_inform(c(
      "Found {n_conversions} recording{?s} requiring conversion",
      "x" = "Non-native formats need conversion before {.fun {funName}}",
      "i" = "Use native formats to avoid this: {.val {nativeFiletypes}}"
    ))
    
    # Get indices needing conversion (0-based for Rcpp)
    conv_indices <- which(needs_conversion | needs_timewindow) - 1L
    
    if(verbose) {
      cli::cli_inform(c(
        "i" = "Converting: {.file {fast_basename(listOfFiles[conv_indices + 1L])}}"
      ))
    }
    
    # Get durations only when needed for time window calculations
    needs_duration <- needs_timewindow[conv_indices + 1L]

    durations <- numeric(length(conv_indices))
    if(any(needs_duration)) {
      # Only call av_media_info for files that need it
      files_need_duration <- which(needs_duration)
      for(i in files_need_duration) {
        idx <- conv_indices[i] + 1L
        info <- av::av_media_info(listOfFiles[idx])
        durations[i] <- if(is.null(info$duration)) 0.0 else info$duration
      }
    }

    # Calculate conversion parameters using Rcpp
    time_params <- fast_calculate_conversion_times(
      beginTime,
      endTime,
      durations,
      windowShift,
      conv_indices
    )
    
    # Build conversion dataframe using Rcpp
    toConvert <- fast_build_conversion_df(
      listOfFiles,
      dsp_input,
      time_params$start_time,
      time_params$total_time,
      conv_indices
    )
    
    # Setup progress reporting
    convert_pb <- if(verbose) {
      list(
        name = "Converting media files",
        format = "Converting to {.field {preferedFiletype}} {cli::pb_bar} {cli::pb_current}/{cli::pb_total}",
        show_after = 1,
        clear = FALSE
      )
    } else {
      NULL
    }
    
    # Perform conversions
    purrr::pwalk(
      toConvert,
      av::av_audio_convert,
      verbose = FALSE,
      channels = 1,
      format = NULL,
      .progress = convert_pb
    )
    
    # Verify conversions
    conversion_success <- file.exists(toConvert$output)
    
    if(!all(conversion_success)) {
      failed_indices <- which(!conversion_success)
      cli::cli_abort(c(
        "Some file conversions failed:",
        "x" = "{.file {fast_basename(toConvert$audio[failed_indices])}}"
      ))
    }
    
    toClear <- toConvert$output
  }
  
  return(list(listOfFilesDF, toConvert, toClear))
}


##' Write SSFF output file
##' 
##' @param ssffobj SSFF data object to write
##' @param filename Character string of original filename
##' @param ext Character string of output extension
##' @param outputDirectory Character string of output directory (optional)
##' @param verbose Logical whether to show progress messages
##' @return Logical indicating success
##' @keywords internal
writeSSFFOutputFile <- function(ssffobj, filename, ext, 
                                outputDirectory = NULL, verbose = FALSE) {
  
  if(!is.null(outputDirectory) && !is.character(outputDirectory)) {
    cli::cli_abort("Invalid output directory")
  }
  
  # Create output directory if needed
  if(!is.null(outputDirectory) && is.character(outputDirectory)) {
    if(!dir.exists(outputDirectory)) {
      dir.create(outputDirectory, recursive = TRUE, showWarnings = verbose)
      if(verbose) {
        cli::cli_inform("Creating output directory {.path {outputDirectory}}")
      }
    }
  }
  
  # Determine output directory
  out_dir <- if(is.null(outputDirectory)) {
    dirname(filename)
  } else {
    outputDirectory
  }
  
  # Use Rcpp for fast path construction
  base_name <- fast_file_path_sans_ext(c(filename))[1]
  base_name <- fast_basename(c(base_name))[1]
  outputfile <- file.path(out_dir, paste0(base_name, ".", ext))
  
  if(verbose) {
    cli::cli_inform(
      "Writing SSFF with tracks {.field {names(ssffobj)}} to {.file {basename(outputfile)}}"
    )
  }
  
  # Write the SSFF file
  write.AsspDataObj(ssffobj, outputfile)
  
  return(file.exists(outputfile))
}


##' Cleanup converted input media files
##' 
##' @param toClear Character vector of files to remove
##' @param keepConverted Logical whether to keep converted files
##' @param verbose Logical whether to show progress messages
##' @keywords internal
cleanupConvertedInputMediaFiles <- function(toClear, keepConverted, verbose) {
  
  # Early exit if nothing to clear or keeping files
  if(keepConverted || length(toClear) == 0) {
    return(invisible(NULL))
  }
  
  if(verbose) {
    cli::cli_inform("Cleaning up {length(toClear)} temporary file{?s}")
  }
  
  # Remove temporary files
  unlink(toClear, recursive = FALSE, force = FALSE, expand = FALSE)
  
  # Verify cleanup
  still_exist <- file.exists(toClear)
  if(any(still_exist) && verbose) {
    cli::cli_warn(c(
      "!" = "Failed to remove {sum(still_exist)} temporary file{?s}",
      "i" = "{.file {fast_basename(toClear[still_exist])}}"
    ))
  }
  
  invisible(NULL)
}


#' Detect creaky voice using the Union Method
#'
#' The Union Method \insertCite{White2022}{superassp} combines two complementary
#' approaches for detecting creaky voice: the AM (antimode) method based on F0
#' distribution analysis, and the CD (creak detector) method using a trained
#' artificial neural network with 36 acoustic features.
#'
#' This function uses the \code{av} package to load audio files in-memory,
#' supporting all formats (WAV, MP3, MP4, etc.). It returns comprehensive
#' track data including both individual methods (AM and CD) and their union,
#' along with intermediate features used in classification.
#'
#' @param listOfFiles Vector of file paths to audio files
#' @param beginTime Start time in seconds (not yet implemented)
#' @param endTime End time in seconds (not yet implemented)
#' @param windowShift Frame shift in milliseconds (default 10ms)
#' @param use.am Use AM (antimode) method for detection
#' @param use.cd Use CD (trained ANN) method for detection
#' @param cd.threshold Classification threshold for CD method (0-1, default 0.3)
#' @param use.reaper Use REAPER for F0 estimation \insertCite{Ghahremani2014}{superassp} (requires pyreaper in Python)
#' @param return.features Include intermediate feature tracks in output
#' @param return.probabilities Include classification probabilities (not just binary decisions)
#' @param explicitExt File extension for output files (default "crk")
#' @param outputDirectory Directory for output files (default NULL = same as input)
#' @param toFile Whether to write results to disk (default TRUE)
#' @param conda.env Conda environment name for Python (default NULL = auto-detect)
#'
#' @return An AsspDataObj (or list of AsspDataObjs) containing multiple tracks:
#' \describe{
#'   \item{AM_creak}{Binary decisions from AM method (if use.am=TRUE)}
#'   \item{CD_creak}{Binary decisions from CD method (if use.cd=TRUE)}
#'   \item{union_creak}{Union of AM and CD decisions}
#'   \item{CD_prob}{Classification probability from ANN (if return.probabilities=TRUE)}
#'   \item{F0}{Fundamental frequency track (if use.am=TRUE)}
#'   \item{H2_H1}{Harmonic amplitude difference (if return.features=TRUE)}
#'   \item{peak_prominence}{Residual peak prominence (if return.features=TRUE)}
#'   \item{IFP}{Intraframe periodicity \insertCite{Ishi2008}{superassp} (if return.features=TRUE)}
#'   \item{IPS}{Interpulse similarity \insertCite{Ishi2008}{superassp} (if return.features=TRUE)}
#'   \item{PwP_rise}{Power pulse waveform rising \insertCite{Ishi2008}{superassp} (if return.features=TRUE)}
#'   \item{PwP_fall}{Power pulse waveform falling \insertCite{Ishi2008}{superassp} (if return.features=TRUE)}
#' }
#' Attributes include antimode value and method details.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage - Union method with both AM and CD
#' creak <- trk_creak_union("speech.wav")
#' 
#' # Access individual methods
#' plot(creak$time, creak$AM_creak, type='l', col='blue')
#' lines(creak$time, creak$CD_creak, col='red')
#' lines(creak$time, creak$union_creak, col='black', lwd=2)
#' 
#' # Include probabilities and features
#' creak_full <- trk_creak_union("speech.wav",
#'                                return.probabilities = TRUE,
#'                                return.features = TRUE)
#' 
#' # Plot CD probability
#' plot(creak_full$time, creak_full$CD_prob, type='l',
#'      xlab="Time (s)", ylab="Creak Probability")
#' 
#' # Access Ishi features
#' plot(creak_full$time, creak_full$IFP, type='l',
#'      xlab="Time (s)", ylab="Intraframe Periodicity")
#' 
#' # Batch processing
#' files <- c("speaker1.wav", "speaker2.wav", "speaker3.wav")
#' results <- trk_creak_union(files, toFile=TRUE, outputDirectory="output")
#' }
#'
#' @references
#'   \insertAllCited{}
#'
trk_creak_union <- function(listOfFiles,
                             beginTime = 0,
                             endTime = 0,
                             windowShift = 10,
                             use.am = TRUE,
                             use.cd = TRUE,
                             cd.threshold = 0.3,
                             use.reaper = TRUE,
                             return.features = FALSE,
                             return.probabilities = FALSE,
                             explicitExt = "crk",
                             outputDirectory = NULL,
                             toFile = TRUE,
                             conda.env = NULL) {
  
  # Validate inputs
  if (!requireNamespace("av", quietly = TRUE)) {
    stop("Package 'av' is required but not installed. Install with: install.packages('av')")
  }
  
  if (!use.am && !use.cd) {
    stop("At least one method (use.am or use.cd) must be TRUE")
  }
  
  # Initialize Python module
  creak_module <- .init_creak_python(conda.env)
  
  # Process each file
  results <- lapply(listOfFiles, function(file) {
    
    # Check file exists
    if (!file.exists(file)) {
      warning("File not found: ", file)
      return(NULL)
    }
    
    # Load audio with av (in-memory)
    audio_info <- av::av_media_info(file)
    
    if (is.null(audio_info$audio)) {
      warning("No audio stream found in: ", file)
      return(NULL)
    }
    
    # Read audio as mono
    audio_data <- av::read_audio_bin(file, channels = 1)
    sample_rate <- audio_info$audio$sample_rate
    
    # Handle beginTime/endTime if specified
    if (beginTime > 0 || endTime > 0) {
      warning("beginTime/endTime not yet implemented - processing entire file")
    }
    
    # Call Python detector with extended output
    py_result <- creak_module$detect_creak_union_extended(
      audio = audio_data,
      sample_rate = as.integer(sample_rate),
      use_am = use.am,
      use_cd = use.cd,
      cd_threshold = cd.threshold,
      use_reaper = use.reaper,
      frame_shift_ms = windowShift,
      return_features = return.features,
      return_probabilities = return.probabilities
    )
    
    # Convert to AsspDataObj with all tracks
    assp_obj <- .create_assp_creak_obj(
      py_result = py_result,
      file_path = file,
      sample_rate = sample_rate,
      use_am = use.am,
      use_cd = use.cd,
      return_features = return.features,
      return_probabilities = return.probabilities
    )
    
    # Write to file if requested
    if (toFile) {
      output_path <- .construct_creak_output_path(file, explicitExt, outputDirectory)
      .write_ssff_creak(assp_obj, output_path)
      message("Saved: ", output_path)
    }
    
    assp_obj
  })
  
  # Filter out NULLs from failed files
  results <- Filter(Negate(is.null), results)
  
  # Return single object or list
  if (length(results) == 0) {
    warning("No files were successfully processed")
    return(NULL)
  } else if (length(results) == 1) {
    return(results[[1]])
  } else {
    return(results)
  }
}


# Internal: Initialize Python creak detection module
.init_creak_python <- function(conda_env = NULL) {
  
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package 'reticulate' is required but not installed. ",
         "Install with: install.packages('reticulate')")
  }
  
  # Set conda env if specified
  if (!is.null(conda_env)) {
    reticulate::use_condaenv(conda_env, required = TRUE)
  }
  
  # Try to import module
  tryCatch({
    creak_mod <- reticulate::import("union_creak_detection_method")
    return(creak_mod)
  }, error = function(e) {
    
    # Try alternative module name (with hyphens)
    tryCatch({
      creak_mod <- reticulate::import("union-creak-detection-method")
      return(creak_mod)
    }, error = function(e2) {
      
      stop("Could not import Python module 'union_creak_detection_method'.\n",
           "Please ensure:\n",
           "  1. Python is available via reticulate\n",
           "  2. Required packages are installed (numpy, scipy, pyreaper)\n",
           "  3. Module is in PYTHONPATH\n",
           "Error: ", e$message)
    })
  })
}


# Internal: Create AsspDataObj with all tracks
.create_assp_creak_obj <- function(py_result, file_path, sample_rate,
                                    use_am, use_cd, return_features,
                                    return_probabilities) {
  
  # Extract time vector
  time <- py_result$time
  n_frames <- length(time)
  frame_rate <- 1000 / (time[2] - time[1])  # Hz
  
  # Build track list
  tracks <- list()
  track_formats <- character()
  
  # AM method tracks
  if (use_am && !is.null(py_result$am_decisions)) {
    tracks$AM_creak <- matrix(as.integer(py_result$am_decisions), ncol = 1)
    track_formats <- c(track_formats, "INT16")
    
    # F0 track
    if (!is.null(py_result$F0)) {
      tracks$F0 <- matrix(as.numeric(py_result$F0), ncol = 1)
      track_formats <- c(track_formats, "REAL32")
    }
  }
  
  # CD method tracks
  if (use_cd && !is.null(py_result$cd_decisions)) {
    tracks$CD_creak <- matrix(as.integer(py_result$cd_decisions), ncol = 1)
    track_formats <- c(track_formats, "INT16")
    
    # CD probability
    if (return_probabilities && !is.null(py_result$CD_prob)) {
      tracks$CD_prob <- matrix(as.numeric(py_result$CD_prob), ncol = 1)
      track_formats <- c(track_formats, "REAL32")
    }
  }
  
  # Union track (always included)
  tracks$union_creak <- matrix(as.integer(py_result$union_decisions), ncol = 1)
  track_formats <- c(track_formats, "INT16")
  
  # Feature tracks (if requested)
  if (return_features && !is.null(py_result$features)) {
    
    feature_names <- c('H2_H1', 'peak_prominence', 'ZCR', 'IFP', 'IPS',
                      'PwP_fall', 'PwP_rise', 'F0_feat', 'F0_mean', 
                      'energy', 'power_std', 'creak_F0')
    
    # Add static features (first 12 columns)
    for (i in 1:12) {
      feat_name <- feature_names[i]
      tracks[[feat_name]] <- matrix(as.numeric(py_result$features[, i]), ncol = 1)
      track_formats <- c(track_formats, "REAL32")
    }
    
    # Optionally add delta features (columns 13-24)
    # Commented out to reduce track count, uncomment if needed
    # for (i in 13:24) {
    #   feat_name <- paste0(feature_names[i-12], "_delta")
    #   tracks[[feat_name]] <- matrix(as.numeric(py_result$features[, i]), ncol = 1)
    #   track_formats <- c(track_formats, "REAL32")
    # }
  }
  
  # Create AsspDataObj
  obj <- new("AsspDataObj",
             trackFormats = track_formats,
             origFreq = sample_rate,
             samplingRate = frame_rate,
             startTime = time[1],
             startRecord = 0,
             endRecord = n_frames - 1,
             filePath = file_path,
             tracks = tracks)
  
  # Add metadata attributes
  attr(obj, "method") <- "Union Method"
  
  if (use_am && !is.null(py_result$antimode)) {
    attr(obj, "antimode") <- py_result$antimode
  }
  
  if (use_cd) {
    attr(obj, "cd_threshold") <- py_result$cd_threshold
  }
  
  attr(obj, "processing_date") <- Sys.time()
  
  # Add track descriptions
  track_desc <- list()
  if ("AM_creak" %in% names(tracks)) {
    track_desc$AM_creak <- "AM method creak detection (0=modal, 1=creak)"
  }
  if ("CD_creak" %in% names(tracks)) {
    track_desc$CD_creak <- "CD method creak detection (0=modal, 1=creak)"
  }
  track_desc$union_creak <- "Union of AM and CD methods (0=modal, 1=creak)"
  
  if ("F0" %in% names(tracks)) {
    track_desc$F0 <- "Fundamental frequency (Hz)"
  }
  if ("CD_prob" %in% names(tracks)) {
    track_desc$CD_prob <- "CD method creak probability (0-1)"
  }
  if ("H2_H1" %in% names(tracks)) {
    track_desc$H2_H1 <- "H2-H1: Second harmonic minus first harmonic (dB)"
  }
  if ("peak_prominence" %in% names(tracks)) {
    track_desc$peak_prominence <- "Residual peak prominence (normalized)"
  }
  if ("IFP" %in% names(tracks)) {
    track_desc$IFP <- "Intraframe periodicity (Ishi et al., 2008)"
  }
  if ("IPS" %in% names(tracks)) {
    track_desc$IPS <- "Interpulse similarity (Ishi et al., 2008)"
  }
  if ("PwP_rise" %in% names(tracks)) {
    track_desc$PwP_rise <- "Power pulse waveform rising (dB, Ishi et al., 2008)"
  }
  if ("PwP_fall" %in% names(tracks)) {
    track_desc$PwP_fall <- "Power pulse waveform falling (dB, Ishi et al., 2008)"
  }
  
  attr(obj, "track_descriptions") <- track_desc
  
  obj
}


# Internal: Construct output file path
.construct_creak_output_path <- function(input_file, ext, output_dir) {
  
  base_name <- tools::file_path_sans_ext(basename(input_file))
  
  if (is.null(output_dir)) {
    output_dir <- dirname(input_file)
  } else {
    # Create output directory if it doesn't exist
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
  }
  
  output_file <- file.path(output_dir, paste0(base_name, ".", ext))
  
  output_file
}


# Internal: Write SSFF file
.write_ssff_creak <- function(assp_obj, output_path) {
  
  # Use wrassp if available, otherwise fall back to basic write
  if (requireNamespace("wrassp", quietly = TRUE)) {
    tryCatch({
      wrassp::write.AsspDataObj(assp_obj, file = output_path)
    }, error = function(e) {
      warning("Could not write SSFF file with wrassp: ", e$message)
      .write_ssff_basic(assp_obj, output_path)
    })
  } else {
    .write_ssff_basic(assp_obj, output_path)
  }
}


# Internal: Basic SSFF writer (fallback)
.write_ssff_basic <- function(assp_obj, output_path) {
  
  # For now, write as CSV if wrassp not available
  csv_path <- sub("\\.crk$", ".csv", output_path)
  
  # Combine all tracks into dataframe
  df <- data.frame(time = seq(assp_obj@startTime,
                              by = 1/assp_obj@samplingRate,
                              length.out = assp_obj@endRecord + 1))
  
  for (track_name in names(assp_obj@tracks)) {
    df[[track_name]] <- as.vector(assp_obj@tracks[[track_name]])
  }
  
  write.csv(df, csv_path, row.names = FALSE)
  message("Note: Written as CSV (install wrassp for SSFF format): ", csv_path)
}

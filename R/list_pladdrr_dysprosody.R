#' Extract Dysprosody Prosodic Features
#'
#' Extracts 193 prosodic features using pladdrr-based dysprosody pipeline including
#' F0 analysis (MOMEL/INTSINT), intensity, formants, and spectral tilt measures
#' with Iseli-Alwan harmonic correction.
#'
#' @param listOfFiles Vector of file paths to audio files
#' @param beginTime Start time in seconds (default: 0.0)
#' @param endTime End time in seconds (default: 0.0 = full duration)
#' @param minF Minimum F0 in Hz (default: 60)
#' @param maxF Maximum F0 in Hz (default: 750)
#' @param windowShift Window shift in milliseconds (default: 10)
#' @param toFile Write output to .dyp files (default: FALSE)
#' @param explicitExt Output file extension (default: "dyp")
#' @param outputDirectory Output directory (default: NULL = same as input)
#' @param verbose Show progress (default: TRUE)
#'
#' @return If toFile=FALSE: list of data frames (one per file) with 193 features.
#'   If toFile=TRUE: invisibly returns vector of output file paths.
#'
#' @details
#' **Requires**: pladdrr >= 4.8.23
#'
#' **Features Extracted** (193 total):
#' - F0 analysis: MOMEL quadratic spline modeling, INTSINT tonal coding
#' - Intensity measures: statistics of intensity contour at INTSINT targets
#' - Spectral tilt: 7 measures including Iseli-Alwan harmonic correction (L2L1, L2cL1c, L1cLF3c, L1LF3, SLF, C1, Spectral Balance)
#' - Time-series statistics: mean, std, variation, IQR, max, min for all measures and their first differences
#' - Global measures: duration, pitch key, pitch range, INTSINT concentration
#'
#' **Performance**: ~10-12 seconds per file (optimized via batch queries)
#'
#' **Output Format**: When toFile=TRUE, writes JSON Track Format (.dyp) files
#' with all 193 features organized by time slice.
#'
#' @references
#' \insertCite{Villarubia2025}{superassp}
#'
#' @examples
#' \dontrun{
#' # Extract dysprosody features from audio
#' result <- lst_dysprosody("speech.wav", toFile = FALSE)
#' 
#' # Access specific features
#' cat("Duration:", result$Duration, "seconds\n")
#' cat("Pitch mean:", result$momelpitchtmean, "Hz\n")
#' cat("Pitch key:", result$PitchKey, "Hz\n")
#' cat("Pitch range:", result$PitchRange, "octaves\n")
#' 
#' # Batch processing with file output
#' files <- c("speech1.wav", "speech2.wav")
#' lst_dysprosody(files, toFile = TRUE, outputDirectory = "output")
#' }
#'
#' @export
lst_dysprosody <- function(listOfFiles,
                           beginTime = 0.0,
                           endTime = 0.0,
                           minF = 60,
                           maxF = 750,
                           windowShift = 10,
                           toFile = FALSE,
                           explicitExt = "dyp",
                           outputDirectory = NULL,
                           verbose = TRUE) {
  
  # Check pladdrr availability
  if (!requireNamespace("pladdrr", quietly = TRUE)) {
    stop("pladdrr package not available. Install with: install.packages('pladdrr')")
  }
  
  # Load pladdrr
  require(pladdrr, quietly = TRUE)
  
  # Validate minimum version
  pladdrr_version <- as.character(packageVersion("pladdrr"))
  if (compareVersion(pladdrr_version, "4.8.23") < 0) {
    stop("pladdrr >= 4.8.23 required (current: ", pladdrr_version, ")")
  }
  
  # Note: dysprosody functions (momel_c, intsint, prosody_measures) are 
  # defined in R/dysprosody_*.R and loaded as part of package namespace
  
  # Validate files
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)
  
  filesEx <- file.exists(listOfFiles)
  if (!all(filesEx)) {
    filesNotExist <- listOfFiles[!filesEx]
    stop("Unable to find file(s): ", paste(filesNotExist, collapse = ", "))
  }
  
  # Progress bar for multiple files
  pb <- NULL
  if (verbose && length(listOfFiles) > 1) {
    pb <- txtProgressBar(min = 0, max = length(listOfFiles), style = 3)
  }
  
  # Process each file
  results_list <- list()
  output_paths <- character(length(listOfFiles))
  
  for (i in seq_along(listOfFiles)) {
    file_path <- listOfFiles[i]
    
    # Load audio with time windowing
    sound <- tryCatch({
      if (beginTime > 0 || endTime > 0) {
        temp_sound <- pladdrr::Sound(file_path)
        duration <- temp_sound$get_duration()
        
        bt <- if (beginTime > 0) beginTime else 0
        et <- if (endTime > 0 && endTime < duration) endTime else duration
        
        temp_sound$extract_part(bt, et, window_shape = "rectangular", relative_width = 1.0)
      } else {
        pladdrr::Sound(file_path)
      }
    }, error = function(e) {
      stop("Failed to load audio file ", file_path, ": ", e$message)
    })
    
    # Run dysprosody pipeline
    result <- tryCatch({
      prosody_measures(
        sound = sound,     # Pass sound object
        minF = minF,
        maxF = maxF,
        windowShift = windowShift
      )
    }, error = function(e) {
      warning("Failed to extract dysprosody features from ", file_path, ": ", e$message)
      NULL
    })
    
    if (is.null(result)) {
      next
    }
    
    # Convert to data frame with time slice metadata
    result_df <- as.data.frame(t(unlist(result)), stringsAsFactors = FALSE)
    result_df$begin_time <- if (beginTime > 0) beginTime else 0
    result_df$end_time <- if (endTime > 0) endTime else sound$get_duration()
    
    if (toFile) {
      # Create JsonTrackObj for JSTF output
      json_obj <- create_json_track_obj(
        results = result,
        function_name = "lst_dysprosody",
        file_path = file_path,
        sample_rate = sound$.cpp$sampling_frequency,
        audio_duration = sound$get_duration(),
        beginTime = result_df$begin_time,
        endTime = result_df$end_time,
        parameters = list(
          minF = minF,
          maxF = maxF,
          windowShift = windowShift
        )
      )
      
      # Write to file
      base_name <- tools::file_path_sans_ext(basename(file_path))
      out_dir <- if (is.null(outputDirectory)) dirname(file_path) else outputDirectory
      output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))
      
      write_json_track(json_obj, output_path)
      output_paths[i] <- output_path
    } else {
      results_list[[i]] <- result_df
    }
    
    if (verbose && !is.null(pb)) {
      setTxtProgressBar(pb, i)
    }
  }
  
  if (verbose && !is.null(pb)) {
    close(pb)
  }
  
  if (toFile) {
    return(invisible(output_paths[output_paths != ""]))
  } else {
    # Remove NULL entries
    results_list <- results_list[!sapply(results_list, is.null)]
    
    # Simplify single file output
    if (length(results_list) == 1) {
      return(results_list[[1]])
    }
    return(results_list)
  }
}

# Set function attributes
attr(lst_dysprosody, "ext") <- "dyp"
attr(lst_dysprosody, "outputType") <- "JSTF"
attr(lst_dysprosody, "format") <- "JSON"

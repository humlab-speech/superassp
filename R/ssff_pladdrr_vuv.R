#' Voiced/Unvoiced Detection using pladdrr
#'
#' Detect voiced and unvoiced segments in speech using two-pass adaptive pitch
#' detection. Implements the algorithm from Al-Tamimi & Khattab (2015, 2018)
#' with speaker-adaptive pitch range estimation.
#'
#' The algorithm uses a two-pass approach:
#' 1. Initial rough pitch detection (50-800 Hz) to estimate speaker range
#' 2. Adaptive pitch range calculation from quartiles (Q1×0.75, Q3×1.5)
#' 3. Refined pitch detection with adaptive range
#' 4. PointProcess generation and voicing decision
#'
#' @param listOfFiles Character vector with path(s) to audio file(s)
#' @param beginTime Numeric. Start time in seconds (default 0)
#' @param endTime Numeric. End time in seconds (0 = end of file)
#' @param timeStep Numeric. Time step for pitch analysis in seconds (default 0.005)
#' @param initialMinPitch Numeric. Minimum pitch for first pass in Hz (default 50)
#' @param initialMaxPitch Numeric. Maximum pitch for first pass in Hz (default 800)
#' @param voicingThreshold Numeric. Voicing threshold 0-1 (default 0.45)
#' @param vuvMaxPeriod Numeric. Maximum period for VUV detection in seconds (default 0.02)
#' @param minPeriod Numeric. Minimum period for mean period calculation (default 0.0001)
#' @param maxPeriod Numeric. Maximum period for mean period calculation (default 0.02)
#' @param maxPeriodFactor Numeric. Maximum period factor (default 1.3)
#' @param windowShape Character. Window shape for extraction (default "Gaussian1")
#' @param relativeWidth Numeric. Relative width for window (default 1.0)
#' @param outputFormat Character. Output format: "textgrid" or "ssff" (default "textgrid")
#' @param toFile Logical. If TRUE, write results to file. Default TRUE.
#' @param explicitExt Character. File extension for output. Default "TextGrid" (or "vuv" for SSFF).
#' @param outputDirectory Character. Output directory path. Default NULL (use input directory).
#' @param verbose Logical. Print progress messages (default TRUE)
#'
#' @return If \code{outputFormat="textgrid"}:
#'   - If \code{toFile=TRUE}: invisibly returns path(s) to TextGrid file(s)
#'   - If \code{toFile=FALSE}: returns pladdrr TextGrid object (or list of TextGrids)
#'
#'   If \code{outputFormat="ssff"}:
#'   - If \code{toFile=TRUE}: invisibly returns path(s) to SSFF file(s)
#'   - If \code{toFile=FALSE}: returns AsspDataObj with binary voicing track (or list)
#'
#' @details
#' The VUV tier contains intervals labeled:
#' - "V" - Voiced
#' - "U" - Unvoiced
#'
#' The algorithm applies a bandpass filter (0-500 Hz) before pitch detection
#' to reduce noise and improve voicing decisions.
#'
#' Adaptive pitch range is calculated as:
#' - Minimum: Q1 × 0.75
#' - Maximum: Q3 × 1.5
#'
#' where Q1 and Q3 are the first and third quartiles of the initial pitch estimate.
#'
#' @references
#' \itemize{
#'   \item Al-Tamimi, J., & Khattab, G. (2015). Acoustic correlates of the voicing contrast in Lebanese Arabic singleton and geminate stops. JASA, 138(1), 344-360.
#'   \item Al-Tamimi, J., & Khattab, G. (2018). Acoustic cue weighting in the singleton vs geminate contrast in Lebanese Arabic: The case of fricative consonants. Journal of Phonetics, 71, 306-325.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate VUV TextGrid
#' trk_vuv("speech.wav", toFile = TRUE)
#' # Creates speech.TextGrid with VUV tier
#' 
#' # Get TextGrid object in memory
#' tg <- trk_vuv("speech.wav", toFile = FALSE)
#' 
#' # Generate binary SSFF track instead
#' result <- trk_vuv("speech.wav", outputFormat = "ssff", toFile = FALSE)
#' plot(result$voicing, type = "l")
#' 
#' # Custom pitch range
#' trk_vuv("speech.wav", initialMinPitch = 75, initialMaxPitch = 500)
#' }
trk_vuv <- function(listOfFiles,
                    beginTime = 0.0,
                    endTime = 0.0,
                    timeStep = 0.005,
                    initialMinPitch = 50,
                    initialMaxPitch = 800,
                    voicingThreshold = 0.45,
                    vuvMaxPeriod = 0.02,
                    minPeriod = 0.0001,
                    maxPeriod = 0.02,
                    maxPeriodFactor = 1.3,
                    windowShape = "Gaussian1",
                    relativeWidth = 1.0,
                    outputFormat = "textgrid",
                    toFile = TRUE,
                    explicitExt = NULL,
                    outputDirectory = NULL,
                    verbose = TRUE) {
  
  # Check pladdrr availability
  if (!pladdrr_available()) {
    stop("pladdrr package not available. Install with: install_pladdrr()")
  }
  
  # Validate output format
  outputFormat <- match.arg(outputFormat, c("textgrid", "ssff"))
  
  # Set default extension based on format
  if (is.null(explicitExt)) {
    explicitExt <- if (outputFormat == "textgrid") "TextGrid" else "vuv"
  }
  
  # Check single file restriction for toFile=FALSE
  if (length(listOfFiles) > 1 && !toFile) {
    stop("length(listOfFiles) > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.")
  }
  
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
    
    # Load sound with pladdrr
    sound <- av_load_for_pladdrr(
      file_path = file_path,
      start_time = beginTime,
      end_time = endTime,
      window_type = windowShape,
      relative_width = relativeWidth
    )
    
    # Apply bandpass filter (0-500 Hz with 20 Hz smoothing)
    sound_filtered <- sound$filter_pass_hann_band(0, 500, 20)
    
    # Two-pass adaptive pitch extraction
    pitch_result <- pladdrr::two_pass_adaptive_pitch(
      sound_filtered,
      initial_floor = initialMinPitch,
      initial_ceiling = initialMaxPitch,
      voicing_threshold = voicingThreshold,
      q1_factor = 0.75,
      q3_factor = 1.5,
      method = "cc"
    )
    
    # Create PointProcess from Sound + Pitch
    point_process_ptr <- pladdrr::to_point_process_from_sound_and_pitch(
      sound_filtered,
      pitch_result$pitch
    )
    point_process <- pladdrr::PointProcess(.xptr = point_process_ptr)
    
    # Get mean period
    mean_period <- point_process$get_mean_period(
      0, 0,
      minPeriod,
      maxPeriod,
      maxPeriodFactor
    )
    
    # Create VUV TextGrid from PointProcess
    textgrid_vuv <- point_process$to_textgrid_vuv(mean_period, vuvMaxPeriod)
    
    if (outputFormat == "textgrid") {
      # TextGrid output
      if (toFile) {
        base_name <- tools::file_path_sans_ext(basename(file_path))
        out_dir <- if (is.null(outputDirectory)) dirname(file_path) else outputDirectory
        
        if (!is.null(outputDirectory) && !dir.exists(outputDirectory)) {
          dir.create(outputDirectory, recursive = TRUE)
        }
        
        output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))
        textgrid_vuv$save(output_path)
        output_paths[i] <- output_path
        
        if (verbose && length(listOfFiles) == 1) {
          message("Wrote VUV TextGrid to: ", output_path)
          message(sprintf("  Adaptive pitch range: %.1f-%.1f Hz", 
                         pitch_result$min_pitch, pitch_result$max_pitch))
        }
      }
      
      results_list[[i]] <- textgrid_vuv
      
    } else {
      # SSFF output - convert TextGrid to binary track
      # Extract intervals from tier 1
      n_intervals <- textgrid_vuv$get_number_of_intervals(1)
      
      # Create time grid
      sample_rate <- sound$.cpp$sampling_frequency
      duration <- sound$.cpp$duration
      num_samples <- as.integer(duration * sample_rate)
      voicing <- rep(0, num_samples)
      
      # Fill in voiced segments
      for (j in 1:n_intervals) {
        label <- textgrid_vuv$get_interval_text(1, j)
        if (label == "V") {
          start_time <- textgrid_vuv$get_interval_start_time(1, j)
          end_time <- textgrid_vuv$get_interval_end_time(1, j)
          start_sample <- max(1, as.integer(start_time * sample_rate) + 1)
          end_sample <- min(num_samples, as.integer(end_time * sample_rate))
          voicing[start_sample:end_sample] <- 1
        }
      }
      
      # Downsample to requested time step
      frame_rate <- 1 / timeStep
      num_frames <- as.integer(duration * frame_rate)
      voicing_ds <- numeric(num_frames)
      
      for (j in 1:num_frames) {
        frame_time <- (j - 1) * timeStep
        sample_idx <- as.integer(frame_time * sample_rate) + 1
        if (sample_idx >= 1 && sample_idx <= num_samples) {
          voicing_ds[j] <- voicing[sample_idx]
        }
      }
      
      # Build AsspDataObj
      assp_obj <- list(
        voicing = voicing_ds
      )
      
      attr(assp_obj, "sampleRate") <- frame_rate
      attr(assp_obj, "startTime") <- 0
      attr(assp_obj, "startRecord") <- 1L
      attr(assp_obj, "endRecord") <- as.integer(num_frames)
      attr(assp_obj, "trackFormats") <- "INT16"
      attr(assp_obj, "origFreq") <- sample_rate
      class(assp_obj) <- "AsspDataObj"
      
      # Write to file if requested
      if (toFile) {
        base_name <- tools::file_path_sans_ext(basename(file_path))
        out_dir <- if (is.null(outputDirectory)) dirname(file_path) else outputDirectory
        
        if (!is.null(outputDirectory) && !dir.exists(outputDirectory)) {
          dir.create(outputDirectory, recursive = TRUE)
        }
        
        output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))
        wrassp::write.AsspDataObj(assp_obj, output_path)
        output_paths[i] <- output_path
        
        if (verbose && length(listOfFiles) == 1) {
          message("Wrote VUV binary track to: ", output_path)
        }
      }
      
      results_list[[i]] <- assp_obj
    }
    
    # Update progress bar
    if (!is.null(pb)) {
      setTxtProgressBar(pb, i)
    }
  }
  
  # Close progress bar
  if (!is.null(pb)) {
    close(pb)
  }
  
  # Return results
  if (toFile) {
    return(invisible(output_paths))
  } else {
    if (length(listOfFiles) == 1) {
      return(results_list[[1]])
    } else {
      return(results_list)
    }
  }
}

# Set function attributes (for SSFF mode)
attr(trk_vuv, "ext") <- "vuv"
attr(trk_vuv, "tracks") <- "voicing"
attr(trk_vuv, "outputType") <- "SSFF"  # or TextGrid
attr(trk_vuv, "nativeFiletypes") <- c("wav")

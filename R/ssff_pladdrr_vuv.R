#' Voiced/unvoiced segmentation via two-pass adaptive pitch detection
#'
#' Classifies each frame as voiced (V) or unvoiced (U) using a two-pass,
#' speaker-adaptive pitch detection strategy (Al-Tamimi & Khattab 2015, 2018).
#' Output is either a Praat TextGrid with a VUV interval tier or an SSFF binary
#' \code{voicing} track (0 = unvoiced, 1 = voiced).
#'
#' @inheritParams trk_acf
#' @param timeStep Numeric. Frame shift for pitch analysis in seconds. Also sets the
#'   SSFF output frame rate (1 / timeStep Hz) when \code{outputFormat = "ssff"}.
#'   Default 0.005 s (200 Hz).
#' @param initialMinPitch Numeric. Lower F0 bound for the first-pass pitch estimate
#'   in Hz. Default 50 Hz.
#' @param initialMaxPitch Numeric. Upper F0 bound for the first-pass pitch estimate
#'   in Hz. Default 800 Hz.
#' @param voicingThreshold Numeric. Minimum pitch strength for a frame to be voiced
#'   (0–1). Default 0.45.
#' @param vuvMaxPeriod Numeric. Maximum glottal period (s) for TextGrid VUV conversion.
#'   Default 0.02 s (minimum 50 Hz).
#' @param minPeriod Numeric. Minimum glottal period (s) accepted when computing mean
#'   period from the PointProcess. Default 0.0001 s.
#' @param maxPeriod Numeric. Maximum glottal period (s) accepted when computing mean
#'   period from the PointProcess. Default 0.02 s.
#' @param maxPeriodFactor Numeric. Maximum ratio between consecutive periods still
#'   treated as periodic. Default 1.3.
#' @param windowShape Character. Window shape applied to audio before loading.
#'   Default \code{"Gaussian1"}.
#' @param relativeWidth Numeric. Relative width of the window. Default 1.0.
#' @param outputFormat Character. Output format: \code{"textgrid"} (Praat TextGrid
#'   with interval tier) or \code{"ssff"} (binary AsspDataObj track). Default
#'   \code{"textgrid"}.
#' @param toFile Logical. If \code{TRUE}, write output files and return paths
#'   (invisibly). If \code{FALSE}, return the in-memory object. Default \code{TRUE}.
#' @param explicitExt Character. Output file extension. Defaults to \code{"TextGrid"}
#'   when \code{outputFormat = "textgrid"} and \code{"vuv"} when \code{"ssff"}.
#'
#' @return Depends on \code{outputFormat} and \code{toFile}:
#'   \describe{
#'     \item{\code{outputFormat = "textgrid"}, \code{toFile = FALSE}}{A pladdrr
#'       TextGrid object (or list) with one interval tier containing V/U labels.}
#'     \item{\code{outputFormat = "textgrid"}, \code{toFile = TRUE}}{Character vector
#'       of output TextGrid paths, returned invisibly.}
#'     \item{\code{outputFormat = "ssff"}, \code{toFile = FALSE}}{An \code{AsspDataObj}
#'       with track \code{voicing} (INT16, binary 0/1, n_frames x 1). Frame rate:
#'       \code{1 / timeStep} Hz (default 200 Hz).}
#'     \item{\code{outputFormat = "ssff"}, \code{toFile = TRUE}}{Character vector
#'       of output SSFF file paths, returned invisibly.}
#'   }
#'
#' @details
#' Pass 1: pitch estimated across \code{initialMinPitch}–\code{initialMaxPitch} Hz
#' on a 0–500 Hz bandpass-filtered signal. Pass 2: adaptive bounds set to
#' Q1 x 0.75 and Q3 x 1.5 of the voiced frames from pass 1. The PointProcess
#' derived from the refined pitch drives TextGrid VUV interval creation.
#'
#' @references
#' \insertCite{AlTamimi2015}{superassp}
#'
#' \insertCite{AlTamimi2018}{superassp}
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

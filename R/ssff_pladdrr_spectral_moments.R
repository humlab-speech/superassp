#' Spectral moments (CoG, SD, skewness, kurtosis)
#'
#' Extracts the four spectral moments — center of gravity, standard deviation,
#' skewness, and kurtosis — as time-series tracks via Praat's spectrogram in
#' pladdrr. Spectral moments characterize spectral shape and are widely used
#' for consonant place-of-articulation and voice quality analysis.
#'
#' @inheritParams trk_acf
#' @param windowLength Numeric. Spectrogram analysis window length in seconds.
#'   Default 0.005 s.
#' @param maximum_frequency Numeric. Highest frequency included in moment calculations
#'   in Hz. Set to 0 to use the Nyquist frequency. Default 0.
#' @param time_step Numeric. Frame shift in seconds; sets output frame rate
#'   (1 / time_step Hz). Default 0.005 s (200 Hz).
#' @param frequency_step Numeric. Frequency resolution of the spectrogram in Hz.
#'   Default 20 Hz.
#' @param power Numeric. Exponent applied to the amplitude spectrum before computing
#'   moments. Default 2 (power spectrum).
#' @param windowShape Character. Window shape for audio extraction. Default
#'   \code{"Gaussian1"}.
#' @param relativeWidth Numeric. Relative width of the extraction window. Default 1.0.
#' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
#'   paths written (invisibly). If \code{FALSE}, return an \code{AsspDataObj}.
#'   Default \code{TRUE}.
#' @param explicitExt Character. Output file extension. Default \code{"spm"}.
#'
#' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
#'   \describe{
#'     \item{\code{cog}}{REAL32, Hz, n_frames x 1. Spectral center of gravity (spectral mean).}
#'     \item{\code{sd}}{REAL32, Hz, n_frames x 1. Spectral standard deviation (spread).}
#'     \item{\code{skewness}}{REAL32, dimensionless, n_frames x 1. Spectral skewness (asymmetry).}
#'     \item{\code{kurtosis}}{REAL32, dimensionless, n_frames x 1. Spectral kurtosis (peakedness).}
#'   }
#'   Frame rate: \code{1 / time_step} Hz (default 200 Hz).
#'   If \code{toFile = TRUE}: character vector of output file paths, returned invisibly.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract spectral moments
#' result <- trk_spectral_moments("speech.wav", toFile = FALSE)
#' 
#' # Access tracks
#' plot(result$cog, type = "l", main = "Center of Gravity")
#' 
#' # Write to file
#' trk_spectral_moments("speech.wav", toFile = TRUE)
#' }
trk_spectral_moments <- function(listOfFiles,
                                  beginTime = 0.0,
                                  endTime = 0.0,
                                  windowLength = 0.005,
                                  maximum_frequency = 0.0,
                                  time_step = 0.005,
                                  frequency_step = 20.0,
                                  power = 2.0,
                                  windowShape = "Gaussian1",
                                  relativeWidth = 1.0,
                                  toFile = TRUE,
                                  explicitExt = "spm",
                                  outputDirectory = NULL,
                                  verbose = TRUE) {
  
  # Check pladdrr availability
  if (!pladdrr_available()) {
    stop("pladdrr package not available. Install with: install.packages('pladdrr')")
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
    stop("Unable to find the sound file(s) ", paste(filesNotExist, collapse = ", "))
  }
  
  # Ensure time vectors match file count
  n_files <- length(listOfFiles)
  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime) == 1) endTime <- rep(endTime, n_files)
  
  # Progress bar
  if (verbose && n_files > 1) {
    pb <- txtProgressBar(min = 0, max = n_files, style = 3)
  }
  
  # Process each file
  results_list <- list()
  output_paths <- character(n_files)
  
  for (i in seq_along(listOfFiles)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[i]
    
    tryCatch({
      # Load sound with pladdrr
      sound <- pladdrr::Sound(file_path)
      
      duration <- sound$.cpp$duration
      sampling_frequency <- sound$.cpp$sampling_frequency
      
      # Set maximum frequency to Nyquist if not specified
      max_freq <- if (maximum_frequency == 0.0) sampling_frequency / 2 else maximum_frequency
      
      # Handle time bounds
      if (bt > 0 || et > 0) {
        actual_end <- if (et > 0) et else duration
        analysis_sound <- sound$extract_part(
          from_time = bt,
          to_time = actual_end,
          window_shape = windowShape,
          relative_width = relativeWidth,
          preserve_times = TRUE
        )
      } else {
        analysis_sound <- sound
      }
      
      # Create spectrogram
      spectrogram <- analysis_sound$to_spectrogram(
        window_length = windowLength,
        max_frequency = max_freq,
        time_step = time_step,
        frequency_step = frequency_step,
        window_shape = "Gaussian"
      )
      
      # Get number of frames
      num_frames <- spectrogram$get_number_of_time_bins()
      
      if (num_frames == 0) {
        warning("No frames extracted for file: ", file_path)
        next
      }
      
      # Pre-allocate results
      times <- numeric(num_frames)
      cogs <- numeric(num_frames)
      sds <- numeric(num_frames)
      skewnesses <- numeric(num_frames)
      kurtoses <- numeric(num_frames)
      
      # Extract spectral moments for each frame
      for (frame in seq_len(num_frames)) {
        # Get time for this frame
        curr_time <- spectrogram$get_time_from_frame(frame)
        times[frame] <- curr_time
        
        # Extract spectrum slice at this time
        spectrum <- spectrogram$to_spectrum(curr_time)
        
        # Get spectral moments
        cogs[frame] <- spectrum$get_centre_of_gravity(power)
        sds[frame] <- spectrum$get_standard_deviation(power)
        skewnesses[frame] <- spectrum$get_skewness(power)
        kurtoses[frame] <- spectrum$get_kurtosis(power)
      }
      
      # Create AsspDataObj
      assp_obj <- list(
        cog = cogs,
        sd = sds,
        skewness = skewnesses,
        kurtosis = kurtoses
      )
      
      # Set attributes
      attr(assp_obj, "sampleRate") <- 1.0 / time_step
      attr(assp_obj, "startTime") <- times[1]
      attr(assp_obj, "startRecord") <- 1L
      attr(assp_obj, "endRecord") <- as.integer(num_frames)
      attr(assp_obj, "trackFormats") <- rep("REAL32", 4)
      attr(assp_obj, "origFreq") <- sampling_frequency
      class(assp_obj) <- "AsspDataObj"
      
      # Handle file output
      if (toFile) {
        base_name <- tools::file_path_sans_ext(basename(file_path))
        out_dir <- if (is.null(outputDirectory)) dirname(file_path) else outputDirectory
        
        # Create output directory if needed
        if (!dir.exists(out_dir)) {
          dir.create(out_dir, recursive = TRUE)
        }
        
        output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))
        wrassp::write.AsspDataObj(assp_obj, output_path)
        output_paths[i] <- output_path
      } else {
        results_list[[i]] <- assp_obj
      }
      
    }, error = function(e) {
      warning("Error processing file ", file_path, ": ", conditionMessage(e))
      if (!toFile) {
        results_list[[i]] <- NULL
      }
    })
    
    if (verbose && n_files > 1) {
      setTxtProgressBar(pb, i)
    }
  }
  
  if (verbose && n_files > 1) {
    close(pb)
  }
  
  # Return results
  if (toFile) {
    return(invisible(output_paths))
  } else {
    if (n_files == 1) {
      return(results_list[[1]])
    } else {
      return(results_list)
    }
  }
}


# Set function attributes
attr(trk_spectral_moments, "ext") <- "spm"
attr(trk_spectral_moments, "tracks") <- c("cog", "sd", "skewness", "kurtosis")
attr(trk_spectral_moments, "outputType") <- "SSFF"
attr(trk_spectral_moments, "nativeFiletypes") <- c("wav")

#' Cepstral Peak Prominence Smoothed (CPPS)
#'
#' Extracts time-series Cepstral Peak Prominence Smoothed (CPPS) via Praat's
#' PowerCepstrogram. CPPS quantifies voice periodicity and is a robust correlate
#' of breathiness and dysphonia. Prefer this over instantaneous CPP when temporal
#' smoothing is desired.
#'
#' @param listOfFiles Character vector of audio file paths. Any format supported by
#'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
#' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
#' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
#' @param minF Numeric. Lower quefrency bound for cepstral peak search, in Hz (as
#'   reciprocal of quefrency). Sets the minimum F0 detectable. Default 60 Hz.
#' @param maxF Numeric. Upper quefrency bound for cepstral peak search, in Hz.
#'   Sets the maximum F0 detectable. Default 333 Hz.
#' @param timeStep Numeric. Frame shift for the PowerCepstrogram in seconds.
#'   Sets output frame rate (1 / timeStep Hz). Default 0.002 s (500 Hz).
#' @param maximumFrequency Numeric. Highest frequency included in the cepstrum in Hz.
#'   Default 5000 Hz. Set to 0 for Nyquist.
#' @param preEmphFrom Numeric. Pre-emphasis onset frequency in Hz. Default 50 Hz.
#' @param windowShape Character. Window shape applied to each analysis frame.
#'   Default \code{"Hanning"}.
#' @param relativeWidth Numeric. Relative width of the analysis window. Default 1.0.
#' @param subtractTilt Logical. If \code{TRUE}, subtract the fitted spectral tilt
#'   trend before measuring peak prominence (gives CPPS rather than CPP). Default \code{TRUE}.
#' @param timeAveragingWindow Numeric. Duration of the smoothing window along the
#'   time axis in seconds. Default 0.02 s.
#' @param quefrencyAveragingWindow Numeric. Width of the smoothing window along the
#'   quefrency axis in seconds. Default 0.0005 s.
#' @param interpolation Character. Peak interpolation method: one of \code{"none"},
#'   \code{"parabolic"}, \code{"cubic"}, \code{"sinc70"}, \code{"sinc700"}.
#'   Default \code{"parabolic"}.
#' @param trendLineQuefrencyMin Numeric. Minimum quefrency (s) for trend line fitting.
#'   Default 0.001 s.
#' @param trendLineQuefrencyMax Numeric. Maximum quefrency (s) for trend line fitting.
#'   Default 0.05 s.
#' @param trendType Character. Shape of the fitted trend: \code{"straight"} or
#'   \code{"exponential decay"}. Default \code{"exponential decay"}.
#' @param fitMethod Character. Regression method: \code{"robust"}, \code{"least squares"},
#'   or \code{"robust slow"}. Default \code{"robust"}.
#' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
#'   paths written (invisibly). If \code{FALSE}, return an \code{AsspDataObj}.
#'   Default \code{TRUE}.
#' @param explicitExt Character. Output file extension. Default \code{"cps"}.
#' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
#'   writes alongside the input file.
#' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
#'
#' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
#'   \describe{
#'     \item{\code{cpp}}{REAL32, dB, n_frames x 1. Cepstral Peak Prominence Smoothed.
#'       Higher values indicate more periodic (healthier) phonation.}
#'   }
#'   Frame rate: \code{1 / timeStep} Hz (default 500 Hz).
#'   If \code{toFile = TRUE}: character vector of output file paths, returned invisibly.
#'
#' @details
#' CPPS is computed via Praat's PowerCepstrogram. Each frame's cepstral peak
#' prominence is measured relative to a fitted trend line (removing spectral tilt),
#' then smoothed over \code{timeAveragingWindow} and \code{quefrencyAveragingWindow}.
#' Typical values: 15–25 dB for normal voice; below 10 dB for breathy or dysphonic voice.
#'
#' @references
#' \insertCite{Hillenbrand.1994.10.1044/jshr.3704.769}{superassp}
#'
#' \insertCite{Hillenbrand1996}{superassp}
#'
#' \insertCite{HemanAckah2003}{superassp}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Extract CPPS from audio file
#' result <- trk_cpps("speech.wav", toFile = FALSE)
#' 
#' # Plot CPPS over time
#' plot(result$cpp, type = "l", main = "CPPS", ylab = "CPP (dB)", xlab = "Frame")
#' 
#' # Custom pitch range for female speaker
#' result <- trk_cpps("speech.wav", minF = 100, maxF = 400, toFile = FALSE)
#' 
#' # Batch process multiple files
#' trk_cpps(c("f1.wav", "f2.wav", "f3.wav"), toFile = TRUE)
#' }
trk_cpps <- function(listOfFiles,
                     beginTime = 0.0,
                     endTime = 0.0,
                     minF = 60,
                     maxF = 333,
                     timeStep = 0.002,
                     maximumFrequency = 5000,
                     preEmphFrom = 50,
                     windowShape = "Hanning",
                     relativeWidth = 1.0,
                     subtractTilt = TRUE,
                     timeAveragingWindow = 0.02,
                     quefrencyAveragingWindow = 0.0005,
                     interpolation = "parabolic",
                     trendLineQuefrencyMin = 0.001,
                     trendLineQuefrencyMax = 0.05,
                     trendType = "exponential decay",
                     fitMethod = "robust",
                     toFile = TRUE,
                     explicitExt = "cps",
                     outputDirectory = NULL,
                     verbose = TRUE) {
  
  # Check pladdrr availability
  if (!pladdrr_available()) {
    stop("pladdrr package not available. Install with: install_pladdrr()")
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
    
    duration <- sound$.cpp$duration
    sample_rate <- sound$.cpp$sampling_frequency
    
    # Create PowerCepstrogram for per-frame analysis
    power_cepstrogram <- tryCatch({
      sound$to_powercepstrogram(
        pitch_floor = minF,
        time_step = timeStep,
        maximum_frequency = maximumFrequency,
        pre_emphasis_frequency = preEmphFrom
      )
    }, error = function(e) {
      stop("Failed to create PowerCepstrogram for ", file_path, ": ", e$message)
    })
    
    # Get matrix to determine frame count and times
    pc_matrix <- power_cepstrogram$to_matrix()
    num_frames <- pc_matrix$get_number_of_rows()
    
    if (num_frames == 0) {
      warning("No frames generated for file: ", file_path)
      next
    }
    
    # Initialize output arrays
    times <- numeric(num_frames)
    cpp_arr <- numeric(num_frames)
    
    # Loop through frames - use get_cpp_at_time (pladdrr v4.8+)
    for (j in seq_len(num_frames)) {
      # Calculate frame time: xmin + (frame_number - 1) * dx
      frame_time <- pc_matrix$get_xmin() + (j - 1) * pc_matrix$get_dx()
      times[j] <- frame_time
      
      # Get CPP value at this time
      cpp <- tryCatch({
        power_cepstrogram$get_cpp_at_time(frame_time)
      }, error = function(e) NaN)
      
      cpp_arr[j] <- ifelse(is.null(cpp) || is.na(cpp), NaN, cpp)
    }
    
    # Build AsspDataObj
    assp_obj <- list(
      cpp = cpp_arr
    )
    
    # Set AsspDataObj attributes
    attr(assp_obj, "sampleRate") <- 1 / timeStep  # Convert time step to Hz
    attr(assp_obj, "startTime") <- times[1]
    attr(assp_obj, "startRecord") <- 1L
    attr(assp_obj, "endRecord") <- as.integer(num_frames)
    attr(assp_obj, "trackFormats") <- "REAL32"
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
        message("Wrote CPPS results to: ", output_path)
      }
    }
    
    results_list[[i]] <- assp_obj
    
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

# Set function attributes
attr(trk_cpps, "ext") <- "cps"
attr(trk_cpps, "tracks") <- "cpp"
attr(trk_cpps, "outputType") <- "SSFF"
attr(trk_cpps, "nativeFiletypes") <- c("wav")

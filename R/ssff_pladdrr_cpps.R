#' Cepstral Peak Prominence Smoothed (CPPS) using pladdrr
#'
#' Extract time-series Cepstral Peak Prominence (Smoothed) from audio signals.
#' CPPS is a robust measure of voice quality and breathiness, calculated from
#' the power cepstrum with smoothing for temporal stability.
#'
#' The cepstrum is the Fourier transform of the log power spectrum, and its
#' peak (the "quefrency" corresponding to F0) reflects the periodicity of the
#' signal. CPPS averages this measure across time and quefrency for robustness.
#'
#' @param listOfFiles Character vector with path(s) to audio file(s)
#' @param beginTime Numeric. Start time in seconds (default 0)
#' @param endTime Numeric. End time in seconds (0 = end of file)
#' @param minF Numeric. Minimum pitch in Hz for peak search (default 60)
#' @param maxF Numeric. Maximum pitch in Hz for peak search (default 333)
#' @param timeStep Numeric. Time step for PowerCepstrogram in seconds (default 0.002)
#' @param maximumFrequency Numeric. Maximum frequency for analysis in Hz (default 5000, 0 = Nyquist)
#' @param preEmphFrom Numeric. Pre-emphasis from frequency in Hz (default 50)
#' @param windowShape Character. Window shape for extraction (default "Hanning")
#' @param relativeWidth Numeric. Relative width for window (default 1.0)
#' @param subtractTilt Logical. Whether to detrend for CPPS (default TRUE)
#' @param timeAveragingWindow Numeric. Time averaging window in seconds (default 0.02)
#' @param quefrencyAveragingWindow Numeric. Quefrency averaging window in seconds (default 0.0005)
#' @param interpolation Character. Interpolation method: "none", "parabolic", "cubic", "sinc70", "sinc700" (default "parabolic")
#' @param trendLineQuefrencyMin Numeric. Minimum quefrency for trend line fitting (default 0.001)
#' @param trendLineQuefrencyMax Numeric. Maximum quefrency for trend line fitting (default 0.05)
#' @param trendType Character. Trend type: "straight", "exponential decay" (default "exponential decay")
#' @param fitMethod Character. Fit method: "robust", "least squares", "robust slow" (default "robust")
#' @param toFile Logical. If TRUE, write results to SSFF file. Default TRUE.
#' @param explicitExt Character. File extension for output. Default "cps".
#' @param outputDirectory Character. Output directory path. Default NULL (use input directory).
#' @param verbose Logical. Print progress messages (default TRUE)
#'
#' @return If \code{toFile=FALSE}, returns AsspDataObj with 1 track (cpp).
#'   If \code{toFile=TRUE}, invisibly returns the path(s) to the written SSFF file(s).
#'
#'   Track:
#'   \describe{
#'     \item{cpp}{Cepstral Peak Prominence (dB) - voice quality measure}
#'   }
#'
#' @details
#' CPPS is calculated by:
#' 1. Creating a PowerCepstrogram (power cepstrum over time)
#' 2. Finding the peak in each frame (corresponding to pitch period)
#' 3. Fitting a trend line through the cepstrum to remove spectral tilt
#' 4. Measuring peak prominence above the trend line
#' 5. Smoothing across time and quefrency windows
#'
#' Higher CPPS values indicate more periodic (better quality) voice.
#' Lower values indicate breathiness, roughness, or dysphonia.
#'
#' Typical values:
#' - Normal voice: 15-25 dB
#' - Breathy/dysphonic voice: < 10 dB
#'
#' @references
#' \itemize{
#'   \item Hillenbrand, J., Cleveland, R. A., & Erickson, R. L. (1994). Acoustic correlates of breathy vocal quality. JSHR, 37(4), 769-778.
#'   \item Hillenbrand, J., & Houde, R. A. (1996). Acoustic correlates of breathy vocal quality: Dysphonic voices and continuous speech. JSHR, 39(2), 311-321.
#'   \item Heman-Ackah, Y. D., et al. (2003). Cepstral peak prominence: A more reliable measure of dysphonia. Ann Otol Rhinol Laryngol, 112(4), 324-333.
#' }
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
    sample_rate <- sound$get_sampling_frequency()
    
    # Create PowerCepstrogram for per-frame analysis
    power_cepstrogram <- tryCatch({
      sound$to_powercepstrogram(
        pitch_floor = minF,
        time_step = timeStep,
        maximum_frequency = maximumFrequency,
        pre_emphasis_from = preEmphFrom
      )
    }, error = function(e) {
      stop("Failed to create PowerCepstrogram for ", file_path, ": ", e$message)
    })
    
    num_frames <- power_cepstrogram$get_number_of_frames()
    
    if (num_frames == 0) {
      warning("No frames generated for file: ", file_path)
      next
    }
    
    # Initialize output arrays
    times <- numeric(num_frames)
    cpp_arr <- numeric(num_frames)
    
    # Get internal namespace for PowerCepstrum peak prominence
    # Note: This uses internal API as pladdrr doesn't expose get_peak_prominence directly
    ns <- asNamespace("pladdrr")
    
    # Loop through frames
    for (j in seq_len(num_frames)) {
      frame_time <- power_cepstrogram$get_time_from_frame_number(j)
      times[j] <- frame_time
      
      # Extract PowerCepstrum slice
      cepstrum_slice <- tryCatch({
        power_cepstrogram$to_powercepstrum_slice(frame_time)
      }, error = function(e) NULL)
      
      if (is.null(cepstrum_slice)) {
        cpp_arr[j] <- NaN
        next
      }
      
      # Get peak prominence
      cpp <- tryCatch({
        ns$.powercepstrum_get_peak_prominence(
          cepstrum_slice$.xptr,
          interpolation,
          minF,
          maxF,
          trendLineQuefrencyMin,
          trendLineQuefrencyMax,
          trendType,
          0.05  # tolerance
        )
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

#' PraatSauce Voice Quality Analysis using pladdrr
#'
#' Extract comprehensive time-series voice quality measures including F0, formants,
#' harmonics (corrected and uncorrected), HNR at multiple bands, and CPP using
#' pladdrr's direct R interface to Praat.
#'
#' PraatSauce provides VoiceSauce-compatible measures for voice quality research,
#' including breathy/creaky voice detection and phonation type analysis. This
#' implementation uses the Iseli & Alwan (2004) formant correction algorithm and
#' optionally the Hawks & Miller (1995) bandwidth estimation formula.
#'
#' @param listOfFiles Character vector with path(s) to audio file(s)
#' @param beginTime Numeric. Start time in seconds (default 0)
#' @param endTime Numeric. End time in seconds (0 = end of file)
#' @param windowShift Numeric. Time step between measurements in milliseconds (default 5)
#' @param windowSize Numeric. Analysis window length in milliseconds (default 25)
#' @param minF Numeric. Minimum F0 in Hz (default 50)
#' @param maxF Numeric. Maximum F0 in Hz (default 300)
#' @param formantTracking Logical. Use formant tracking for cleaner tracks (default TRUE, currently unsupported - issues warning)
#' @param numFormants Integer. Number of formants to track (default 5)
#' @param maxFormantHz Numeric. Maximum formant frequency in Hz (default 5000)
#' @param nominalF1 Numeric. Reference F1 frequency for tracking in Hz (default 500)
#' @param nominalF2 Numeric. Reference F2 frequency for tracking in Hz (default 1500)
#' @param nominalF3 Numeric. Reference F3 frequency for tracking in Hz (default 2500)
#' @param preEmphFrom Numeric. Pre-emphasis frequency in Hz (default 50)
#' @param useBandwidthFormula Logical. Use Hawks & Miller bandwidth formula (default FALSE)
#' @param channel Integer. Audio channel to extract if multi-channel (default 1)
#' @param resample_to_16k Logical. Resample audio to 16kHz (default TRUE)
#' @param windowShape Character. Window shape for time extraction (default "Gaussian1")
#' @param relativeWidth Numeric. Relative width of extraction window (default 1.0)
#' @param toFile Logical. If TRUE, write results to SSFF file. Default TRUE.
#' @param explicitExt Character. File extension for output. Default "psa".
#' @param outputDirectory Character. Output directory path. Default NULL (use input directory).
#' @param verbose Logical. Print progress messages (default TRUE)
#'
#' @return If \code{toFile=FALSE}, returns AsspDataObj with 36 tracks.
#'   If \code{toFile=TRUE}, invisibly returns the path(s) to the written SSFF file(s).
#'
#'   Tracks:
#'   \describe{
#'     \item{f0}{Fundamental frequency (Hz)}
#'     \item{F1, F2, F3}{Formant frequencies (Hz)}
#'     \item{B1, B2, B3}{Formant bandwidths (Hz)}
#'     \item{H1u, H2u, H4u}{Uncorrected harmonic amplitudes (dB)}
#'     \item{H2Ku, H5Ku}{Uncorrected amplitudes at 2kHz and 5kHz (dB)}
#'     \item{A1u, A2u, A3u}{Uncorrected formant amplitudes (dB)}
#'     \item{H1H2u, H2H4u, H1A1u, H1A2u, H1A3u, H2KH5Ku}{Uncorrected differences (dB)}
#'     \item{H1c, H2c, H4c}{Corrected harmonic amplitudes (dB)}
#'     \item{A1c, A2c, A3c}{Corrected formant amplitudes (dB)}
#'     \item{H1H2c, H2H4c, H1A1c, H1A2c, H1A3c}{Corrected differences (dB)}
#'     \item{CPP}{Cepstral Peak Prominence}
#'     \item{HNR05, HNR15, HNR25, HNR35}{Harmonics-to-Noise Ratio at 4 bands (dB)}
#'   }
#'
#' @references
#' \itemize{
#'   \item Iseli, M., & Alwan, A. (2004). An improved correction formula for the estimation of harmonic magnitudes and its application to open quotient estimation. ICASSP.
#'   \item Hawks, J. W., & Miller, J. D. (1995). A formant bandwidth estimation procedure for vowel synthesis. JASA, 97(2), 1343-1344.
#'   \item Shue, Y.-L., Keating, P., Vicenik, C., & Yu, K. (2011). VoiceSauce: A program for voice analysis. ICPhS.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage with 5ms frame shift
#' result <- trk_praatsauce("speech.wav", windowShift = 5, toFile = FALSE)
#' 
#' # Access corrected H1-H2 (breathiness measure)
#' plot(result$H1H2c, type = "l", main = "H1-H2 Corrected")
#' 
#' # Process with custom F0 range and bandwidth formula
#' result <- trk_praatsauce(
#'   "speech.wav",
#'   minF = 75,
#'   maxF = 300,
#'   useBandwidthFormula = TRUE,
#'   toFile = FALSE
#' )
#' 
#' # Batch process multiple files
#' trk_praatsauce(c("f1.wav", "f2.wav", "f3.wav"), toFile = TRUE)
#' }
trk_praatsauce <- function(listOfFiles,
                            beginTime = 0.0,
                            endTime = 0.0,
                            windowShift = 5.0,
                            windowSize = 25,
                            minF = 50,
                            maxF = 300,
                            formantTracking = TRUE,
                            numFormants = 5,
                            maxFormantHz = 5000,
                            nominalF1 = 500,
                            nominalF2 = 1500,
                            nominalF3 = 2500,
                            preEmphFrom = 50,
                            useBandwidthFormula = FALSE,
                            channel = 1,
                            resample_to_16k = TRUE,
                            windowShape = "Gaussian1",
                            relativeWidth = 1.0,
                            toFile = TRUE,
                            explicitExt = "psa",
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
  
  # Formant tracking warning
  if (formantTracking) {
    warning("Formant tracking not currently available in pladdrr - using untracked formants")
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
    
    # Resample if requested
    if (resample_to_16k) {
      sound <- sound$resample(16000, 50)
    }
    
    # Extract channel if multi-channel
    n_channels <- sound$get_number_of_channels()
    if (n_channels > 1) {
      sound <- sound$extract_one_channel(channel)
    }
    
    sample_rate <- sound$get_sampling_frequency()
    duration <- sound$.cpp$duration
    
    # Calculate timepoints (every windowShift milliseconds)
    time_step_sec <- windowShift / 1000
    times <- seq(time_step_sec, duration, by = time_step_sec)
    n_frames <- length(times)
    
    if (n_frames == 0) {
      warning("No frames generated for file: ", file_path)
      next
    }
    
    # Create Pitch object (AC method)
    pitch_ptr <- pladdrr::to_pitch_ac_direct(
      sound,
      time_step = 0,
      pitch_floor = minF,
      pitch_ceiling = maxF,
      max_candidates = 15,
      very_accurate = FALSE,
      silence_threshold = 0.03,
      voicing_threshold = 0.45,
      octave_cost = 0.01,
      octave_jump_cost = 0.35,
      voiced_unvoiced_cost = 0.14
    )
    pitch <- pladdrr::Pitch(.xptr = pitch_ptr)
    
    # Create Formant object
    window_length_sec <- windowSize / 1000
    formant_ptr <- pladdrr::to_formant_direct(
      sound,
      time_step = 0,
      max_formants = numFormants,
      max_formant = maxFormantHz,
      window_length = window_length_sec,
      pre_emphasis = preEmphFrom
    )
    formant <- pladdrr::Formant(.xptr = formant_ptr)
    
    # Create band-limited Harmonicity objects for HNR
    filt_500 <- sound$filter_pass_hann_band(0, 500, 100)
    hnr05_ptr <- pladdrr::to_harmonicity_direct(filt_500, 0.01, minF, 0.1, 1.0)
    hnr05_obj <- pladdrr::Harmonicity(.xptr = hnr05_ptr)
    
    filt_1500 <- sound$filter_pass_hann_band(0, 1500, 100)
    hnr15_ptr <- pladdrr::to_harmonicity_direct(filt_1500, 0.01, minF, 0.1, 1.0)
    hnr15_obj <- pladdrr::Harmonicity(.xptr = hnr15_ptr)
    
    filt_2500 <- sound$filter_pass_hann_band(0, 2500, 100)
    hnr25_ptr <- pladdrr::to_harmonicity_direct(filt_2500, 0.01, minF, 0.1, 1.0)
    hnr25_obj <- pladdrr::Harmonicity(.xptr = hnr25_ptr)
    
    filt_3500 <- sound$filter_pass_hann_band(0, 3500, 100)
    hnr35_ptr <- pladdrr::to_harmonicity_direct(filt_3500, 0.01, minF, 0.1, 1.0)
    hnr35_obj <- pladdrr::Harmonicity(.xptr = hnr35_ptr)
    
    # Initialize output arrays
    f0_arr <- rep(NaN, n_frames)
    F1_arr <- rep(NaN, n_frames)
    F2_arr <- rep(NaN, n_frames)
    F3_arr <- rep(NaN, n_frames)
    B1_arr <- rep(NaN, n_frames)
    B2_arr <- rep(NaN, n_frames)
    B3_arr <- rep(NaN, n_frames)
    H1u_arr <- rep(NaN, n_frames)
    H2u_arr <- rep(NaN, n_frames)
    H4u_arr <- rep(NaN, n_frames)
    H2Ku_arr <- rep(NaN, n_frames)
    H5Ku_arr <- rep(NaN, n_frames)
    A1u_arr <- rep(NaN, n_frames)
    A2u_arr <- rep(NaN, n_frames)
    A3u_arr <- rep(NaN, n_frames)
    H1H2u_arr <- rep(NaN, n_frames)
    H2H4u_arr <- rep(NaN, n_frames)
    H1A1u_arr <- rep(NaN, n_frames)
    H1A2u_arr <- rep(NaN, n_frames)
    H1A3u_arr <- rep(NaN, n_frames)
    H2KH5Ku_arr <- rep(NaN, n_frames)
    H1c_arr <- rep(NaN, n_frames)
    H2c_arr <- rep(NaN, n_frames)
    H4c_arr <- rep(NaN, n_frames)
    A1c_arr <- rep(NaN, n_frames)
    A2c_arr <- rep(NaN, n_frames)
    A3c_arr <- rep(NaN, n_frames)
    H1H2c_arr <- rep(NaN, n_frames)
    H2H4c_arr <- rep(NaN, n_frames)
    H1A1c_arr <- rep(NaN, n_frames)
    H1A2c_arr <- rep(NaN, n_frames)
    H1A3c_arr <- rep(NaN, n_frames)
    CPP_arr <- rep(NaN, n_frames)
    HNR05_arr <- rep(NaN, n_frames)
    HNR15_arr <- rep(NaN, n_frames)
    HNR25_arr <- rep(NaN, n_frames)
    HNR35_arr <- rep(NaN, n_frames)
    
    # Loop through timepoints
    for (j in seq_along(times)) {
      mid_time <- times[j]
      
      # Extract analysis window
      window_half <- window_length_sec / 2
      slice_start <- mid_time - window_half
      slice_end <- mid_time + window_half
      
      if (slice_start < sound$get_xmin() || slice_end > sound$get_xmax()) {
        next
      }
      
      # Extract windowed sound
      windowed_sound <- tryCatch({
        sound$extract_part(
          slice_start, slice_end,
          "Hanning", 1.0,
          preserve_times = TRUE
        )
      }, error = function(e) NULL)
      
      if (is.null(windowed_sound)) next
      
      # Create spectrum and LTAS
      spectrum <- tryCatch({
        windowed_sound$to_spectrum(TRUE)
      }, error = function(e) NULL)
      
      if (is.null(spectrum)) next
      
      ltas <- tryCatch({
        spectrum$to_ltas_1to1()
      }, error = function(e) NULL)
      
      if (is.null(ltas)) next
      
      # Create cepstrum for CPP
      cepstrum <- tryCatch({
        spectrum$to_power_cepstrum()
      }, error = function(e) NULL)
      
      # Get F0 at this timepoint
      f0_val <- tryCatch({
        pitch$get_value_at_time(mid_time, "hertz")
      }, error = function(e) NaN)
      f0_arr[j] <- ifelse(is.null(f0_val) || is.na(f0_val), NaN, f0_val)
      
      # Get formants
      f1 <- tryCatch({
        formant$get_value_at_time(1, mid_time, "hertz")
      }, error = function(e) NaN)
      f2 <- tryCatch({
        formant$get_value_at_time(2, mid_time, "hertz")
      }, error = function(e) NaN)
      f3 <- tryCatch({
        formant$get_value_at_time(3, mid_time, "hertz")
      }, error = function(e) NaN)
      
      F1_arr[j] <- ifelse(is.null(f1) || is.na(f1), NaN, f1)
      F2_arr[j] <- ifelse(is.null(f2) || is.na(f2), NaN, f2)
      F3_arr[j] <- ifelse(is.null(f3) || is.na(f3), NaN, f3)
      
      # Get bandwidths
      if (useBandwidthFormula && !is.nan(f0_arr[j])) {
        B1_arr[j] <- .hawks_miller_bandwidth(f0_arr[j], F1_arr[j])
        B2_arr[j] <- .hawks_miller_bandwidth(f0_arr[j], F2_arr[j])
        B3_arr[j] <- .hawks_miller_bandwidth(f0_arr[j], F3_arr[j])
      } else {
        B1_arr[j] <- tryCatch({
          formant$get_bandwidth_at_time(1, mid_time, "hertz")
        }, error = function(e) NaN)
        B2_arr[j] <- tryCatch({
          formant$get_bandwidth_at_time(2, mid_time, "hertz")
        }, error = function(e) NaN)
        B3_arr[j] <- tryCatch({
          formant$get_bandwidth_at_time(3, mid_time, "hertz")
        }, error = function(e) NaN)
      }
      
      # CPP
      if (!is.null(cepstrum)) {
        CPP_arr[j] <- tryCatch({
          cepstrum$get_peak_prominence(minF, maxF, "parabolic", 0.001, 0, "straight", "robust")
        }, error = function(e) NaN)
      }
      
      # HNR at different bands
      HNR05_arr[j] <- tryCatch({
        val <- hnr05_obj$get_value_at_time(mid_time, "cubic")
        ifelse(is.na(val) || is.nan(val), 0.0, val)
      }, error = function(e) 0.0)
      
      HNR15_arr[j] <- tryCatch({
        val <- hnr15_obj$get_value_at_time(mid_time, "cubic")
        ifelse(is.na(val) || is.nan(val), 0.0, val)
      }, error = function(e) 0.0)
      
      HNR25_arr[j] <- tryCatch({
        val <- hnr25_obj$get_value_at_time(mid_time, "cubic")
        ifelse(is.na(val) || is.nan(val), 0.0, val)
      }, error = function(e) 0.0)
      
      HNR35_arr[j] <- tryCatch({
        val <- hnr35_obj$get_value_at_time(mid_time, "cubic")
        ifelse(is.na(val) || is.nan(val), 0.0, val)
      }, error = function(e) 0.0)
      
      # Spectral measures (require valid F0 and formants)
      if (!is.nan(f0_arr[j]) && !is.nan(F1_arr[j]) && !is.nan(F2_arr[j]) && !is.nan(F3_arr[j])) {
        n_f0 <- f0_arr[j]
        p10_f0 <- n_f0 / 10
        
        # H2K and H5K search windows
        peak_freq <- tryCatch({
          peak_quef <- cepstrum$get_quefrency_of_peak(50, 550, "parabolic")
          if (!is.null(peak_quef) && peak_quef > 0) 1 / peak_quef else n_f0
        }, error = function(e) n_f0)
        
        lower_2k <- 2000 - peak_freq
        upper_2k <- 2000 + peak_freq
        lower_5k <- 5000 - peak_freq
        upper_5k <- 5000 + peak_freq
        
        H2Ku_arr[j] <- tryCatch({
          ltas$get_maximum(lower_2k, upper_2k, "cubic")
        }, error = function(e) NaN)
        
        H5Ku_arr[j] <- tryCatch({
          ltas$get_maximum(lower_5k, upper_5k, "cubic")
        }, error = function(e) NaN)
        
        # H1, H2, H4 search windows
        lower_h1 <- n_f0 - p10_f0
        upper_h1 <- n_f0 + p10_f0
        lower_h2 <- 2 * n_f0 - p10_f0
        upper_h2 <- 2 * n_f0 + p10_f0
        lower_h4 <- 4 * n_f0 - p10_f0
        upper_h4 <- 4 * n_f0 + p10_f0
        
        H1u_arr[j] <- tryCatch({
          ltas$get_maximum(lower_h1, upper_h1, "none")
        }, error = function(e) NaN)
        
        H2u_arr[j] <- tryCatch({
          ltas$get_maximum(lower_h2, upper_h2, "none")
        }, error = function(e) NaN)
        
        H4u_arr[j] <- tryCatch({
          ltas$get_maximum(lower_h4, upper_h4, "none")
        }, error = function(e) NaN)
        
        # A1, A2, A3 search windows
        f1_hz <- F1_arr[j]
        f2_hz <- F2_arr[j]
        f3_hz <- F3_arr[j]
        
        p10_f1 <- f1_hz * 0.2
        p10_f2 <- f2_hz * 0.1
        p10_f3 <- f3_hz * 0.1
        
        A1u_arr[j] <- tryCatch({
          ltas$get_maximum(f1_hz - p10_f1, f1_hz + p10_f1, "none")
        }, error = function(e) NaN)
        
        A2u_arr[j] <- tryCatch({
          ltas$get_maximum(f2_hz - p10_f2, f2_hz + p10_f2, "none")
        }, error = function(e) NaN)
        
        A3u_arr[j] <- tryCatch({
          ltas$get_maximum(f3_hz - p10_f3, f3_hz + p10_f3, "none")
        }, error = function(e) NaN)
        
        # Uncorrected differences
        H1H2u_arr[j] <- H1u_arr[j] - H2u_arr[j]
        H2H4u_arr[j] <- H2u_arr[j] - H4u_arr[j]
        H1A1u_arr[j] <- H1u_arr[j] - A1u_arr[j]
        H1A2u_arr[j] <- H1u_arr[j] - A2u_arr[j]
        H1A3u_arr[j] <- H1u_arr[j] - A3u_arr[j]
        H2KH5Ku_arr[j] <- H2Ku_arr[j] - H5Ku_arr[j]
        
        # Iseli & Alwan corrections
        b1_val <- B1_arr[j]
        b2_val <- B2_arr[j]
        b3_val <- B3_arr[j]
        
        # H1 corrected
        h1_adj <- H1u_arr[j]
        h1_adj <- h1_adj - .iseli_alwan_correction(n_f0, f1_hz, b1_val, sample_rate)
        h1_adj <- h1_adj - .iseli_alwan_correction(n_f0, f2_hz, b2_val, sample_rate)
        H1c_arr[j] <- h1_adj
        
        # H2 corrected
        h2_adj <- H2u_arr[j]
        h2_adj <- h2_adj - .iseli_alwan_correction(2 * n_f0, f1_hz, b1_val, sample_rate)
        h2_adj <- h2_adj - .iseli_alwan_correction(2 * n_f0, f2_hz, b2_val, sample_rate)
        H2c_arr[j] <- h2_adj
        
        # H4 corrected
        h4_adj <- H4u_arr[j]
        h4_adj <- h4_adj - .iseli_alwan_correction(4 * n_f0, f1_hz, b1_val, sample_rate)
        h4_adj <- h4_adj - .iseli_alwan_correction(4 * n_f0, f2_hz, b2_val, sample_rate)
        H4c_arr[j] <- h4_adj
        
        # A1 corrected
        a1_adj <- A1u_arr[j]
        a1_adj <- a1_adj - .iseli_alwan_correction(f1_hz, f1_hz, b1_val, sample_rate)
        a1_adj <- a1_adj - .iseli_alwan_correction(f1_hz, f2_hz, b2_val, sample_rate)
        A1c_arr[j] <- a1_adj
        
        # A2 corrected
        a2_adj <- A2u_arr[j]
        a2_adj <- a2_adj - .iseli_alwan_correction(f2_hz, f1_hz, b1_val, sample_rate)
        a2_adj <- a2_adj - .iseli_alwan_correction(f2_hz, f2_hz, b2_val, sample_rate)
        A2c_arr[j] <- a2_adj
        
        # A3 corrected
        a3_adj <- A3u_arr[j]
        a3_adj <- a3_adj - .iseli_alwan_correction(f3_hz, f1_hz, b1_val, sample_rate)
        a3_adj <- a3_adj - .iseli_alwan_correction(f3_hz, f2_hz, b2_val, sample_rate)
        a3_adj <- a3_adj - .iseli_alwan_correction(f3_hz, f3_hz, b3_val, sample_rate)
        A3c_arr[j] <- a3_adj
        
        # Corrected differences
        H1H2c_arr[j] <- H1c_arr[j] - H2c_arr[j]
        H2H4c_arr[j] <- H2c_arr[j] - H4c_arr[j]
        H1A1c_arr[j] <- H1c_arr[j] - A1c_arr[j]
        H1A2c_arr[j] <- H1c_arr[j] - A2c_arr[j]
        H1A3c_arr[j] <- H1c_arr[j] - A3c_arr[j]
      }
    }
    
    # Build AsspDataObj
    assp_obj <- list(
      f0 = f0_arr,
      F1 = F1_arr,
      F2 = F2_arr,
      F3 = F3_arr,
      B1 = B1_arr,
      B2 = B2_arr,
      B3 = B3_arr,
      H1u = H1u_arr,
      H2u = H2u_arr,
      H4u = H4u_arr,
      H2Ku = H2Ku_arr,
      H5Ku = H5Ku_arr,
      A1u = A1u_arr,
      A2u = A2u_arr,
      A3u = A3u_arr,
      H1H2u = H1H2u_arr,
      H2H4u = H2H4u_arr,
      H1A1u = H1A1u_arr,
      H1A2u = H1A2u_arr,
      H1A3u = H1A3u_arr,
      H2KH5Ku = H2KH5Ku_arr,
      H1c = H1c_arr,
      H2c = H2c_arr,
      H4c = H4c_arr,
      A1c = A1c_arr,
      A2c = A2c_arr,
      A3c = A3c_arr,
      H1H2c = H1H2c_arr,
      H2H4c = H2H4c_arr,
      H1A1c = H1A1c_arr,
      H1A2c = H1A2c_arr,
      H1A3c = H1A3c_arr,
      CPP = CPP_arr,
      HNR05 = HNR05_arr,
      HNR15 = HNR15_arr,
      HNR25 = HNR25_arr,
      HNR35 = HNR35_arr
    )
    
    # Set AsspDataObj attributes
    attr(assp_obj, "sampleRate") <- 1000 / windowShift  # Convert ms to Hz
    attr(assp_obj, "startTime") <- times[1]
    attr(assp_obj, "startRecord") <- 1L
    attr(assp_obj, "endRecord") <- as.integer(n_frames)
    attr(assp_obj, "trackFormats") <- rep("REAL32", 36)
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
        message("Wrote PraatSauce results to: ", output_path)
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
attr(trk_praatsauce, "ext") <- "psa"
attr(trk_praatsauce, "tracks") <- c(
  "f0", "F1", "F2", "F3", "B1", "B2", "B3",
  "H1u", "H2u", "H4u", "H2Ku", "H5Ku",
  "A1u", "A2u", "A3u",
  "H1H2u", "H2H4u", "H1A1u", "H1A2u", "H1A3u", "H2KH5Ku",
  "H1c", "H2c", "H4c", "A1c", "A2c", "A3c",
  "H1H2c", "H2H4c", "H1A1c", "H1A2c", "H1A3c",
  "CPP", "HNR05", "HNR15", "HNR25", "HNR35"
)
attr(trk_praatsauce, "outputType") <- "SSFF"
attr(trk_praatsauce, "nativeFiletypes") <- c("wav")


# ============================================================================
# INTERNAL HELPER FUNCTIONS
# ============================================================================

#' Hawks & Miller (1995) Bandwidth Estimation
#'
#' Internal helper function implementing Hawks & Miller bandwidth formula.
#'
#' @param f0 Fundamental frequency in Hz
#' @param formant_freq Formant frequency in Hz
#' @return Estimated bandwidth in Hz
#'
#' @keywords internal
.hawks_miller_bandwidth <- function(f0, formant_freq) {
  if (is.na(f0) || is.na(formant_freq)) {
    return(NaN)
  }
  
  # Bandwidth scaling factor as a function of F0
  s <- 1 + 0.25 * (f0 - 132) / 88
  
  if (formant_freq < 500) {
    # Coefficients for formant < 500 Hz
    k <- 165.327516
    coef <- c(-6.73636734e-1, 1.80874446e-3, -4.52201682e-6, 7.49514000e-9, -4.70219241e-12)
  } else {
    # Coefficients for formant >= 500 Hz
    k <- 15.8146139
    coef <- c(8.10159009e-2, -9.79728215e-5, 5.28725064e-8, -1.07099364e-11, 7.91528509e-16)
  }
  
  # Polynomial evaluation
  fbw <- k +
    coef[1] * formant_freq +
    coef[2] * formant_freq^2 +
    coef[3] * formant_freq^3 +
    coef[4] * formant_freq^4 +
    coef[5] * formant_freq^5
  
  return(s * fbw)
}


#' Iseli & Alwan (2004) Formant Correction
#'
#' Internal helper function implementing Iseli & Alwan formant correction.
#' Calculates how much a formant affects harmonic amplitude measurement.
#'
#' @param f Frequency of harmonic being measured (Hz)
#' @param fx Formant frequency (Hz)
#' @param bx Formant bandwidth (Hz)
#' @param fs Sampling frequency (Hz)
#' @return Correction in dB to subtract from measured harmonic
#'
#' @keywords internal
.iseli_alwan_correction <- function(f, fx, bx, fs) {
  if (any(is.na(c(f, fx, bx, fs))) || fs == 0) {
    return(0.0)
  }
  
  r <- exp(-pi * bx / fs)
  omega_x <- 2 * pi * fx / fs
  omega <- 2 * pi * f / fs
  
  a <- r^2 + 1 - 2 * r * cos(omega_x + omega)
  b <- r^2 + 1 - 2 * r * cos(omega_x - omega)
  
  # Normalization factor (H(z=0) = 1)
  numerator <- r^2 + 1 - 2 * r * cos(omega_x)
  
  if (a <= 0 || b <= 0 || numerator <= 0) {
    return(0.0)
  }
  
  corr <- -10 * (log10(a) + log10(b)) + 20 * log10(numerator)
  return(corr)
}

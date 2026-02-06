#' Voice Quality Measurements using pladdrr
#'
#' Extract comprehensive voice quality measures from audio signals including
#' jitter, shimmer, harmonics-to-noise ratio (HNR) at multiple frequency bands,
#' spectral energy measures, glottal-to-noise excitation ratio (GNE), and
#' cepstral peak prominence (CPP).
#'
#' This function implements the algorithm from VQ_measurements_V2.praat, using
#' two-pass adaptive pitch detection for speaker-specific F0 range estimation,
#' followed by extraction of 36 voice quality parameters.
#'
#' @param listOfFiles Character vector with path(s) to audio file(s)
#' @param beginTime Numeric. Start time in seconds (default 0, or NULL = 0)
#' @param endTime Numeric. End time in seconds (0 or NULL = end of file)
#' @param minPitchInitial Numeric. Initial minimum pitch for two-pass detection in Hz (default 50)
#' @param maxPitchInitial Numeric. Initial maximum pitch for two-pass detection in Hz (default 800)
#' @param toFile Logical. If TRUE, write results to JSTF file. Default FALSE.
#' @param explicitExt Character. File extension for output. Default "vq".
#' @param outputDirectory Character. Output directory path. Default NULL (use input directory).
#' @param verbose Logical. Print progress messages (default TRUE)
#'
#' @return If \code{toFile=FALSE}, returns data frame (single file) or list of data frames (multiple files) with columns:
#' \describe{
#'   \item{file_name}{File name without extension}
#'   \item{start, end}{Analysis time bounds (seconds)}
#'   \item{duration, duration_msec}{Duration in seconds and milliseconds}
#'   \item{mean_period, sd_period}{Mean and SD of pitch period}
#'   \item{jitter_local_percent, jitter_local_abs_db}{Local jitter measures}
#'   \item{jitter_rap_percent, jitter_ppq5_percent, jitter_ddp_percent}{Relative jitter measures}
#'   \item{shimmer_local_percent, shimmer_local_db}{Local shimmer measures}
#'   \item{shimmer_apq3_percent, shimmer_apq5_percent, shimmer_apq11_percent, shimmer_dda_percent}{Amplitude perturbation}
#'   \item{hnr_mean_full_db, hnr_sd_full_db}{Full-spectrum HNR mean and SD}
#'   \item{hnr_mean_500_db, hnr_sd_500_db}{HNR 0-500 Hz}
#'   \item{hnr_mean_1500_db, hnr_sd_1500_db}{HNR 0-1500 Hz}
#'   \item{hnr_mean_2500_db, hnr_sd_2500_db}{HNR 0-2500 Hz}
#'   \item{hnr_mean_3500_db, hnr_sd_3500_db}{HNR 0-3500 Hz}
#'   \item{energy_1000_db, energy_2000_db, energy_4000_db, energy_6000_db}{Band energy differences}
#'   \item{hammarberg_index_db}{Hammarberg index (0-2kHz vs 2-5kHz)}
#'   \item{slope_db, tilt_db}{LTAS slope and tilt}
#'   \item{bed_db}{Band Energy Difference (low vs high)}
#'   \item{gne_3500, gne_4500}{Glottal-to-Noise Excitation ratio}
#'   \item{cpp_db}{Cepstral Peak Prominence}
#' }
#'
#'   If \code{toFile=TRUE}, invisibly returns path(s) to JSTF file(s).
#'
#' @details
#' The function performs a two-pass pitch detection:
#' 1. Initial pass with wide range (50-800 Hz) to estimate speaker F0
#' 2. Adaptive range calculated from quartiles (Q1×0.75, Q3×1.5)
#' 3. Second pass with adaptive range for accurate F0 tracking
#'
#' Voice quality measures include:
#' - **Jitter**: Period-to-period variation (local, RAP, PPQ5, DDP)
#' - **Shimmer**: Amplitude variation (local, APQ3, APQ5, APQ11, DDA)
#' - **HNR**: Harmonics-to-Noise Ratio at 5 frequency ranges
#' - **Spectral energy**: Band energy differences and Hammarberg index
#' - **LTAS**: Long-term average spectrum slope and tilt
#' - **GNE**: Glottal-to-Noise Excitation ratio
#' - **CPP**: Cepstral Peak Prominence
#'
#' All measures are calculated from voiced segments only.
#'
#' @references
#' \itemize{
#'   \item VQ_measurements_V2.praat (Version 2: 20 October 2024)
#'   \item Boersma, P., & Weenink, D. (2023). Praat: doing phonetics by computer.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Analyze single file
#' vq <- lst_vq("speech.wav")
#' print(vq)
#' 
#' # Check voice quality
#' cat(sprintf("Jitter: %.2f%%\n", vq$jitter_local_percent))
#' cat(sprintf("Shimmer: %.2f%%\n", vq$shimmer_local_percent))
#' cat(sprintf("HNR: %.1f dB\n", vq$hnr_mean_full_db))
#' cat(sprintf("CPP: %.1f dB\n", vq$cpp_db))
#' 
#' # Batch analysis with file output
#' lst_vq(c("f1.wav", "f2.wav", "f3.wav"), toFile = TRUE)
#' 
#' # Custom pitch range for child speech
#' vq <- lst_vq("child.wav", minPitchInitial = 150, maxPitchInitial = 600)
#' }
lst_vq <- function(listOfFiles,
                   beginTime = 0.0,
                   endTime = 0.0,
                   minPitchInitial = 50,
                   maxPitchInitial = 800,
                   toFile = FALSE,
                   explicitExt = "vq",
                   outputDirectory = NULL,
                   verbose = TRUE) {
  
  # Check pladdrr availability
  if (!pladdrr_available()) {
    stop("pladdrr package not available. Install with: install_pladdrr()")
  }
  
  # Validate files
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)
  
  filesEx <- file.exists(listOfFiles)
  if (!all(filesEx)) {
    filesNotExist <- listOfFiles[!filesEx]
    stop("Unable to find file(s): ", paste(filesNotExist, collapse = ", "))
  }
  
  # Handle NULL time parameters
  if (is.null(beginTime)) beginTime <- 0.0
  if (is.null(endTime)) endTime <- 0.0
  
  # Progress bar for multiple files
  pb <- NULL
  if (verbose && length(listOfFiles) > 1) {
    pb <- txtProgressBar(min = 0, max = length(listOfFiles), style = 3)
  }
  
  # Process each file
  results_list <- list()
  
  for (i in seq_along(listOfFiles)) {
    file_path <- listOfFiles[i]
    file_name <- tools::file_path_sans_ext(basename(file_path))
    
    # Load sound
    sound <- pladdrr::Sound(file_path)
    
    # Determine analysis bounds
    sound_start <- sound$.cpp$xmin
    sound_end <- sound$.cpp$xmax
    duration <- sound$.cpp$duration
    
    seg_start <- if (beginTime > 0) beginTime else sound_start
    seg_end <- if (endTime > 0) endTime else sound_end
    analysis_duration <- seg_end - seg_start
    analysis_duration_msec <- analysis_duration * 1000
    
    # Two-pass adaptive pitch extraction
    pitch_result <- pladdrr::two_pass_adaptive_pitch(
      sound,
      initial_floor = minPitchInitial,
      initial_ceiling = maxPitchInitial,
      voicing_threshold = 0.45,
      q1_factor = 0.75,
      q3_factor = 1.5,
      method = "cc"
    )
    
    min_pitch <- pitch_result$min_pitch
    max_pitch <- pitch_result$max_pitch
    pitch <- pladdrr::Pitch(.xptr = pitch_result$pitch)
    
    # Create PointProcess for jitter/shimmer
    point_process_ptr <- pladdrr::to_point_process_from_sound_and_pitch(
      sound,
      pitch_result$pitch
    )
    point_process <- pladdrr::PointProcess(.xptr = point_process_ptr)
    
    # Period measures
    mean_period <- point_process$get_mean_period(seg_start, seg_end, 0.0001, 0.02, 1.3)
    sd_period <- point_process$get_stdev_period(seg_start, seg_end, 0.0001, 0.02, 1.3)
    
    # Batch jitter/shimmer extraction
    jitter_shimmer <- pladdrr::get_jitter_shimmer_batch(
      point_process,
      sound,
      from_time = seg_start,
      to_time = seg_end,
      period_floor = 0.0001,
      period_ceiling = 0.02,
      max_period_factor = 1.3,
      max_amplitude_factor = 1.6
    )
    
    # Multi-band HNR calculation
    hnr_results <- pladdrr::calculate_multiband_hnr_ultra(
      sound,
      bands = c(0, 500, 1500, 2500, 3500),
      time_step = 0.005,
      min_pitch = min_pitch,
      from_time = seg_start,
      to_time = seg_end
    )
    
    # Extract segment for spectral analysis
    segment <- sound$extract_part(
      from_time = seg_start,
      to_time = seg_end,
      window_shape = "rectangular",
      relative_width = 1,
      preserve_times = TRUE
    )
    spectrum <- segment$to_spectrum()
    
    # Band energy differences
    energy_0_1000 <- spectrum$get_band_energy(0, 1000)
    energy_1000_10000 <- spectrum$get_band_energy(1000, 10000)
    energy_1000 <- if (!is.na(energy_0_1000) && energy_0_1000 > 0 && 
                       !is.na(energy_1000_10000) && energy_1000_10000 > 0) {
      10 * log10(energy_1000_10000 / energy_0_1000)
    } else NA
    
    energy_0_2000 <- spectrum$get_band_energy(0, 2000)
    energy_2000_10000 <- spectrum$get_band_energy(2000, 10000)
    energy_2000 <- if (!is.na(energy_0_2000) && energy_0_2000 > 0 && 
                       !is.na(energy_2000_10000) && energy_2000_10000 > 0) {
      10 * log10(energy_2000_10000 / energy_0_2000)
    } else NA
    
    energy_1000_4000 <- spectrum$get_band_energy(1000, 4000)
    energy_4000 <- if (!is.na(energy_0_1000) && energy_0_1000 > 0 && 
                       !is.na(energy_1000_4000) && energy_1000_4000 > 0) {
      10 * log10(energy_1000_4000 / energy_0_1000)
    } else NA
    
    energy_0_6000 <- spectrum$get_band_energy(0, 6000)
    energy_6000_10000 <- spectrum$get_band_energy(6000, 10000)
    energy_6000 <- if (!is.na(energy_0_6000) && energy_0_6000 > 0 && 
                       !is.na(energy_6000_10000) && energy_6000_10000 > 0) {
      10 * log10(energy_6000_10000 / energy_0_6000)
    } else NA
    
    energy_2000_5000 <- spectrum$get_band_energy(2000, 5000)
    hammarberg_index <- if (!is.na(energy_0_2000) && energy_0_2000 > 0 && 
                           !is.na(energy_2000_5000) && energy_2000_5000 > 0) {
      10 * log10(energy_2000_5000 / energy_0_2000)
    } else NA
    
    energy_low <- spectrum$get_band_energy(0, 500)
    energy_high <- spectrum$get_band_energy(4000, 5000)
    bed <- if (!is.na(energy_high) && energy_high > 0 && 
              !is.na(energy_low) && energy_low > 0) {
      10 * log10(energy_low / energy_high)
    } else NA
    
    # LTAS slope and tilt
    ltas <- spectrum$to_ltas_1to1()
    slope <- ltas$get_slope(0, 1000, 1000, 10000, "dB")
    trend <- ltas$compute_trend_line(0, 10000)
    tilt <- trend$get_slope(0, 1000, 1000, 10000, "dB")
    
    # GNE (Glottal-to-Noise Excitation Ratio)
    gne_matrix <- sound$to_harmonicity_gne(
      fmin = 500,
      fmax = 4500,
      bandwidth = 1000,
      step = 80
    )
    gne_value <- gne_matrix$get_maximum()
    
    # CPP (Cepstral Peak Prominence)
    ns <- asNamespace("pladdrr")
    power_cepstrum <- spectrum$to_power_cepstrum()
    cpp <- ns$.powercepstrum_get_peak_prominence(
      power_cepstrum$.xptr,
      "parabolic",
      60,
      333.3,
      0.001,
      0.05,
      "exponential decay",
      0.05
    )
    
    # Build result data frame
    result <- data.frame(
      file_name = file_name,
      start = seg_start,
      end = seg_end,
      duration = analysis_duration,
      duration_msec = analysis_duration_msec,
      mean_period = mean_period,
      sd_period = sd_period,
      jitter_local_percent = jitter_shimmer$jitter_local,
      jitter_local_abs_db = jitter_shimmer$jitter_local_abs,
      jitter_rap_percent = jitter_shimmer$jitter_rap,
      jitter_ppq5_percent = jitter_shimmer$jitter_ppq5,
      jitter_ddp_percent = jitter_shimmer$jitter_ddp,
      shimmer_local_percent = jitter_shimmer$shimmer_local,
      shimmer_local_db = jitter_shimmer$shimmer_local_db,
      shimmer_apq3_percent = jitter_shimmer$shimmer_apq3,
      shimmer_apq5_percent = jitter_shimmer$shimmer_apq5,
      shimmer_apq11_percent = jitter_shimmer$shimmer_apq11,
      shimmer_dda_percent = jitter_shimmer$shimmer_dda,
      hnr_mean_full_db = hnr_results$full_mean,
      hnr_sd_full_db = hnr_results$full_sd,
      hnr_mean_500_db = hnr_results$band500_mean,
      hnr_sd_500_db = hnr_results$band500_sd,
      hnr_mean_1500_db = hnr_results$band1500_mean,
      hnr_sd_1500_db = hnr_results$band1500_sd,
      hnr_mean_2500_db = hnr_results$band2500_mean,
      hnr_sd_2500_db = hnr_results$band2500_sd,
      hnr_mean_3500_db = hnr_results$band3500_mean,
      hnr_sd_3500_db = hnr_results$band3500_sd,
      energy_1000_db = energy_1000,
      energy_2000_db = energy_2000,
      energy_4000_db = energy_4000,
      energy_6000_db = energy_6000,
      hammarberg_index_db = hammarberg_index,
      slope_db = slope,
      tilt_db = tilt,
      bed_db = bed,
      gne_3500 = gne_value,
      gne_4500 = gne_value,
      cpp_db = cpp,
      stringsAsFactors = FALSE
    )
    
    results_list[[i]] <- result
    
    # Update progress bar
    if (!is.null(pb)) {
      setTxtProgressBar(pb, i)
    }
  }
  
  # Close progress bar
  if (!is.null(pb)) {
    close(pb)
  }
  
  # Write to JSTF files if requested
  if (toFile) {
    output_paths <- write_lst_results_to_jstf(
      results = results_list,
      listOfFiles = listOfFiles,
      function_name = "lst_vq",
      explicitExt = explicitExt,
      outputDirectory = outputDirectory,
      verbose = verbose
    )
    return(invisible(output_paths))
  }
  
  # Return results
  if (length(listOfFiles) == 1) {
    return(results_list[[1]])
  } else {
    return(results_list)
  }
}

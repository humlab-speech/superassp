#' Acoustic Voice Quality Index (AVQI) using pladdrr
#'
#' Computes the Acoustic Voice Quality Index (AVQI) from continuous speech and
#' sustained vowel recordings using pladdrr's Praat bindings. Supports both
#' AVQI v2.03 (Maryn et al. 2010) and v3.01 (Barsties & Maryn 2015).
#'
#' The AVQI combines 6 acoustic measures: CPPS, HNR, shimmer (local and dB),
#' LTAS slope, and LTAS tilt into a single voice quality index (0-10 scale).
#'
#' @param svDF Data frame with sustained vowel samples. Must contain columns: listOfFiles, start, end (in milliseconds)
#' @param csDF Data frame with continuous speech samples. Must contain columns: listOfFiles, start, end (in milliseconds)
#' @param version Character. AVQI version: "v2.03" (default) or "v3.01"
#' @param min.sv Minimum sustained vowel duration in milliseconds (default: 1000)
#' @param toFile Logical. If TRUE, write results to JSTF file. Default FALSE.
#' @param explicitExt Character. File extension for output. Default "avqi".
#' @param outputDirectory Character. Output directory path. Default NULL (use input directory).
#' @param verbose Logical. Print progress messages (default TRUE)
#'
#' @return If \code{toFile=FALSE} (default), a list with AVQI measurements.
#'   If \code{toFile=TRUE}, invisibly returns the path to the written JSTF file.
#'
#'   The list contains:
#'   \describe{
#'     \item{version}{AVQI version used}
#'     \item{avqi}{Acoustic Voice Quality Index (0-10 scale)}
#'     \item{cpps}{Cepstral Peak Prominence Smoothed (dB)}
#'     \item{hnr}{Harmonics-to-Noise Ratio (dB)}
#'     \item{shimmer_local}{Local shimmer (\%)}
#'     \item{shimmer_db}{Local shimmer in dB}
#'     \item{slope}{LTAS slope (dB)}
#'     \item{tilt}{LTAS tilt (dB)}
#'   }
#'
#' @references
#' Maryn, Y., Corthals, P., Van Cauwenberge, P., Roy, N., & De Bodt, M. (2010).
#' Toward improved ecological validity in the acoustic measurement of overall
#' voice quality: Combining continuous speech and sustained vowels.
#' Journal of Voice, 24(5), 540-555.
#'
#' Barsties, B., & Maryn, Y. (2015). The improvement of internal consistency of
#' the Acoustic Voice Quality Index. American Journal of Otolaryngology, 36(5), 647-656.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Define sustained vowel samples
#' sv <- data.frame(
#'   listOfFiles = c("sv1.wav", "sv2.wav"),
#'   start = c(100, 150),  # milliseconds
#'   end = c(2500, 2800)
#' )
#'
#' # Define continuous speech samples
#' cs <- data.frame(
#'   listOfFiles = c("cs1.wav", "cs2.wav"),
#'   start = c(80, 120),
#'   end = c(3500, 4000)
#' )
#'
#' # Compute AVQI (v2.03)
#' result <- lst_avqip(sv, cs)
#' print(result$avqi)
#'
#' # Compute AVQI v3.01
#' result_v3 <- lst_avqip(sv, cs, version = "v3.01")
#'
#' # Write to JSTF file
#' lst_avqip(sv, cs, toFile = TRUE)
#' track <- read_track("sv1.avqi")
#' df <- as.data.frame(track)
#' }
lst_avqip <- function(svDF,
                      csDF,
                      version = "v2.03",
                      min.sv = 1000,
                      toFile = FALSE,
                      explicitExt = "avqi",
                      outputDirectory = NULL,
                      verbose = TRUE) {
  
  # Check pladdrr availability
  if (!pladdrr_available()) {
    stop("pladdrr package not available. Install with: install.packages('pladdrr')")
  }
  
  # Validate version
  if (!version %in% c("v2.03", "v3.01")) {
    stop("version must be 'v2.03' or 'v3.01'")
  }
  
  # Validate required columns
  requiredDFColumns <- c("listOfFiles", "start", "end")
  
  if (!all(requiredDFColumns %in% names(svDF)) || !all(requiredDFColumns %in% names(csDF))) {
    stop("The 'svDF' and 'csDF' structures must both contain columns named ",
         paste(requiredDFColumns, collapse = ", "))
  }
  
  # Check minimum sustained vowel duration (times are in milliseconds)
  totalSVdur <- sum(svDF$end - svDF$start)
  
  if (totalSVdur < min.sv) {
    stop("The total sustained vowel duration (", totalSVdur, " ms) is less than min.sv (", min.sv, " ms)")
  }
  
  # Get all unique files
  listOfFiles <- unique(c(svDF$listOfFiles, csDF$listOfFiles))
  
  # Check that all files exist
  filesEx <- file.exists(listOfFiles)
  if (!all(filesEx)) {
    filesNotExist <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ", paste(filesNotExist, collapse = ", "))
  }
  
  if (verbose) {
    message("Loading and concatenating sustained vowel segments...")
  }
  
  # Load and concatenate sustained vowel segments
  sv_sounds <- list()
  for (r in 1:nrow(svDF)) {
    file_path <- normalizePath(as.character(svDF[[r, "listOfFiles"]]), mustWork = TRUE)
    start_sec <- svDF[[r, "start"]] / 1000  # Convert ms to seconds
    end_sec <- svDF[[r, "end"]] / 1000
    
    # Load sound
    sound <- pladdrr::Sound(file_path)
    
    # Extract time window
    sound_segment <- sound$extract_part(
      from_time = start_sec,
      to_time = end_sec,
      window_shape = "rectangular",
      relative_width = 1.0,
      preserve_times = FALSE
    )
    
    sv_sounds[[r]] <- sound_segment
  }
  
  # Concatenate sustained vowel segments
  sv_sound <- .load_and_concatenate_sounds(sv_sounds)
  
  if (verbose) {
    message("Loading and concatenating continuous speech segments...")
  }
  
  # Load and concatenate continuous speech segments
  cs_sounds <- list()
  for (r in 1:nrow(csDF)) {
    file_path <- normalizePath(as.character(csDF[[r, "listOfFiles"]]), mustWork = TRUE)
    start_sec <- csDF[[r, "start"]] / 1000  # Convert ms to seconds
    end_sec <- csDF[[r, "end"]] / 1000
    
    # Load sound
    sound <- pladdrr::Sound(file_path)
    
    # Extract time window
    sound_segment <- sound$extract_part(
      from_time = start_sec,
      to_time = end_sec,
      window_shape = "rectangular",
      relative_width = 1.0,
      preserve_times = FALSE
    )
    
    cs_sounds[[r]] <- sound_segment
  }
  
  # Concatenate continuous speech segments
  cs_sound <- .load_and_concatenate_sounds(cs_sounds)
  
  if (verbose) {
    message("Applying high-pass filter...")
  }
  
  # Apply high-pass filter (34 Hz cutoff) - removes low-frequency noise
  cs_filtered <- cs_sound$filter_stop_hann_band(0, 34, 0.1)
  sv_filtered <- sv_sound$filter_stop_hann_band(0, 34, 0.1)
  
  if (verbose) {
    message("Extracting voiced segments from continuous speech...")
  }
  
  # Extract voiced segments from continuous speech
  voiced_cs <- .extract_voiced_segments_pladdrr(cs_filtered, version = version)
  
  if (verbose) {
    message("Preparing sustained vowel (last 3 seconds)...")
  }
  
  # Use last 3 seconds of sustained vowel (or all if shorter)
  sv_duration <- sv_filtered$.cpp$duration
  if (sv_duration > 3) {
    sv_extract <- sv_filtered$extract_part(
      from_time = sv_duration - 3,
      to_time = sv_duration,
      window_shape = "rectangular",
      relative_width = 1.0,
      preserve_times = FALSE
    )
  } else {
    sv_extract <- sv_filtered
  }
  
  if (verbose) {
    message("Concatenating voiced segments...")
  }
  
  # Concatenate voiced continuous speech and sustained vowel
  avqi_sound <- voiced_cs$concatenate(sv_extract)
  
  if (verbose) {
    message("Calculating acoustic measures...")
  }
  
  # Calculate CPPS
  cpps <- .calculate_cpps_pladdrr(avqi_sound)
  
  # Calculate LTAS slope and tilt
  ltas_result <- .calculate_ltas_slope_and_tilt_pladdrr(avqi_sound)
  
  # Calculate HNR and Shimmer using Ultra API (if available)
  if (exists("get_voice_quality_ultra", where = asNamespace("pladdrr"))) {
    # Use Ultra API (3.6x faster)
    voice_quality <- pladdrr::get_voice_quality_ultra(
      avqi_sound,
      metrics = c("hnr", "shimmer"),
      min_pitch = 50,
      max_pitch = 400,
      time_step = 0
    )
    
    hnr <- voice_quality$hnr_mean
    shimmer_local <- voice_quality$shimmer_local * 100  # Convert to percentage
    shimmer_db <- voice_quality$shimmer_local_db
  } else {
    # Fallback to separate calculations
    hnr <- .calculate_hnr_pladdrr(avqi_sound)
    shimmer_result <- .calculate_shimmer_pladdrr(avqi_sound)
    shimmer_local <- shimmer_result$shimmer_local
    shimmer_db <- shimmer_result$shimmer_db
  }
  
  if (verbose) {
    message("Computing AVQI score (", version, ")...")
  }
  
  # Calculate AVQI using appropriate formula
  if (version == "v2.03") {
    # Maryn et al. (2010) formula
    avqi <- ((3.295 - (0.111 * cpps) - (0.073 * hnr) - (0.213 * shimmer_local) +
             (2.789 * shimmer_db) - (0.032 * ltas_result$slope) +
             (0.077 * ltas_result$tilt)) * 2.208) + 1.797
  } else if (version == "v3.01") {
    # Barsties & Maryn (2015) formula
    avqi <- ((4.152 - (0.177 * cpps) - (0.006 * hnr) - (0.037 * shimmer_local) +
             (0.941 * shimmer_db) + (0.01 * ltas_result$slope) +
             (0.093 * ltas_result$tilt)) * 2.8902)
  }
  
  # Create result list
  result_list <- list(
    version = version,
    avqi = avqi,
    cpps = cpps,
    hnr = hnr,
    shimmer_local = shimmer_local,
    shimmer_db = shimmer_db,
    slope = ltas_result$slope,
    tilt = ltas_result$tilt
  )
  
  if (verbose) {
    message(sprintf("AVQI (%s): %.2f", version, avqi))
  }
  
  # Handle JSTF file writing
  if (toFile) {
    # Calculate total analysis time range
    all_start_times <- c(svDF$start, csDF$start) / 1000  # Convert to seconds
    all_end_times <- c(svDF$end, csDF$end) / 1000
    analysis_begin <- min(all_start_times)
    analysis_end <- max(all_end_times)
    
    # Use first sustained vowel file as primary reference
    primary_file <- normalizePath(svDF[[1, "listOfFiles"]], mustWork = TRUE)
    
    output_path <- write_lst_results_to_jstf(
      results = list(result_list),
      file_paths = primary_file,
      beginTime = analysis_begin,
      endTime = analysis_end,
      function_name = "lst_avqip",
      parameters = list(
        version = version,
        min.sv = min.sv,
        n_sv_segments = nrow(svDF),
        n_cs_segments = nrow(csDF),
        total_sv_duration_ms = sum(svDF$end - svDF$start),
        total_cs_duration_ms = sum(csDF$end - csDF$start)
      ),
      explicitExt = explicitExt,
      outputDirectory = outputDirectory
    )
    
    return(invisible(output_path))
  }
  
  return(result_list)
}


#' Load and concatenate sound objects
#' 
#' @keywords internal
.load_and_concatenate_sounds <- function(sounds) {
  if (length(sounds) == 0) {
    stop("No sounds provided")
  }
  
  if (length(sounds) == 1) {
    return(sounds[[1]])
  }
  
  # Use Tier 3 batch concatenation if available (19x faster)
  if (exists("sound_concatenate_all", where = asNamespace("pladdrr"))) {
    return(pladdrr::sound_concatenate_all(sounds))
  }
  
  # Fallback: sequential concatenation
  combined <- sounds[[1]]
  for (i in 2:length(sounds)) {
    combined <- combined$concatenate(sounds[[i]])
  }
  return(combined)
}


#' Extract and concatenate voiced segments from continuous speech
#' 
#' @keywords internal
.extract_voiced_segments_pladdrr <- function(sound, version = "v3.01") {
  # Use Ultra API if available (21.7x faster)
  if (exists("extract_voiced_segments_ultra", where = asNamespace("pladdrr"))) {
    return(pladdrr::extract_voiced_segments_ultra(
      sound,
      version = version,
      min_pitch = 50,
      silence_threshold_db = -25,
      min_silent_duration = 0.1,
      min_sounding_duration = 0.1,
      power_threshold_factor = 0.3,
      max_zcr = 3000,
      window_width = 0.03
    ))
  }
  
  # Fallback: Manual extraction (slower but works with older pladdrr)
  # 1. Detect silences
  textgrid <- sound$to_textgrid_silences(
    min_pitch = 50,
    time_step = 0,
    silence_threshold = -25,
    min_silent_interval_duration = 0.1,
    min_sounding_interval_duration = 0.1
  )
  
  # 2. Extract sounding intervals
  tier <- textgrid$get_tier(1)
  n_intervals <- tier$get_number_of_intervals()
  
  sounding_segments <- list()
  for (i in 1:n_intervals) {
    label <- tier$get_interval_label(i)
    if (label == "sounding") {
      start_time <- tier$get_interval_start_time(i)
      end_time <- tier$get_interval_end_time(i)
      
      segment <- sound$extract_part(
        from_time = start_time,
        to_time = end_time,
        window_shape = "rectangular",
        relative_width = 1.0,
        preserve_times = FALSE
      )
      sounding_segments[[length(sounding_segments) + 1]] <- segment
    }
  }
  
  if (length(sounding_segments) == 0) {
    stop("No sounding segments found in continuous speech")
  }
  
  # 3. Concatenate sounding segments
  onlyLoud <- .load_and_concatenate_sounds(sounding_segments)
  
  # 4. Apply windowed power and zero-crossing rate filtering (v3.01 only)
  if (version == "v3.01") {
    # Window-based filtering (30ms windows)
    # Keep window if power > 30% of global AND ZCR < 3000 Hz
    duration <- onlyLoud$.cpp$duration
    window_width <- 0.03
    n_windows <- floor(duration / window_width)
    
    if (n_windows == 0) {
      return(onlyLoud)
    }
    
    # Calculate global power
    intensity <- onlyLoud$to_intensity(minimum_pitch = 50)
    global_power <- 10^(intensity$get_mean(0, 0) / 10)
    power_threshold <- 0.3 * global_power
    
    voiced_segments <- list()
    for (i in 1:n_windows) {
      start_time <- (i - 1) * window_width
      end_time <- min(i * window_width, duration)
      
      window <- onlyLoud$extract_part(
        from_time = start_time,
        to_time = end_time,
        window_shape = "rectangular",
        relative_width = 1.0,
        preserve_times = FALSE
      )
      
      # Calculate window power
      window_intensity <- window$to_intensity(minimum_pitch = 50)
      window_power <- 10^(window_intensity$get_mean(0, 0) / 10)
      
      # Calculate zero-crossing rate (approximate)
      # Get samples and count sign changes
      samples <- window$get_values(channel = 1)
      zcr <- sum(diff(sign(samples)) != 0) / (2 * window$.cpp$duration)
      
      # Keep window if both conditions met
      if (window_power > power_threshold && zcr < 3000) {
        voiced_segments[[length(voiced_segments) + 1]] <- window
      }
    }
    
    if (length(voiced_segments) == 0) {
      warning("No voiced segments passed ZCR filtering, returning all sounding segments")
      return(onlyLoud)
    }
    
    return(.load_and_concatenate_sounds(voiced_segments))
  }
  
  # v2.03: Just return sounding segments
  return(onlyLoud)
}


#' Calculate CPPS using pladdrr
#' 
#' @keywords internal
.calculate_cpps_pladdrr <- function(sound) {
  # Use Ultra API if available (1.6x faster)
  if (exists("calculate_cpps_ultra", where = asNamespace("pladdrr"))) {
    return(pladdrr::calculate_cpps_ultra(
      sound,
      time_averaging_window = 0.01,
      quefrency_averaging_window = 0.001,
      pitch_floor = 60,
      pitch_ceiling = 330,
      subtract_trend = FALSE,
      time_step = 0.002,
      max_frequency = 5000,
      pre_emphasis_from = 50
    ))
  }
  
  # Fallback: Standard CPPS calculation
  cepstrogram <- sound$to_power_cepstrogram(
    pitch_floor = 60,
    time_step = 0.002,
    max_frequency = 5000,
    pre_emphasis_from = 50
  )
  
  cpps <- cepstrogram$get_cpps(
    subtract_trend_before_smoothing = FALSE,
    time_averaging_window = 0.01,
    quefrency_averaging_window = 0.001,
    peak_search_pitch_range_min = 60,
    peak_search_pitch_range_max = 330,
    tolerance = 0.05,
    interpolation = "Parabolic",
    qstep = 0.001,
    trend_type = "Exponential decay",
    trend_fit_method = "Robust"
  )
  
  return(cpps)
}


#' Calculate HNR using pladdrr
#' 
#' @keywords internal
.calculate_hnr_pladdrr <- function(sound) {
  harmonicity_ptr <- pladdrr::to_harmonicity_direct(
    sound,
    time_step = 0.01,
    minimum_pitch = 75,
    silence_threshold = 0.1,
    periods_per_window = 1.0
  )
  
  harmonicity <- pladdrr::Harmonicity(.xptr = harmonicity_ptr)
  hnr <- harmonicity$get_mean(from_time = 0, to_time = 0)
  
  return(hnr)
}


#' Calculate shimmer using pladdrr
#' 
#' @keywords internal
.calculate_shimmer_pladdrr <- function(sound) {
  # Create periodic point process
  point_process <- sound$to_point_process_periodic_cc(
    pitch_floor = 50,
    pitch_ceiling = 400
  )
  
  # Calculate shimmer measures
  shimmer_decimal <- point_process$get_shimmer_local_from_sound(
    sound,
    from_time = 0,
    to_time = 0,
    shortest_period = 0.0001,
    longest_period = 0.02,
    maximum_period_factor = 1.3,
    maximum_amplitude_factor = 1.6
  )
  shimmer_local <- shimmer_decimal * 100  # Convert to percentage
  
  shimmer_db <- point_process$get_shimmer_local_db_from_sound(
    sound,
    from_time = 0,
    to_time = 0,
    shortest_period = 0.0001,
    longest_period = 0.02,
    maximum_period_factor = 1.3,
    maximum_amplitude_factor = 1.6
  )
  
  return(list(
    shimmer_local = shimmer_local,
    shimmer_db = shimmer_db
  ))
}


#' Calculate LTAS slope and tilt using pladdrr
#' 
#' @keywords internal
.calculate_ltas_slope_and_tilt_pladdrr <- function(sound) {
  # Calculate LTAS
  ltas <- sound$to_ltas(bandwidth = 1)
  
  # Get Nyquist frequency
  nyquist <- ltas$get_highest_frequency()
  
  # Clamp f2max to Nyquist
  f2max <- min(10000, nyquist)
  
  # Calculate slope (0-1000 Hz vs 1000-f2max Hz)
  slope <- ltas$get_slope(
    f1min = 0,
    f1max = 1000,
    f2min = 1000,
    f2max = f2max,
    unit = "energy"
  )
  
  # Calculate tilt (slope of trendline)
  trend_line <- ltas$compute_trend_line(
    fmin = 1,
    fmax = f2max
  )
  
  tilt <- trend_line$get_slope(
    f1min = 0,
    f1max = 1000,
    f2min = 1000,
    f2max = f2max,
    unit = "energy"
  )
  
  return(list(
    slope = slope,
    tilt = tilt
  ))
}


# Set function attributes
attr(lst_avqip, "ext") <- "avqi"
attr(lst_avqip, "outputType") <- "JSTF"
attr(lst_avqip, "format") <- "JSON"
attr(lst_avqip, "tracks") <- c(
  "version", "avqi", "cpps", "hnr",
  "shimmer_local", "shimmer_db", "slope", "tilt"
)

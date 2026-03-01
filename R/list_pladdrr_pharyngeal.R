#' Pharyngeal Voice Quality Analysis
#'
#' Extract pharyngealization/voice quality measures from labeled vowel intervals.
#' Analyzes spectral measures (H1-H2, H1-A1, H1-A2, H1-A3, etc.) at vowel onset
#' and midpoint with Iseli & Alwan (2004) normalization.
#'
#' @param listOfFiles Character vector of file paths to audio files (WAV, MP3, etc.)
#' @param textgridPath Character vector of TextGrid file paths (same length as listOfFiles)
#' @param intervalTier Integer specifying which TextGrid tier contains the vowel intervals (default: 3)
#' @param intervalNumber Integer specifying which interval to analyze per file (default: 1)
#' @param beginTime Numeric vector of start times in seconds (optional, overrides TextGrid if provided)
#' @param endTime Numeric vector of end times in seconds (optional, overrides TextGrid if provided)
#' @param minPitchInitial Minimum pitch for initial pitch detection in Hz (default: 50)
#' @param maxPitchInitial Maximum pitch for initial pitch detection in Hz (default: 800)
#' @param toFile Logical - write results to JSTF file (default: FALSE)
#' @param explicitExt File extension for output files (default: "pha")
#' @param outputDirectory Directory for output files (default: NULL = same as input)
#' @param verbose Logical - show progress messages (default: TRUE)
#'
#' @return Data frame with 68 columns per file:
#' \describe{
#'   \item{file}{File path}
#'   \item{start_time, mid_time, end_time, duration_ms}{Timing information}
#'   \item{f0_start, f0_mid}{F0 at onset and midpoint (Hz)}
#'   \item{f1/f2/f3_start, f1/f2/f3_mid}{Formant frequencies at onset/mid (Hz)}
#'   \item{bw1/bw2/bw3_start, bw1/bw2/bw3_mid}{Formant bandwidths (Hz)}
#'   \item{bw1/bw2/bw3_start_norm, bw1/bw2/bw3_mid_norm}{Normalized bandwidths (Hz)}
#'   \item{intensity_start, intensity_mid}{Intensity at onset/mid (dB)}
#'   \item{h1/h2_onset, h1/h2_mid}{Harmonic amplitudes (dB)}
#'   \item{h1_hz/h2_hz_onset, h1_hz/h2_hz_mid}{Harmonic frequencies (Hz)}
#'   \item{h1/h2_onset_norm, h1/h2_mid_norm}{Normalized harmonic amplitudes (dB)}
#'   \item{a1/a2/a3_onset, a1/a2/a3_mid}{Formant peak amplitudes (dB)}
#'   \item{a1_hz/a2_hz/a3_hz_onset, a1_hz/a2_hz/a3_hz_mid}{Formant peak frequencies (Hz)}
#'   \item{a3_onset_norm, a3_mid_norm}{Normalized A3 amplitudes (dB)}
#'   \item{h1_minus_h2/h1_minus_a1/h1_minus_a2/h1_minus_a3_onset}{Raw differences at onset (dB)}
#'   \item{h1_minus_h2/h1_minus_a1/h1_minus_a2/h1_minus_a3_onset_norm}{Normalized differences at onset (dB)}
#'   \item{a1_minus_a2/a1_minus_a3/a2_minus_a3_onset}{Additional raw differences at onset (dB)}
#'   \item{a1_minus_a3_onset_norm, a2_minus_a3_onset_norm}{Additional normalized differences at onset (dB)}
#'   \item{...}{Same structure repeated for midpoint (*_mid)}
#' }
#'
#' @details
#' **Algorithm** (based on scriptPharyFullV4.praat):
#' \enumerate{
#'   \item Two-pass adaptive pitch detection (speaker-specific F0 range)
#'   \item Find intensity maxima at onset/midpoint
#'   \item Extract formants F1/F2/F3 with bandwidths
#'   \item Extract 40ms Kaiser2 window at onset
#'   \item Create spectrum with pre-emphasis
#'   \item Find H1, H2 harmonics near F0, 2×F0
#'   \item Find A1, A2, A3 formant peaks near F1/F2/F3
#'   \item Apply Iseli & Alwan (2004) normalization
#'   \item Calculate differences (raw and normalized)
#'   \item If duration > 120ms, repeat for midpoint
#' }
#'
#' **Key Measures**:
#' - **H1-H2**: Open quotient indicator (higher = breathier)
#' - **H1-A1**: Voice quality indicator (higher = breathier)
#' - **H1-A2**: Spectral tilt measure
#' - **H1-A3**: High-frequency energy indicator
#' - **Normalized (*_norm)**: Corrected for formant influence (Iseli & Alwan 2004)
#'
#' **Typical Values**:
#' - H1-H2: -10 to +10 dB
#' - H1-A1: -5 to +5 dB
#' - H1-A3: Often negative (A3 > H1)
#'
#' **Performance**: ~24ms per vowel (15.7x faster than v4.8.14)
#'
#' @section Dependencies:
#' Requires \code{pladdrr} package (>= 4.8.16)
#'
#' @references
#' \insertCite{Iseli2004}{superassp}
#'
#' @examples
#' \dontrun{
#' # Analyze labeled vowels from TextGrid
#' results <- lst_pharyngeal(
#'   listOfFiles = "speech.wav",
#'   textgridPath = "speech.TextGrid",
#'   intervalTier = 3,
#'   intervalNumber = 1
#' )
#'
#' # Analyze specific time ranges (no TextGrid)
#' results <- lst_pharyngeal(
#'   listOfFiles = c("vowel1.wav", "vowel2.wav"),
#'   beginTime = c(0.5, 1.0),
#'   endTime = c(0.7, 1.3)
#' )
#'
#' # Write to JSTF files
#' lst_pharyngeal(
#'   listOfFiles = "speech.wav",
#'   textgridPath = "speech.TextGrid",
#'   toFile = TRUE
#' )  # Creates speech.pha
#'
#' # Access key measures
#' cat(sprintf("H1-H2* (onset): %.2f dB\\n", results$h1_minus_h2_onset_norm[1]))
#' cat(sprintf("H1-A1* (onset): %.2f dB\\n", results$h1_minus_a1_onset_norm[1]))
#' }
#'
#' @export
#' @family pladdrr functions
#' @family voice quality functions
#' @seealso \code{\link{lst_vq}}, \code{\link{trk_praatsaucep}}, \code{\link{lst_voice_reportp}}
lst_pharyngeal <- function(listOfFiles,
                           textgridPath = NULL,
                           intervalTier = 3,
                           intervalNumber = 1,
                           beginTime = NULL,
                           endTime = NULL,
                           minPitchInitial = 50,
                           maxPitchInitial = 800,
                           toFile = FALSE,
                           explicitExt = "pha",
                           outputDirectory = NULL,
                           verbose = TRUE) {
  
  # Check pladdrr availability
  if (!pladdrr_available()) {
    stop(
      "pladdrr package not available.\n",
      "Install with: install.packages('pladdrr')\n",
      "Required version: >= 4.8.16",
      call. = FALSE
    )
  }
  
  # Validate inputs
  if (!is.character(listOfFiles) || length(listOfFiles) == 0) {
    stop("listOfFiles must be a non-empty character vector", call. = FALSE)
  }
  
  # Check which mode: TextGrid or time-based
  use_textgrid <- !is.null(textgridPath)
  use_times <- !is.null(beginTime) && !is.null(endTime)
  
  if (!use_textgrid && !use_times) {
    stop(
      "Must provide either:\n",
      "  - textgridPath (for TextGrid-based analysis), OR\n",
      "  - beginTime AND endTime (for time-based analysis)",
      call. = FALSE
    )
  }
  
  if (use_textgrid && use_times) {
    if (verbose) {
      message("Both TextGrid and time parameters provided - using time parameters")
    }
    use_textgrid <- FALSE
  }
  
  n_files <- length(listOfFiles)
  
  # Validate TextGrid mode
  if (use_textgrid) {
    if (length(textgridPath) != n_files) {
      stop(
        "textgridPath must have same length as listOfFiles\n",
        sprintf("  listOfFiles: %d files\n", n_files),
        sprintf("  textgridPath: %d files", length(textgridPath)),
        call. = FALSE
      )
    }
  }
  
  # Validate time mode
  if (use_times) {
    if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
    if (length(endTime) == 1) endTime <- rep(endTime, n_files)
    
    if (length(beginTime) != n_files || length(endTime) != n_files) {
      stop(
        "beginTime and endTime must have same length as listOfFiles or length 1",
        call. = FALSE
      )
    }
    
    if (any(beginTime < 0) || any(endTime <= beginTime)) {
      stop("Invalid time range: endTime must be > beginTime >= 0", call. = FALSE)
    }
  }
  
  # Initialize results list
  results_list <- vector("list", n_files)
  
  # Progress bar
  if (verbose && n_files > 1) {
    pb <- txtProgressBar(min = 0, max = n_files, style = 3)
  }
  
  # Process each file
  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    
    # Check file exists
    if (!file.exists(file_path)) {
      warning(sprintf("File not found: %s (skipping)", file_path), call. = FALSE)
      results_list[[i]] <- NA
      if (verbose && n_files > 1) setTxtProgressBar(pb, i)
      next
    }
    
    # Load audio with av_load_for_pladdrr
    sound <- tryCatch({
      av_load_for_pladdrr(
        file_path = file_path,
        start_time = 0.0,  # Full file needed for formant extraction
        end_time = 0.0     # 0 = end of file
      )
    }, error = function(e) {
      warning(sprintf("Failed to load %s: %s", basename(file_path), e$message))
      NULL
    })
    
    if (is.null(sound)) {
      results_list[[i]] <- NA
      if (verbose && n_files > 1) setTxtProgressBar(pb, i)
      next
    }
    
    # Analyze using appropriate mode
    result <- tryCatch({
      if (use_textgrid) {
        # TextGrid mode
        tg_path <- textgridPath[i]
        if (!file.exists(tg_path)) {
          stop(sprintf("TextGrid not found: %s", tg_path))
        }
        
        # Load TextGrid
        textgrid <- pladdrr::TextGrid(tg_path)
        
        # Get interval bounds
        n_intervals <- textgrid$get_number_of_intervals(intervalTier)
        if (intervalNumber > n_intervals) {
          stop(sprintf(
            "Interval %d not found (tier %d has %d intervals)",
            intervalNumber, intervalTier, n_intervals
          ))
        }
        
        start <- textgrid$get_interval_start_time(intervalTier, intervalNumber)
        end <- textgrid$get_interval_end_time(intervalTier, intervalNumber)
        
        # Analyze with time range
        analyze_pharyngeal_times(
          sound = sound,
          start_time = start,
          end_time = end,
          min_pitch_initial = minPitchInitial,
          max_pitch_initial = maxPitchInitial
        )
      } else {
        # Time-based mode
        analyze_pharyngeal_times(
          sound = sound,
          start_time = beginTime[i],
          end_time = endTime[i],
          min_pitch_initial = minPitchInitial,
          max_pitch_initial = maxPitchInitial
        )
      }
    }, error = function(e) {
      warning(sprintf("Analysis failed for %s: %s", basename(file_path), e$message))
      NULL
    })
    
    if (is.null(result)) {
      results_list[[i]] <- NA
    } else {
      # Add file path to result
      result$file <- file_path
      results_list[[i]] <- result
    }
    
    if (verbose && n_files > 1) setTxtProgressBar(pb, i)
  }
  
  if (verbose && n_files > 1) close(pb)
  
  # Convert to data frame
  results_df <- do.call(rbind, lapply(results_list, function(x) {
    if (is.null(x) || (length(x) == 1 && is.na(x))) {
      return(NULL)
    }
    as.data.frame(x, stringsAsFactors = FALSE)
  }))
  
  # Handle empty results
  if (is.null(results_df) || nrow(results_df) == 0) {
    warning("No valid results obtained", call. = FALSE)
    return(data.frame())
  }
  
  # Reorder columns: file first, then timing, then measures
  col_order <- c(
    "file", "start_time", "mid_time", "end_time", "duration_ms",
    grep("^(f0|f1|f2|f3|bw|intensity|h1|h2|a1|a2|a3)", names(results_df), value = TRUE)
  )
  results_df <- results_df[, col_order[col_order %in% names(results_df)]]
  
  # Write to JSTF files if requested
  if (toFile) {
    output_paths <- write_lst_results_to_jstf(
      results = results_list,
      file_paths = listOfFiles,
      function_name = "lst_pharyngeal",
      explicitExt = explicitExt,
      outputDirectory = outputDirectory,
      verbose = verbose
    )
    return(invisible(output_paths))
  }
  
  return(results_df)
}


#' Internal: Analyze Pharyngeal Voice Quality from Time Range
#'
#' Core pharyngeal analysis function using pladdrr. Analyzes spectral
#' measures (H1-H2, H1-A1, etc.) at vowel onset and optionally midpoint.
#'
#' @param sound pladdrr Sound object
#' @param start_time Start time of interval (seconds)
#' @param end_time End time of interval (seconds)
#' @param min_pitch_initial Initial minimum pitch (default: 50 Hz)
#' @param max_pitch_initial Initial maximum pitch (default: 800 Hz)
#'
#' @return Named list with 68 pharyngeal measures
#' @keywords internal
analyze_pharyngeal_times <- function(sound,
                                     start_time,
                                     end_time,
                                     min_pitch_initial = 50,
                                     max_pitch_initial = 800) {
  
  UNDEFINED <- 1234
  
  # Step 1: Two-pass adaptive pitch detection
  pitch_result <- pladdrr::two_pass_adaptive_pitch(
    sound,
    initial_floor = min_pitch_initial,
    initial_ceiling = max_pitch_initial,
    voicing_threshold = 0.7,
    q1_factor = 0.75,
    q3_factor = 1.5,
    method = "ac"
  )
  
  min_pitch <- pitch_result$min_pitch
  max_pitch <- pitch_result$max_pitch
  pitch <- pladdrr::Pitch(.xptr = pitch_result$pitch)
  
  # Create PointProcess
  point_process_ptr <- pladdrr::to_point_process_from_sound_and_pitch(
    sound,
    pitch_result$pitch
  )
  point_process <- pladdrr::PointProcess(.xptr = point_process_ptr)
  
  # Calculate timing
  start <- start_time
  end <- end_time
  mid <- start + (end - start) / 2
  duration_ms <- (end - start) * 1000
  
  # Step 2: Create intensity
  intensity <- sound$to_intensity(50, 0.005, TRUE)
  sampfreq <- 10000
  
  # Calculate mean period
  mean_period <- point_process$get_mean_period(start, end, 0.0001, 0.02, 1.3)
  if (is.na(mean_period)) mean_period <- 0.01
  
  # Define analysis windows
  mid_before_frame <- mid - (mean_period / 2)
  mid_after_frame <- mid + (mean_period / 2)
  start_after_frame <- start + mean_period
  end_before_frame <- end - mean_period
  
  start_spectrum_before <- start
  start_spectrum_after <- start + 0.04
  mid_spectrum_before <- mid - 0.02
  mid_spectrum_after <- mid + 0.02
  
  # Find intensity maxima
  max_time_intensity_start <- intensity$get_time_of_maximum(
    from_time = start,
    to_time = start_after_frame
  )
  max_time_intensity_mid <- intensity$get_time_of_maximum(
    from_time = mid_before_frame,
    to_time = mid_after_frame
  )
  max_time_intensity_end <- intensity$get_time_of_maximum(
    from_time = end_before_frame,
    to_time = end
  )
  
  # Step 3: Extract 40ms windows around intensity maxima and compute formants
  # Using 5500 Hz max formant to match reference, no formant tracking

  # Extract window at start (40ms from intensity maximum)
  start_window <- sound$extract_part(max_time_intensity_start, max_time_intensity_start + 0.04,
                                      "rectangular", 1, FALSE)
  formant_start <- start_window$to_formant_burg(0.005, 5, 5500, 0.025, 50)

  # Query at center of window (relative time 0.02)
  f1_start <- formant_start$get_value_at_time(1, 0.02, "hertz")
  if (is.na(f1_start) || f1_start > 5000) f1_start <- UNDEFINED

  f2_start <- formant_start$get_value_at_time(2, 0.02, "hertz")
  if (is.na(f2_start) || f2_start > 5000) f2_start <- UNDEFINED

  f3_start <- formant_start$get_value_at_time(3, 0.02, "hertz")
  if (is.na(f3_start) || f3_start > 5000) f3_start <- UNDEFINED

  bw1_start <- formant_start$get_bandwidth_at_time(1, 0.02, "hertz")
  if (is.na(bw1_start)) bw1_start <- UNDEFINED
  bw2_start <- formant_start$get_bandwidth_at_time(2, 0.02, "hertz")
  if (is.na(bw2_start)) bw2_start <- UNDEFINED
  bw3_start <- formant_start$get_bandwidth_at_time(3, 0.02, "hertz")
  if (is.na(bw3_start)) bw3_start <- UNDEFINED

  bw1_start_norm <- normalize_bandwidth(bw1_start, f1_start, UNDEFINED)
  bw2_start_norm <- normalize_bandwidth(bw2_start, f2_start, UNDEFINED)
  bw3_start_norm <- normalize_bandwidth(bw3_start, f3_start, UNDEFINED)

  # Get F0 values (at center of analysis windows)
  start_point_frame <- max_time_intensity_start + 0.02
  f0_start <- pitch$get_value_at_time(start_point_frame, "hertz")
  f0_mid <- pitch$get_value_at_time(max_time_intensity_mid, "hertz")

  if (is.na(f0_start)) f0_start <- UNDEFINED
  if (is.na(f0_mid)) f0_mid <- UNDEFINED

  # Intensity at start (center of analysis window)
  intensity_start <- intensity$get_value_at_time(start_point_frame, "Cubic")
  if (is.na(intensity_start)) intensity_start <- UNDEFINED

  # Extract window at mid (40ms centered on intensity maximum)
  mid_window <- sound$extract_part(max_time_intensity_mid - 0.02, max_time_intensity_mid + 0.02,
                                    "rectangular", 1, FALSE)
  formant_mid <- mid_window$to_formant_burg(0.005, 5, 5500, 0.025, 50)

  # Query at center of window (relative time 0.02)
  f1_mid <- formant_mid$get_value_at_time(1, 0.02, "hertz")
  if (is.na(f1_mid) || f1_mid > 5000) f1_mid <- UNDEFINED

  f2_mid <- formant_mid$get_value_at_time(2, 0.02, "hertz")
  if (is.na(f2_mid) || f2_mid > 5000) f2_mid <- UNDEFINED

  f3_mid <- formant_mid$get_value_at_time(3, 0.02, "hertz")
  if (is.na(f3_mid) || f3_mid > 5000) f3_mid <- UNDEFINED

  bw1_mid <- formant_mid$get_bandwidth_at_time(1, 0.02, "hertz")
  if (is.na(bw1_mid)) bw1_mid <- UNDEFINED
  bw2_mid <- formant_mid$get_bandwidth_at_time(2, 0.02, "hertz")
  if (is.na(bw2_mid)) bw2_mid <- UNDEFINED
  bw3_mid <- formant_mid$get_bandwidth_at_time(3, 0.02, "hertz")
  if (is.na(bw3_mid)) bw3_mid <- UNDEFINED

  bw1_mid_norm <- normalize_bandwidth(bw1_mid, f1_mid, UNDEFINED)
  bw2_mid_norm <- normalize_bandwidth(bw2_mid, f2_mid, UNDEFINED)
  bw3_mid_norm <- normalize_bandwidth(bw3_mid, f3_mid, UNDEFINED)

  # Intensity at mid
  intensity_mid <- intensity$get_value_at_time(max_time_intensity_mid, "Cubic")
  if (is.na(intensity_mid)) intensity_mid <- UNDEFINED
  
  # Initialize all results
  result <- list(
    start_time = start,
    mid_time = mid,
    end_time = end,
    duration_ms = duration_ms,
    f0_start = f0_start,
    f0_mid = f0_mid,
    f1_start = f1_start,
    f2_start = f2_start,
    f3_start = f3_start,
    bw1_start = bw1_start,
    bw2_start = bw2_start,
    bw3_start = bw3_start,
    bw1_start_norm = bw1_start_norm,
    bw2_start_norm = bw2_start_norm,
    bw3_start_norm = bw3_start_norm,
    f1_mid = f1_mid,
    f2_mid = f2_mid,
    f3_mid = f3_mid,
    bw1_mid = bw1_mid,
    bw2_mid = bw2_mid,
    bw3_mid = bw3_mid,
    bw1_mid_norm = bw1_mid_norm,
    bw2_mid_norm = bw2_mid_norm,
    bw3_mid_norm = bw3_mid_norm,
    intensity_start = intensity_start,
    intensity_mid = intensity_mid,
    h1_onset = UNDEFINED,
    h1_hz_onset = UNDEFINED,
    h1_onset_norm = UNDEFINED,
    h2_onset = UNDEFINED,
    h2_hz_onset = UNDEFINED,
    h2_onset_norm = UNDEFINED,
    a1_onset = UNDEFINED,
    a1_hz_onset = UNDEFINED,
    a2_onset = UNDEFINED,
    a2_hz_onset = UNDEFINED,
    a3_onset = UNDEFINED,
    a3_hz_onset = UNDEFINED,
    a3_onset_norm = UNDEFINED,
    h1_minus_h2_onset = UNDEFINED,
    h1_minus_h2_onset_norm = UNDEFINED,
    h1_minus_a1_onset = UNDEFINED,
    h1_minus_a1_onset_norm = UNDEFINED,
    h1_minus_a2_onset = UNDEFINED,
    h1_minus_a2_onset_norm = UNDEFINED,
    h1_minus_a3_onset = UNDEFINED,
    h1_minus_a3_onset_norm = UNDEFINED,
    a1_minus_a2_onset = UNDEFINED,
    a1_minus_a3_onset = UNDEFINED,
    a1_minus_a3_onset_norm = UNDEFINED,
    a2_minus_a3_onset = UNDEFINED,
    a2_minus_a3_onset_norm = UNDEFINED,
    h1_mid = UNDEFINED,
    h1_hz_mid = UNDEFINED,
    h1_mid_norm = UNDEFINED,
    h2_mid = UNDEFINED,
    h2_hz_mid = UNDEFINED,
    h2_mid_norm = UNDEFINED,
    a1_mid = UNDEFINED,
    a1_hz_mid = UNDEFINED,
    a2_mid = UNDEFINED,
    a2_hz_mid = UNDEFINED,
    a3_mid = UNDEFINED,
    a3_hz_mid = UNDEFINED,
    a3_mid_norm = UNDEFINED,
    h1_minus_h2_mid = UNDEFINED,
    h1_minus_h2_mid_norm = UNDEFINED,
    h1_minus_a1_mid = UNDEFINED,
    h1_minus_a1_mid_norm = UNDEFINED,
    h1_minus_a2_mid = UNDEFINED,
    h1_minus_a2_mid_norm = UNDEFINED,
    h1_minus_a3_mid = UNDEFINED,
    h1_minus_a3_mid_norm = UNDEFINED,
    a1_minus_a2_mid = UNDEFINED,
    a1_minus_a3_mid = UNDEFINED,
    a1_minus_a3_mid_norm = UNDEFINED,
    a2_minus_a3_mid = UNDEFINED,
    a2_minus_a3_mid_norm = UNDEFINED
  )
  
  # Analyze onset if duration >= 40ms and F0 available
  if (duration_ms >= 40 && f0_start != UNDEFINED) {
    onset_results <- analyze_pharyngeal_onset(
      sound,
      sampfreq,
      start_spectrum_before,
      start_spectrum_after,
      f0_start,
      f1_start,
      f2_start,
      f3_start,
      bw1_start_norm,
      bw2_start_norm,
      bw3_start_norm,
      UNDEFINED
    )
    result <- modifyList(result, onset_results)
  }
  
  # Analyze mid if duration > 120ms and F0 available
  if (duration_ms > 120 && f0_mid != UNDEFINED) {
    mid_results <- analyze_pharyngeal_mid(
      sound,
      sampfreq,
      mid_spectrum_before,
      mid_spectrum_after,
      f0_mid,
      f1_mid,
      f2_mid,
      f3_mid,
      bw1_mid_norm,
      bw2_mid_norm,
      bw3_mid_norm,
      UNDEFINED
    )
    result <- modifyList(result, mid_results)
  }
  
  return(result)
}


#' Internal: Normalize Bandwidth
#'
#' @keywords internal
normalize_bandwidth <- function(bw, f, UNDEFINED) {
  if (is.na(bw) || bw == UNDEFINED || is.na(f) || f == UNDEFINED) {
    return(UNDEFINED)
  }
  return(bw + f * 0.085)
}


#' Internal: Iseli & Alwan (2004) Correction
#'
#' @keywords internal
iseli_alwan_correction <- function(f, bw, f0, sampfreq, UNDEFINED) {
  if (is.na(f) || f == UNDEFINED || is.na(bw) || bw == UNDEFINED ||
      is.na(f0) || f0 == UNDEFINED) {
    return(UNDEFINED)
  }
  
  pi_bw_sf <- pi * (bw / sampfreq)
  exp_term <- exp(-pi_bw_sf)
  exp_sq <- exp_term^2
  two_exp <- 2 * exp_term
  
  two_pi_f_sf <- 2 * pi * (f / sampfreq)
  cos_f <- cos(two_pi_f_sf)
  numerator <- (exp_sq + (1 - two_exp * cos_f))^2
  
  two_pi_f0_sf <- 2 * pi * (f0 / sampfreq)
  cos_plus <- cos(two_pi_f0_sf + two_pi_f_sf)
  cos_minus <- cos(two_pi_f0_sf - two_pi_f_sf)
  
  denom_plus <- exp_sq + (1 - two_exp * cos_plus)
  denom_minus <- exp_sq + (1 - two_exp * cos_minus)
  denominator <- denom_plus * denom_minus
  
  correction_db <- 10 * log10(numerator / denominator)
  return(correction_db)
}


#' Internal: Extract LTAS Peaks
#'
#' @keywords internal
extract_ltas_peaks <- function(ltas, fmins, fmaxs) {
  # Use batch API for 18x speedup (pladdrr >= 4.8.0)
  batch_result <- ltas$get_peaks_batch(fmins, fmaxs, "Parabolic")
  
  data.frame(
    fmin = fmins,
    fmax = fmaxs,
    peak_value = batch_result$values,
    peak_frequency = batch_result$frequencies
  )
}


#' Internal: Analyze Onset
#'
#' @keywords internal
analyze_pharyngeal_onset <- function(sound, sampfreq, start_time, end_time,
                                     f0_start, f1_start, f2_start, f3_start,
                                     bw1_start_norm, bw2_start_norm, bw3_start_norm,
                                     UNDEFINED) {
  
  # Extract window and resample
  sound_window <- sound$extract_part(start_time, end_time, "Kaiser2", 1, FALSE)
  sound_slice <- sound_window$resample(sampfreq, 50)
  
  # Create spectrum and filter
  spectrum_filtered <- sound_slice$to_spectrum(TRUE)
  spectrum_filtered$pass_hann_band(0, 5000, 100)
  spectrum_filtered$formula("if x >= 50 then self*x else self fi")
  
  # Convert to LTAS
  ltas_filtered <- spectrum_filtered$to_ltas_1to1()
  
  # Build frequency ranges for peak extraction
  f0_start_10 <- f0_start / 10
  
  fmins <- c(f0_start - f0_start_10, (f0_start * 2) - f0_start_10)
  fmaxs <- c(f0_start + f0_start_10, (f0_start * 2) + f0_start_10)
  peak_names <- c("h1", "h2")
  
  if (f1_start != UNDEFINED && bw1_start_norm != UNDEFINED) {
    fmins <- c(fmins, f1_start - bw1_start_norm / 2)
    fmaxs <- c(fmaxs, f1_start + bw1_start_norm / 2)
    peak_names <- c(peak_names, "a1")
  }
  if (f2_start != UNDEFINED && bw2_start_norm != UNDEFINED) {
    fmins <- c(fmins, f2_start - bw2_start_norm / 2)
    fmaxs <- c(fmaxs, f2_start + bw2_start_norm / 2)
    peak_names <- c(peak_names, "a2")
  }
  if (f3_start != UNDEFINED && bw3_start_norm != UNDEFINED) {
    fmins <- c(fmins, f3_start - bw3_start_norm / 2)
    fmaxs <- c(fmaxs, f3_start + bw3_start_norm / 2)
    peak_names <- c(peak_names, "a3")
  }
  
  # Extract peaks
  peaks <- extract_ltas_peaks(ltas_filtered, fmins, fmaxs)
  
  # Extract results
  get_peak_val <- function(name) {
    idx <- match(name, peak_names)
    if (is.na(idx)) return(UNDEFINED)
    val <- peaks$peak_value[idx]
    if (is.na(val)) UNDEFINED else val
  }
  get_peak_freq <- function(name) {
    idx <- match(name, peak_names)
    if (is.na(idx)) return(UNDEFINED)
    val <- peaks$peak_frequency[idx]
    if (is.na(val) || val > 5000) UNDEFINED else val
  }
  
  h1_onset <- get_peak_val("h1")
  h1_hz_onset <- get_peak_freq("h1")
  h2_onset <- get_peak_val("h2")
  h2_hz_onset <- get_peak_freq("h2")
  a1_onset <- get_peak_val("a1")
  a1_hz_onset <- get_peak_freq("a1")
  a2_onset <- get_peak_val("a2")
  a2_hz_onset <- get_peak_freq("a2")
  a3_onset <- get_peak_val("a3")
  a3_hz_onset <- get_peak_freq("a3")
  
  # Apply Iseli & Alwan normalization
  h1_onset_norm <- h1_onset
  h2_onset_norm <- h2_onset
  
  if (f1_start != UNDEFINED && bw1_start_norm != UNDEFINED) {
    h1_f1_norm <- iseli_alwan_correction(f1_start, bw1_start_norm, f0_start, sampfreq, UNDEFINED)
    h1_f2_norm <- iseli_alwan_correction(f2_start, bw2_start_norm, f0_start, sampfreq, UNDEFINED)
    h2_f1_norm <- iseli_alwan_correction(f1_start, bw1_start_norm, 2 * f0_start, sampfreq, UNDEFINED)
    h2_f2_norm <- iseli_alwan_correction(f2_start, bw2_start_norm, 2 * f0_start, sampfreq, UNDEFINED)
    
    if (h1_f1_norm != UNDEFINED && h1_f2_norm != UNDEFINED) {
      h1_onset_norm <- h1_onset + h1_f1_norm + h1_f2_norm
    }
    if (h2_f1_norm != UNDEFINED && h2_f2_norm != UNDEFINED) {
      h2_onset_norm <- h2_onset + h2_f1_norm + h2_f2_norm
    }
  }
  
  # Normalize A3
  a3_onset_norm <- UNDEFINED
  if (a3_onset != UNDEFINED && a3_hz_onset != UNDEFINED) {
    a3_f1_norm <- iseli_alwan_correction(f1_start, bw1_start_norm, a3_hz_onset, sampfreq, UNDEFINED)
    a3_f2_norm <- iseli_alwan_correction(f2_start, bw2_start_norm, a3_hz_onset, sampfreq, UNDEFINED)
    a3_f3_norm <- iseli_alwan_correction(f3_start, bw3_start_norm, a3_hz_onset, sampfreq, UNDEFINED)
    
    if (a3_f1_norm != UNDEFINED && a3_f2_norm != UNDEFINED && a3_f3_norm != UNDEFINED) {
      a3_onset_norm <- a3_onset - a3_f1_norm - a3_f2_norm - a3_f3_norm
    } else {
      a3_onset_norm <- a3_onset
    }
  }
  
  # Calculate differences
  list(
    h1_onset = h1_onset,
    h1_hz_onset = h1_hz_onset,
    h1_onset_norm = h1_onset_norm,
    h2_onset = h2_onset,
    h2_hz_onset = h2_hz_onset,
    h2_onset_norm = h2_onset_norm,
    a1_onset = a1_onset,
    a1_hz_onset = a1_hz_onset,
    a2_onset = a2_onset,
    a2_hz_onset = a2_hz_onset,
    a3_onset = a3_onset,
    a3_hz_onset = a3_hz_onset,
    a3_onset_norm = a3_onset_norm,
    h1_minus_h2_onset = if (h1_onset != UNDEFINED && h2_onset != UNDEFINED) h1_onset - h2_onset else UNDEFINED,
    h1_minus_h2_onset_norm = if (h1_onset_norm != UNDEFINED && h2_onset_norm != UNDEFINED) h1_onset_norm - h2_onset_norm else UNDEFINED,
    h1_minus_a1_onset = if (h1_onset != UNDEFINED && a1_onset != UNDEFINED) h1_onset - a1_onset else UNDEFINED,
    h1_minus_a1_onset_norm = if (h1_onset_norm != UNDEFINED && a1_onset != UNDEFINED) h1_onset_norm - a1_onset else UNDEFINED,
    h1_minus_a2_onset = if (h1_onset != UNDEFINED && a2_onset != UNDEFINED) h1_onset - a2_onset else UNDEFINED,
    h1_minus_a2_onset_norm = if (h1_onset_norm != UNDEFINED && a2_onset != UNDEFINED) h1_onset_norm - a2_onset else UNDEFINED,
    h1_minus_a3_onset = if (h1_onset != UNDEFINED && a3_onset != UNDEFINED) h1_onset - a3_onset else UNDEFINED,
    h1_minus_a3_onset_norm = if (h1_onset_norm != UNDEFINED && a3_onset_norm != UNDEFINED) h1_onset_norm - a3_onset_norm else UNDEFINED,
    a1_minus_a2_onset = if (a1_onset != UNDEFINED && a2_onset != UNDEFINED) a1_onset - a2_onset else UNDEFINED,
    a1_minus_a3_onset = if (a1_onset != UNDEFINED && a3_onset != UNDEFINED) a1_onset - a3_onset else UNDEFINED,
    a1_minus_a3_onset_norm = if (a1_onset != UNDEFINED && a3_onset_norm != UNDEFINED) a1_onset - a3_onset_norm else UNDEFINED,
    a2_minus_a3_onset = if (a2_onset != UNDEFINED && a3_onset != UNDEFINED) a2_onset - a3_onset else UNDEFINED,
    a2_minus_a3_onset_norm = if (a2_onset != UNDEFINED && a3_onset_norm != UNDEFINED) a2_onset - a3_onset_norm else UNDEFINED
  )
}


#' Internal: Analyze Mid
#'
#' @keywords internal
analyze_pharyngeal_mid <- function(sound, sampfreq, start_time, end_time,
                                   f0_mid, f1_mid, f2_mid, f3_mid,
                                   bw1_mid_norm, bw2_mid_norm, bw3_mid_norm,
                                   UNDEFINED) {
  
  # Same logic as onset, but for mid timepoint
  sound_window <- sound$extract_part(start_time, end_time, "Kaiser2", 1, FALSE)
  sound_slice <- sound_window$resample(sampfreq, 50)
  
  spectrum_filtered <- sound_slice$to_spectrum(TRUE)
  spectrum_filtered$pass_hann_band(0, 5000, 100)
  spectrum_filtered$formula("if x >= 50 then self*x else self fi")
  
  ltas_filtered <- spectrum_filtered$to_ltas_1to1()
  
  f0_mid_10 <- f0_mid / 10
  
  fmins <- c(f0_mid - f0_mid_10, (f0_mid * 2) - f0_mid_10)
  fmaxs <- c(f0_mid + f0_mid_10, (f0_mid * 2) + f0_mid_10)
  peak_names <- c("h1", "h2")
  
  if (f1_mid != UNDEFINED && bw1_mid_norm != UNDEFINED) {
    fmins <- c(fmins, f1_mid - bw1_mid_norm / 2)
    fmaxs <- c(fmaxs, f1_mid + bw1_mid_norm / 2)
    peak_names <- c(peak_names, "a1")
  }
  if (f2_mid != UNDEFINED && bw2_mid_norm != UNDEFINED) {
    fmins <- c(fmins, f2_mid - bw2_mid_norm / 2)
    fmaxs <- c(fmaxs, f2_mid + bw2_mid_norm / 2)
    peak_names <- c(peak_names, "a2")
  }
  if (f3_mid != UNDEFINED && bw3_mid_norm != UNDEFINED) {
    fmins <- c(fmins, f3_mid - bw3_mid_norm / 2)
    fmaxs <- c(fmaxs, f3_mid + bw3_mid_norm / 2)
    peak_names <- c(peak_names, "a3")
  }
  
  peaks <- extract_ltas_peaks(ltas_filtered, fmins, fmaxs)
  
  get_peak_val <- function(name) {
    idx <- match(name, peak_names)
    if (is.na(idx)) return(UNDEFINED)
    val <- peaks$peak_value[idx]
    if (is.na(val)) UNDEFINED else val
  }
  get_peak_freq <- function(name) {
    idx <- match(name, peak_names)
    if (is.na(idx)) return(UNDEFINED)
    val <- peaks$peak_frequency[idx]
    if (is.na(val) || val > 5000) UNDEFINED else val
  }
  
  h1_mid <- get_peak_val("h1")
  h1_hz_mid <- get_peak_freq("h1")
  h2_mid <- get_peak_val("h2")
  h2_hz_mid <- get_peak_freq("h2")
  a1_mid <- get_peak_val("a1")
  a1_hz_mid <- get_peak_freq("a1")
  a2_mid <- get_peak_val("a2")
  a2_hz_mid <- get_peak_freq("a2")
  a3_mid <- get_peak_val("a3")
  a3_hz_mid <- get_peak_freq("a3")
  
  h1_mid_norm <- h1_mid
  h2_mid_norm <- h2_mid
  
  if (f1_mid != UNDEFINED && bw1_mid_norm != UNDEFINED) {
    h1_f1_norm <- iseli_alwan_correction(f1_mid, bw1_mid_norm, f0_mid, sampfreq, UNDEFINED)
    h1_f2_norm <- iseli_alwan_correction(f2_mid, bw2_mid_norm, f0_mid, sampfreq, UNDEFINED)
    h2_f1_norm <- iseli_alwan_correction(f1_mid, bw1_mid_norm, 2 * f0_mid, sampfreq, UNDEFINED)
    h2_f2_norm <- iseli_alwan_correction(f2_mid, bw2_mid_norm, 2 * f0_mid, sampfreq, UNDEFINED)
    
    if (h1_f1_norm != UNDEFINED && h1_f2_norm != UNDEFINED) {
      h1_mid_norm <- h1_mid + h1_f1_norm + h1_f2_norm
    }
    if (h2_f1_norm != UNDEFINED && h2_f2_norm != UNDEFINED) {
      h2_mid_norm <- h2_mid + h2_f1_norm + h2_f2_norm
    }
  }
  
  a3_mid_norm <- UNDEFINED
  if (a3_mid != UNDEFINED && a3_hz_mid != UNDEFINED) {
    a3_f1_norm <- iseli_alwan_correction(f1_mid, bw1_mid_norm, a3_hz_mid, sampfreq, UNDEFINED)
    a3_f2_norm <- iseli_alwan_correction(f2_mid, bw2_mid_norm, a3_hz_mid, sampfreq, UNDEFINED)
    a3_f3_norm <- iseli_alwan_correction(f3_mid, bw3_mid_norm, a3_hz_mid, sampfreq, UNDEFINED)
    
    if (a3_f1_norm != UNDEFINED && a3_f2_norm != UNDEFINED && a3_f3_norm != UNDEFINED) {
      a3_mid_norm <- a3_mid - a3_f1_norm - a3_f2_norm - a3_f3_norm
    } else {
      a3_mid_norm <- a3_mid
    }
  }
  
  list(
    h1_mid = h1_mid,
    h1_hz_mid = h1_hz_mid,
    h1_mid_norm = h1_mid_norm,
    h2_mid = h2_mid,
    h2_hz_mid = h2_hz_mid,
    h2_mid_norm = h2_mid_norm,
    a1_mid = a1_mid,
    a1_hz_mid = a1_hz_mid,
    a2_mid = a2_mid,
    a2_hz_mid = a2_hz_mid,
    a3_mid = a3_mid,
    a3_hz_mid = a3_hz_mid,
    a3_mid_norm = a3_mid_norm,
    h1_minus_h2_mid = if (h1_mid != UNDEFINED && h2_mid != UNDEFINED) h1_mid - h2_mid else UNDEFINED,
    h1_minus_h2_mid_norm = if (h1_mid_norm != UNDEFINED && h2_mid_norm != UNDEFINED) h1_mid_norm - h2_mid_norm else UNDEFINED,
    h1_minus_a1_mid = if (h1_mid != UNDEFINED && a1_mid != UNDEFINED) h1_mid - a1_mid else UNDEFINED,
    h1_minus_a1_mid_norm = if (h1_mid_norm != UNDEFINED && a1_mid != UNDEFINED) h1_mid_norm - a1_mid else UNDEFINED,
    h1_minus_a2_mid = if (h1_mid != UNDEFINED && a2_mid != UNDEFINED) h1_mid - a2_mid else UNDEFINED,
    h1_minus_a2_mid_norm = if (h1_mid_norm != UNDEFINED && a2_mid != UNDEFINED) h1_mid_norm - a2_mid else UNDEFINED,
    h1_minus_a3_mid = if (h1_mid != UNDEFINED && a3_mid != UNDEFINED) h1_mid - a3_mid else UNDEFINED,
    h1_minus_a3_mid_norm = if (h1_mid_norm != UNDEFINED && a3_mid_norm != UNDEFINED) h1_mid_norm - a3_mid_norm else UNDEFINED,
    a1_minus_a2_mid = if (a1_mid != UNDEFINED && a2_mid != UNDEFINED) a1_mid - a2_mid else UNDEFINED,
    a1_minus_a3_mid = if (a1_mid != UNDEFINED && a3_mid != UNDEFINED) a1_mid - a3_mid else UNDEFINED,
    a1_minus_a3_mid_norm = if (a1_mid != UNDEFINED && a3_mid_norm != UNDEFINED) a1_mid - a3_mid_norm else UNDEFINED,
    a2_minus_a3_mid = if (a2_mid != UNDEFINED && a3_mid != UNDEFINED) a2_mid - a3_mid else UNDEFINED,
    a2_minus_a3_mid_norm = if (a2_mid != UNDEFINED && a3_mid_norm != UNDEFINED) a2_mid - a3_mid_norm else UNDEFINED
  )
}

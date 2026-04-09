#' Voice Report Analysis (pladdrr)
#'
#' Compute comprehensive voice quality measures from audio files using pladdrr.
#' This function replicates Praat's "Voice report" functionality and returns
#' 30 voice quality measures including pitch, jitter, shimmer, and harmonicity.
#'
#' @param listOfFiles Character vector of audio file paths
#' @param beginTime Numeric. Start time for analysis (seconds). Default 0.0
#' @param endTime Numeric. End time for analysis (seconds). 0 = end of file. Default 0.0
#' @param selectionOffset Numeric. Offset from start for selection window (seconds). Default 0.0
#' @param selectionLength Numeric. Length of selection window (seconds). 0 = entire file. Default 0.0
#' @param minF Numeric. Minimum pitch for tracking (Hz). Default 75
#' @param maxF Numeric. Maximum pitch for tracking (Hz). Default 600
#' @param windowShape Character. Window type for extraction. Default "Gaussian1"
#' @param relativeWidth Numeric. Relative window width. Default 1.0
#' @param max_period_factor Numeric. Max period factor for jitter/shimmer. Default 1.3
#' @param max_ampl_factor Numeric. Max amplitude factor for shimmer. Default 1.6
#' @param silence_threshold Numeric. Silence threshold for pitch. Default 0.03
#' @param voicing_threshold Numeric. Voicing threshold for pitch. Default 0.45
#' @param octave_cost Numeric. Octave cost for pitch. Default 0.01
#' @param octave_jump_cost Numeric. Octave jump cost for pitch. Default 0.35
#' @param voiced_unvoiced_cost Numeric. Voiced/unvoiced cost for pitch. Default 0.14
#' @param toFile Logical. Write to JSTF file? Default FALSE
#' @param explicitExt Character. Output file extension. Default "pvr"
#' @param outputDirectory Character. Output directory. NULL = input directory. Default NULL
#' @param verbose Logical. Show progress? Default TRUE
#'
#' @return If toFile=FALSE, data.frame with 30 voice measures per file.
#'   If toFile=TRUE, invisibly returns output file path(s).
#'
#'   The data.frame contains:
#'   \describe{
#'     \item{file}{Input filename}
#'     \item{Timing (4)}{start_time, end_time, selection_start, selection_end}
#'     \item{Pitch (5)}{median_pitch, mean_pitch, sd_pitch, min_pitch, max_pitch}
#'     \item{Pulses (4)}{num_pulses, num_periods, mean_period, sd_period}
#'     \item{Voicing (3)}{fraction_unvoiced, num_voice_breaks, degree_voice_breaks}
#'     \item{Jitter (5)}{jitter_local_percent, jitter_local_abs, jitter_rap_percent, jitter_ppq5_percent, jitter_ddp_percent}
#'     \item{Shimmer (6)}{shimmer_local_percent, shimmer_local_db, shimmer_apq3_percent, shimmer_apq5_percent, shimmer_apq11_percent, shimmer_dda_percent}
#'     \item{Harmonicity (3)}{mean_autocorrelation, mean_nhr, mean_hnr}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Single file
#' test_file <- system.file("samples/sustained/a1.wav", package = "superassp")
#' result <- lst_voice_report(test_file, verbose = FALSE)
#' print(result)
#'
#' # Multiple files
#' files <- c("vowel1.wav", "vowel2.wav")
#' results <- lst_voice_report(files)
#'
#' # With time windowing
#' result <- lst_voice_report(
#'   test_file,
#'   beginTime = 1.0,
#'   endTime = 3.0
#' )
#'
#' # Write to JSTF file
#' lst_voice_report(test_file, toFile = TRUE)  # Creates a1.pvr
#' }
lst_voice_report <- function(listOfFiles,
                              beginTime = 0.0,
                              endTime = 0.0,
                              selectionOffset = 0.0,
                              selectionLength = 0.0,
                              minF = 75,
                              maxF = 600,
                              windowShape = "Gaussian1",
                              relativeWidth = 1.0,
                              max_period_factor = 1.3,
                              max_ampl_factor = 1.6,
                              silence_threshold = 0.03,
                              voicing_threshold = 0.45,
                              octave_cost = 0.01,
                              octave_jump_cost = 0.35,
                              voiced_unvoiced_cost = 0.14,
                              toFile = FALSE,
                              explicitExt = "pvr",
                              outputDirectory = NULL,
                              verbose = TRUE) {

  # Check pladdrr availability
  if (!pladdrr_available()) {
    stop("pladdrr not available. Install with: install_pladdrr()")
  }

  # Validate window shape
  valid_windows <- c("rectangular", "triangular", "parabolic", "Hanning", "Hamming",
                     "Gaussian1", "Gaussian2", "Gaussian3", "Gaussian4", "Gaussian5",
                     "Kaiser1", "Kaiser2")
  if (!windowShape %in% valid_windows) {
    stop("Invalid windowShape: ", windowShape, ". Must be one of: ",
         paste(valid_windows, collapse = ", "))
  }

  # Validate files
  missing_files <- listOfFiles[!file.exists(listOfFiles)]
  if (length(missing_files) > 0) {
    stop("Files not found: ", paste(missing_files, collapse = ", "))
  }

  n_files <- length(listOfFiles)

  # Progress bar for batch processing
  if (verbose && n_files > 1) {
    pb <- txtProgressBar(min = 0, max = n_files, style = 3)
    on.exit(close(pb), add = TRUE)
  }

  # Process each file
  results <- vector("list", n_files)

  for (i in seq_along(listOfFiles)) {
    file <- listOfFiles[i]

    # Load sound with pladdrr
    sound <- pladdrr::Sound(file)
    sound_duration <- sound$.cpp$duration

    # Calculate analysis window
    actual_start <- beginTime
    actual_end <- if (endTime > 0) endTime else sound_duration

    start_at <- beginTime + selectionOffset
    end_at <- if (endTime > 0) endTime else sound_duration

    if (selectionLength > 0) {
      sel_end <- start_at + selectionLength
      if (sel_end < end_at) {
        end_at <- sel_end
      }
    }

    actual_sel_start <- start_at
    actual_sel_end <- end_at

    # Extract analysis part
    sound_part <- sound$extract_part(
      start_at, end_at,
      windowShape, relativeWidth,
      preserve_times = FALSE
    )

    # Create Pitch object (cross-correlation method)
    pitch_ptr <- pladdrr::to_pitch_cc_direct(
      sound_part,
      time_step = 0.0,  # auto
      pitch_floor = minF,
      pitch_ceiling = maxF,
      max_candidates = 15,
      very_accurate = TRUE,
      silence_threshold = silence_threshold,
      voicing_threshold = voicing_threshold,
      octave_cost = octave_cost,
      octave_jump_cost = octave_jump_cost,
      voiced_unvoiced_cost = voiced_unvoiced_cost
    )
    pitch <- pladdrr::Pitch(.xptr = pitch_ptr)

    # Create PointProcess from Sound and Pitch
    pp_ptr <- pladdrr::to_point_process_from_sound_and_pitch(sound_part, pitch_ptr)
    point_process <- pladdrr::PointProcess(.xptr = pp_ptr)

    # Analysis range (entire extracted part)
    analysis_start <- 0.0
    analysis_end <- 0.0  # 0 means entire duration

    # ========== PITCH MEASURES (5) ==========
    median_pitch <- pitch$get_quantile(analysis_start, analysis_end, 0.5, "hertz")
    mean_pitch <- pitch$get_mean(analysis_start, analysis_end, "hertz")
    sd_pitch <- pitch$get_standard_deviation(analysis_start, analysis_end, "hertz")
    min_pitch_val <- pitch$get_minimum(analysis_start, analysis_end, "hertz", "parabolic")
    max_pitch_val <- pitch$get_maximum(analysis_start, analysis_end, "hertz", "parabolic")

    # Handle NA/NaN
    median_pitch <- if (is.na(median_pitch)) NaN else median_pitch
    mean_pitch <- if (is.na(mean_pitch)) NaN else mean_pitch
    sd_pitch <- if (is.na(sd_pitch)) NaN else sd_pitch
    min_pitch_val <- if (is.na(min_pitch_val)) NaN else min_pitch_val
    max_pitch_val <- if (is.na(max_pitch_val)) NaN else max_pitch_val

    # ========== PULSE/PERIOD MEASURES (4) ==========
    num_pulses <- point_process$get_number_of_points()
    num_periods <- point_process$get_number_of_periods(
      analysis_start, analysis_end, 0.0001, 0.02, max_period_factor
    )
    mean_period <- point_process$get_mean_period(
      analysis_start, analysis_end, 0.0001, 0.02, max_period_factor
    )
    sd_period <- point_process$get_stdev_period(
      analysis_start, analysis_end, 0.0001, 0.02, max_period_factor
    )

    # ========== VOICING MEASURES (3) ==========
    voiced_frames <- pitch$count_voiced_frames()
    n_frames <- pitch$get_number_of_frames()
    fraction_unvoiced <- if (n_frames > 0) {
      (n_frames - voiced_frames) / n_frames
    } else {
      NaN
    }

    # Voice breaks (gaps in voiced regions)
    num_voice_breaks <- 0
    degree_voice_breaks <- 0.0

    duration <- sound_part$.cpp$duration
    if (duration > 0 && n_frames > 0) {
      # Find voiced intervals
      voiced_intervals <- list()
      in_voiced <- FALSE
      voiced_start <- 0

      for (j in 1:n_frames) {
        time_j <- pitch$get_time_from_frame(j)
        val <- tryCatch(
          pitch$get_value_at_time(time_j, "hertz", "linear"),
          error = function(e) NaN
        )

        if (!is.na(val) && !is.nan(val) && val > 0) {
          if (!in_voiced) {
            in_voiced <- TRUE
            voiced_start <- time_j
          }
        } else {
          if (in_voiced) {
            in_voiced <- FALSE
            voiced_intervals[[length(voiced_intervals) + 1]] <- c(voiced_start, time_j)
          }
        }
      }
      if (in_voiced) {
        voiced_intervals[[length(voiced_intervals) + 1]] <- c(voiced_start, pitch$get_end_time())
      }

      # Count breaks between voiced intervals
      if (length(voiced_intervals) > 1) {
        num_voice_breaks <- length(voiced_intervals) - 1
        break_duration <- 0
        for (k in 2:length(voiced_intervals)) {
          break_duration <- break_duration + (voiced_intervals[[k]][1] - voiced_intervals[[k-1]][2])
        }
        degree_voice_breaks <- break_duration / duration
      }
    }

    # ========== JITTER/SHIMMER MEASURES (11) ==========
    jitter_shimmer <- pladdrr::get_jitter_shimmer_batch(
      point_process,
      sound_part,
      from_time = analysis_start,
      to_time = analysis_end,
      period_floor = 0.0001,
      period_ceiling = 0.02,
      max_period_factor = max_period_factor,
      max_amplitude_factor = max_ampl_factor
    )

    # Jitter (5) - convert decimal to percentage
    jitter_local <- jitter_shimmer$jitter_local * 100
    jitter_local_abs <- jitter_shimmer$jitter_local_abs
    jitter_rap <- jitter_shimmer$jitter_rap * 100
    jitter_ppq5 <- jitter_shimmer$jitter_ppq5 * 100
    jitter_ddp <- jitter_shimmer$jitter_ddp * 100

    # Shimmer (6) - convert decimal to percentage
    shimmer_local <- jitter_shimmer$shimmer_local * 100
    shimmer_local_db <- jitter_shimmer$shimmer_local_db
    shimmer_apq3 <- jitter_shimmer$shimmer_apq3 * 100
    shimmer_apq5 <- jitter_shimmer$shimmer_apq5 * 100
    shimmer_apq11 <- jitter_shimmer$shimmer_apq11 * 100
    shimmer_dda <- jitter_shimmer$shimmer_dda * 100

    # ========== HARMONICITY MEASURES (3) ==========
    harmonicity_ptr <- pladdrr::to_harmonicity_direct(
      sound_part, 0.01, minF, 0.1, 1.0
    )
    harmonicity <- pladdrr::Harmonicity(.xptr = harmonicity_ptr)

    mean_hnr <- harmonicity$get_mean(analysis_start, analysis_end)
    if (is.na(mean_hnr)) mean_hnr <- NaN

    # NHR is inverse of HNR in linear scale
    mean_nhr <- if (!is.nan(mean_hnr) && mean_hnr != 0) {
      1 / (10^(mean_hnr / 10))
    } else {
      NaN
    }

    # Mean autocorrelation from pitch strength
    mean_autocorr <- 0.0
    autocorr_sum <- 0
    autocorr_count <- 0
    for (m in 1:n_frames) {
      strength <- tryCatch(
        pitch$get_strength_at_frame(m),
        error = function(e) NaN
      )
      if (!is.na(strength) && !is.nan(strength)) {
        autocorr_sum <- autocorr_sum + strength
        autocorr_count <- autocorr_count + 1
      }
    }
    mean_autocorr <- if (autocorr_count > 0) autocorr_sum / autocorr_count else NaN

    # Build result for this file
    results[[i]] <- data.frame(
      file = basename(file),
      start_time = actual_start,
      end_time = actual_end,
      selection_start = actual_sel_start,
      selection_end = actual_sel_end,
      median_pitch = median_pitch,
      mean_pitch = mean_pitch,
      sd_pitch = sd_pitch,
      min_pitch = min_pitch_val,
      max_pitch = max_pitch_val,
      num_pulses = num_pulses,
      num_periods = num_periods,
      mean_period = mean_period,
      sd_period = sd_period,
      fraction_unvoiced = fraction_unvoiced,
      num_voice_breaks = num_voice_breaks,
      degree_voice_breaks = degree_voice_breaks,
      jitter_local_percent = jitter_local,
      jitter_local_abs = jitter_local_abs,
      jitter_rap_percent = jitter_rap,
      jitter_ppq5_percent = jitter_ppq5,
      jitter_ddp_percent = jitter_ddp,
      shimmer_local_percent = shimmer_local,
      shimmer_local_db = shimmer_local_db,
      shimmer_apq3_percent = shimmer_apq3,
      shimmer_apq5_percent = shimmer_apq5,
      shimmer_apq11_percent = shimmer_apq11,
      shimmer_dda_percent = shimmer_dda,
      mean_autocorrelation = mean_autocorr,
      mean_nhr = mean_nhr,
      mean_hnr = mean_hnr,
      stringsAsFactors = FALSE
    )

    if (verbose && n_files > 1) {
      setTxtProgressBar(pb, i)
    }
  }

  # Handle file output
  if (toFile) {
    # Convert each data.frame to named list (remove 'file' column)
    results_list <- lapply(results, function(df) {
      df$file <- NULL
      as.list(df[1, ])  # Convert single row to named list
    })

    output_paths <- write_lst_results_to_jstf(
      results = results_list,
      file_paths = listOfFiles,
      beginTime = beginTime,
      endTime = endTime,
      function_name = "lst_voice_report",
      parameters = list(
        selectionOffset = selectionOffset,
        selectionLength = selectionLength,
        minF = minF,
        maxF = maxF,
        windowShape = windowShape,
        relativeWidth = relativeWidth,
        max_period_factor = max_period_factor,
        max_ampl_factor = max_ampl_factor,
        silence_threshold = silence_threshold,
        voicing_threshold = voicing_threshold,
        octave_cost = octave_cost,
        octave_jump_cost = octave_jump_cost,
        voiced_unvoiced_cost = voiced_unvoiced_cost
      ),
      explicitExt = explicitExt,
      outputDirectory = outputDirectory
    )

    return(invisible(output_paths))
  }

  # Combine results into single data.frame for in-memory return
  result_df <- do.call(rbind, results)
  rownames(result_df) <- NULL

  return(result_df)
}

# Function attributes
attr(lst_voice_report, "ext") <- "pvr"
attr(lst_voice_report, "outputType") <- "JSTF"
attr(lst_voice_report, "format") <- "JSON"
attr(lst_voice_report, "tracks") <- c(
  "start_time", "end_time", "selection_start", "selection_end",
  "median_pitch", "mean_pitch", "sd_pitch", "min_pitch", "max_pitch",
  "num_pulses", "num_periods", "mean_period", "sd_period",
  "fraction_unvoiced", "num_voice_breaks", "degree_voice_breaks",
  "jitter_local_percent", "jitter_local_abs", "jitter_rap_percent",
  "jitter_ppq5_percent", "jitter_ddp_percent",
  "shimmer_local_percent", "shimmer_local_db", "shimmer_apq3_percent",
  "shimmer_apq5_percent", "shimmer_apq11_percent", "shimmer_dda_percent",
  "mean_autocorrelation", "mean_nhr", "mean_hnr"
)

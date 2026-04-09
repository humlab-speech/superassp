#' Dysphonia Severity Index (DSI) Analysis (pladdrr)
#'
#' Compute the Dysphonia Severity Index (DSI) using pladdrr. DSI is a multiparametric
#' voice quality assessment combining maximum phonation time, softest intensity,
#' highest fundamental frequency, and jitter.
#'
#' DSI Formula (Wuyts et al., 2000):
#' DSI = 1.127 + 0.164*MPT - 0.038*I-low + 0.0053*F0-high - 5.30*Jitter
#'
#' @param softDF Data frame with soft phonation samples. Columns: absolute_file_path, start, end
#' @param highpitchDF Data frame with high pitch samples. Columns: absolute_file_path, start, end
#' @param maxprolongedDF Data frame with maximally prolonged vowel samples. Columns: absolute_file_path, start, end
#' @param stableDF Data frame with stable vowel samples for jitter (optional, defaults to maxprolongedDF)
#' @param use.calibration Logical. Apply calibration to intensity measurements. Default FALSE
#' @param db.calibration Numeric. Calibration factor in dB. Default 10
#' @param speaker.name Character. Speaker name (optional, for identification)
#' @param speaker.ID Character. Speaker ID (optional, for file naming)
#' @param speaker.dob Character. Speaker date of birth (optional)
#' @param session.datetime Character. Session datetime (optional)
#' @param pdf.path Character. Path for PDF output (not implemented in pladdrr version)
#' @param simple.output Logical. Simplified output (not used in pladdrr version)
#' @param overwrite.pdfs Logical. Overwrite PDFs (not used in pladdrr version)
#' @param praat_path Character. Praat path (not used in pladdrr version)
#' @param toFile Logical. Write to JSTF file? Default FALSE
#' @param explicitExt Character. Output file extension. Default "dsi"
#' @param outputDirectory Character. Output directory. NULL = first file's directory. Default NULL
#'
#' @return If toFile=FALSE, list with 6 elements:
#'   \describe{
#'     \item{ID}{Speaker ID}
#'     \item{Maximum_phonation_time}{MPT in seconds}
#'     \item{Softest_intensity_of_voiced_speech}{Minimum intensity in dB}
#'     \item{Maximum_fundamental_frequency}{Maximum F0 in Hz}
#'     \item{Jitter_ppq5}{5-point period perturbation quotient (%)}
#'     \item{Dysphonia_Severity_Index}{DSI composite score}
#'   }
#'   If toFile=TRUE, invisibly returns output file path.
#'
#' @references
#' Wuyts, F. L., De Bodt, M. S., Molenberghs, G., Remacle, M., Heylen, L.,
#' Millet, B., ... & Heyning, P. H. V. D. (2000). The dysphonia severity index:
#' an objective measure of vocal quality based on a multiparameter approach.
#' Journal of Speech, Language, and Hearing Research, 43(3), 796-809.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example DSI analysis
#' soft <- data.frame(
#'   absolute_file_path = "soft1.wav",
#'   start = 0,
#'   end = 2
#' )
#'
#' high <- data.frame(
#'   absolute_file_path = "high1.wav",
#'   start = 0,
#'   end = 1.5
#' )
#'
#' prolonged <- data.frame(
#'   absolute_file_path = c("vowel1.wav", "vowel2.wav"),
#'   start = c(0, 0),
#'   end = c(5, 6)
#' )
#'
#' result <- lst_dsi(
#'   softDF = soft,
#'   highpitchDF = high,
#'   maxprolongedDF = prolonged,
#'   speaker.ID = "P001"
#' )
#'
#' print(result)
#' }
lst_dsi <- function(softDF,
                     highpitchDF,
                     maxprolongedDF,
                     stableDF = NULL,
                     use.calibration = FALSE,
                     db.calibration = 10,
                     speaker.name = NULL,
                     speaker.ID = NULL,
                     speaker.dob = NULL,
                     session.datetime = NULL,
                     pdf.path = NULL,
                     simple.output = FALSE,
                     overwrite.pdfs = FALSE,
                     praat_path = NULL,
                     toFile = FALSE,
                     explicitExt = "dsi",
                     outputDirectory = NULL) {

  # Check pladdrr availability
  if (!pladdrr_available()) {
    stop("pladdrr not available. Install with: install_pladdrr()")
  }

  # Validate required columns
  required_cols <- c("absolute_file_path", "start", "end")

  # Default stableDF to maxprolongedDF if not provided
  if (is.null(stableDF)) {
    stableDF <- maxprolongedDF
  }

  # Check all dataframes are non-empty
  if (nrow(softDF) < 1 || nrow(highpitchDF) < 1 ||
      nrow(maxprolongedDF) < 1 || nrow(stableDF) < 1) {
    stop("All dataframes must be non-empty")
  }

  # Validate columns
  for (df_name in c("softDF", "highpitchDF", "maxprolongedDF", "stableDF")) {
    df <- get(df_name)
    if (!all(required_cols %in% names(df))) {
      stop(df_name, " must contain columns: ", paste(required_cols, collapse = ", "))
    }
  }

  # Validate file paths
  all_files <- unique(c(
    softDF$absolute_file_path,
    highpitchDF$absolute_file_path,
    maxprolongedDF$absolute_file_path,
    stableDF$absolute_file_path
  ))
  missing_files <- all_files[!file.exists(all_files)]
  if (length(missing_files) > 0) {
    stop("Files not found: ", paste(missing_files, collapse = ", "))
  }

  # Calculate MPT (Maximum Phonation Time)
  mpt <- calculate_mpt_pladdrr(maxprolongedDF)

  # Calculate minimum intensity (from soft phonation samples)
  minimum_intensity <- calculate_minimum_intensity_pladdrr(
    softDF, use.calibration, db.calibration
  )

  # Calculate maximum F0 (from high pitch samples)
  maximum_f0 <- calculate_maximum_f0_pladdrr(highpitchDF)

  # Calculate jitter ppq5 (from stable vowel samples)
  jitter_ppq5 <- calculate_jitter_ppq5_pladdrr(stableDF)

  # Calculate DSI
  # Formula: DSI = 1.127 + 0.164*MPT - 0.038*I-low + 0.0053*F0-high - 5.30*Jitter
  dsi <- 1.127 +
    0.164 * mpt -
    0.038 * minimum_intensity +
    0.0053 * maximum_f0 -
    5.30 * jitter_ppq5

  # Build result
  result_list <- list(
    ID = if (!is.null(speaker.ID)) speaker.ID else NA_character_,
    Maximum_phonation_time = mpt,
    Softest_intensity_of_voiced_speech = minimum_intensity,
    Maximum_fundamental_frequency = maximum_f0,
    Jitter_ppq5 = jitter_ppq5,
    Dysphonia_Severity_Index = dsi
  )

  # Handle file output
  if (toFile) {
    # Determine output file path
    primary_file <- softDF$absolute_file_path[1]
    analysis_begin <- 0.0
    analysis_end <- 0.0  # Full duration

    # Get audio metadata
    audio_info <- av::av_media_info(primary_file)
    sample_rate <- audio_info$audio$sample_rate
    audio_duration <- audio_info$duration

    # Create JSTF object
    json_obj <- create_json_track_obj(
      results = result_list,
      function_name = "lst_dsi",
      file_path = primary_file,
      sample_rate = sample_rate,
      audio_duration = audio_duration,
      beginTime = analysis_begin,
      endTime = audio_duration,
      parameters = list(
        use_calibration = use.calibration,
        db_calibration = db.calibration
      )
    )

    # Determine output filename (use speaker.ID if provided)
    if (!is.null(speaker.ID) && nchar(speaker.ID) > 0) {
      base_name <- speaker.ID
    } else {
      base_name <- tools::file_path_sans_ext(basename(primary_file))
    }

    # Determine output directory
    out_dir <- if (is.null(outputDirectory)) dirname(primary_file) else outputDirectory
    output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))

    # Write JSTF file
    write_json_track(json_obj, output_path)

    return(invisible(output_path))
  }

  return(result_list)
}

# Function attributes
attr(lst_dsi, "ext") <- "dsi"
attr(lst_dsi, "outputType") <- "JSTF"
attr(lst_dsi, "format") <- "JSON"
attr(lst_dsi, "tracks") <- c("ID", "Maximum_phonation_time",
                               "Softest_intensity_of_voiced_speech",
                               "Maximum_fundamental_frequency",
                               "Jitter_ppq5",
                               "Dysphonia_Severity_Index")


# ==================== HELPER FUNCTIONS ====================

#' Calculate Maximum Phonation Time
#' @keywords internal
calculate_mpt_pladdrr <- function(df) {
  files <- df$absolute_file_path

  # Use ultra-fast batch duration reading if all are file paths
  if (all(sapply(files, is.character))) {
    if (exists("get_durations_batch", where = asNamespace("pladdrr"))) {
      durations <- pladdrr::get_durations_batch(files)
      return(max(durations, na.rm = TRUE))
    }
  }

  # Fallback: Load each file and get duration
  durations <- numeric(length(files))
  for (i in seq_along(files)) {
    sound <- pladdrr::Sound(files[i])
    durations[i] <- sound$.cpp$duration
  }

  return(max(durations))
}

#' Calculate Minimum Intensity (voiced speech only)
#' @keywords internal
calculate_minimum_intensity_pladdrr <- function(df, apply_calibration, calibration) {
  # Load and concatenate all soft phonation samples
  sounds <- lapply(seq_len(nrow(df)), function(i) {
    sound <- pladdrr::Sound(df$absolute_file_path[i])
    # Extract time window if specified
    if (df$start[i] > 0 || df$end[i] > 0) {
      duration <- sound$.cpp$duration
      end_time <- if (df$end[i] > 0) df$end[i] else duration
      sound <- sound$extract_part(df$start[i], end_time, "rectangular", 1.0, FALSE)
    }
    return(sound)
  })

  # Concatenate sounds
  combined <- if (length(sounds) == 1) {
    sounds[[1]]
  } else if (exists("sound_concatenate_all", where = asNamespace("pladdrr"))) {
    pladdrr::sound_concatenate_all(sounds)
  } else {
    # Fallback: iterative concatenation
    result <- sounds[[1]]
    for (i in 2:length(sounds)) {
      result <- result$concatenate(sounds[[i]])
    }
    result
  }

  # Use ultra API if available
  if (exists("calculate_minimum_intensity_ultra", where = asNamespace("pladdrr"))) {
    minimum_intensity <- pladdrr::calculate_minimum_intensity_ultra(
      combined,
      min_pitch = 70,
      max_pitch = 600,
      time_step = 0,
      subtract_mean = TRUE
    )
  } else {
    # Fallback: manual calculation
    # Extract pitch → point process → voiced intervals → intensity
    pitch_ptr <- pladdrr::to_pitch_cc_direct(combined, 0, 70, 600, 15, TRUE, 0.03, 0.8, 0.01, 0.35, 0.14)
    pitch <- pladdrr::Pitch(.xptr = pitch_ptr)

    pp_ptr <- pladdrr::to_point_process_from_sound_and_pitch(combined, pitch_ptr)

    # Calculate intensity on voiced regions
    intensity_ptr <- pladdrr::to_intensity_direct(combined, 70, 0, TRUE)
    intensity <- pladdrr::Intensity(.xptr = intensity_ptr)

    minimum_intensity <- intensity$get_minimum(0, 0, "parabolic")
  }

  if (is.na(minimum_intensity)) {
    stop("No voiced intervals found in soft phonation samples")
  }

  # Apply calibration if requested
  if (apply_calibration) {
    minimum_intensity <- minimum_intensity + calibration
  }

  return(minimum_intensity)
}

#' Calculate Maximum Fundamental Frequency
#' @keywords internal
calculate_maximum_f0_pladdrr <- function(df) {
  # Load and concatenate all high pitch samples
  sounds <- lapply(seq_len(nrow(df)), function(i) {
    sound <- pladdrr::Sound(df$absolute_file_path[i])
    # Extract time window
    if (df$start[i] > 0 || df$end[i] > 0) {
      duration <- sound$.cpp$duration
      end_time <- if (df$end[i] > 0) df$end[i] else duration
      sound <- sound$extract_part(df$start[i], end_time, "rectangular", 1.0, FALSE)
    }
    return(sound)
  })

  # Concatenate
  combined <- if (length(sounds) == 1) {
    sounds[[1]]
  } else if (exists("sound_concatenate_all", where = asNamespace("pladdrr"))) {
    pladdrr::sound_concatenate_all(sounds)
  } else {
    result <- sounds[[1]]
    for (i in 2:length(sounds)) {
      result <- result$concatenate(sounds[[i]])
    }
    result
  }

  # Use ultra API if available
  if (exists("calculate_f0_stats_ultra", where = asNamespace("pladdrr"))) {
    maximum_f0 <- pladdrr::calculate_f0_stats_ultra(
      combined,
      stat = "max",
      min_pitch = 70,
      max_pitch = 1300,  # DSI uses 1300 Hz for high voices
      time_step = 0,
      voicing_threshold = 0.8  # CRITICAL: 0.8, not default 0.45
    )
  } else {
    # Fallback: manual calculation
    pitch_ptr <- pladdrr::to_pitch_cc_direct(
      combined, 0, 70, 1300, 15, TRUE, 0.03, 0.8, 0.01, 0.35, 0.14
    )
    pitch <- pladdrr::Pitch(.xptr = pitch_ptr)
    maximum_f0 <- pitch$get_maximum(0, 0, "hertz", "parabolic")
  }

  return(maximum_f0)
}

#' Calculate Jitter PPQ5 from Stable Vowel
#' @keywords internal
calculate_jitter_ppq5_pladdrr <- function(df) {
  # Load and concatenate all stable vowel samples
  sounds <- lapply(seq_len(nrow(df)), function(i) {
    sound <- pladdrr::Sound(df$absolute_file_path[i])
    # Extract time window
    if (df$start[i] > 0 || df$end[i] > 0) {
      duration <- sound$.cpp$duration
      end_time <- if (df$end[i] > 0) df$end[i] else duration
      sound <- sound$extract_part(df$start[i], end_time, "rectangular", 1.0, FALSE)
    }
    return(sound)
  })

  # Concatenate
  combined <- if (length(sounds) == 1) {
    sounds[[1]]
  } else if (exists("sound_concatenate_all", where = asNamespace("pladdrr"))) {
    pladdrr::sound_concatenate_all(sounds)
  } else {
    result <- sounds[[1]]
    for (i in 2:length(sounds)) {
      result <- result$concatenate(sounds[[i]])
    }
    result
  }

  # Use last 3 seconds or entire duration if shorter (DSI protocol)
  duration <- combined$.cpp$duration
  if (duration > 3) {
    combined <- combined$extract_part(duration - 3, duration, "rectangular", 1.0, FALSE)
  }

  # Use ultra API if available
  if (exists("get_voice_quality_ultra", where = asNamespace("pladdrr"))) {
    vq_metrics <- pladdrr::get_voice_quality_ultra(
      combined,
      metrics = "jitter",
      min_pitch = 70,
      max_pitch = 600,
      time_step = 0.005
    )
    jitter_ppq5 <- vq_metrics$jitter_ppq5 * 100
  } else {
    # Fallback: manual calculation
    pitch_ptr <- pladdrr::to_pitch_cc_direct(combined, 0, 70, 600, 15, TRUE, 0.03, 0.45, 0.01, 0.35, 0.14)
    pp_ptr <- pladdrr::to_point_process_from_sound_and_pitch(combined, pitch_ptr)
    pp <- pladdrr::PointProcess(.xptr = pp_ptr)

    jitter_ppq5 <- pp$get_jitter_ppq5(0, 0, 0.0001, 0.02, 1.3) * 100
  }

  return(jitter_ppq5)
}

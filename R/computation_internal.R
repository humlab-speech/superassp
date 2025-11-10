#' Internal Computation Functions (Pure DSP)
#'
#' @description
#' This file provides internal computation functions that perform pure DSP
#' processing without I/O operations. These functions are used internally
#' by public DSP functions and are designed to be testable independently.
#'
#' All functions in this file:
#' - Accept audio data in memory (not file paths)
#' - Return computation results (no file writing)
#' - Are deterministic (same input = same output)
#' - Have no side effects (no messages, warnings, or file operations)
#'
#' @name computation_internal
#' @keywords internal
NULL

#' Compute COVAREP voice quality parameters (internal)
#'
#' Pure computation function for COVAREP voice quality extraction.
#' Separates DSP computation from I/O for testability.
#'
#' @param audio_samples Numeric vector of audio samples
#' @param sample_rate Integer sample rate in Hz
#' @param f0 Numeric scalar or vector with F0 values (optional)
#' @param gci Numeric vector with GCI timestamps (optional)
#' @param gci_in_samples Logical indicating if GCI is in samples vs seconds
#' @param covarep_module Python module object (reticulate)
#'
#' @return List with voice quality parameters or NULL on failure
#' @keywords internal
#' @noRd
.compute_covarep_vq_internal <- function(audio_samples, sample_rate,
                                         f0 = NULL, gci = NULL,
                                         gci_in_samples = FALSE,
                                         covarep_module) {

  # Compute IAIF to get glottal flow
  iaif_result <- covarep_module$glottal$iaif_optimized$iaif_optimized(
    x = audio_samples,
    fs = as.integer(sample_rate)
  )

  glottal_flow <- as.numeric(iaif_result[[1]])
  glottal_derivative <- as.numeric(iaif_result[[2]])

  # Check IAIF success
  if (length(glottal_flow) == 0) {
    return(NULL)
  }

  # Prepare F0 parameter for Python
  py_f0 <- NULL
  if (!is.null(f0)) {
    if (length(f0) == 1) {
      py_f0 <- as.numeric(f0)
    } else {
      np <- reticulate::import("numpy", delay_load = FALSE)
      py_f0 <- np$array(as.numeric(f0))
    }
  }

  # Prepare GCI parameter for Python
  py_gci <- NULL
  if (!is.null(gci)) {
    gci_vec <- as.numeric(gci)
    if (!gci_in_samples) {
      gci_vec <- gci_vec * sample_rate
    }
    np <- reticulate::import("numpy", delay_load = FALSE)
    py_gci <- np$array(gci_vec, dtype = "int32")
  }

  # Compute voice quality parameters
  vq_result <- covarep_module$voice_quality$compute_vq_params$compute_vq_params(
    glottal_flow = glottal_flow,
    glottal_derivative = glottal_derivative,
    fs = as.integer(sample_rate),
    f0 = py_f0,
    gci = py_gci
  )

  # Convert to R list
  param_list <- list(
    NAQ = as.numeric(vq_result$NAQ),
    QOQ = as.numeric(vq_result$QOQ),
    H1H2 = as.numeric(vq_result$H1H2),
    HRF = as.numeric(vq_result$HRF),
    PSP = as.numeric(vq_result$PSP)
  )

  return(param_list)
}

#' Compute Phonet phonological posteriors (internal)
#'
#' Pure computation function for Phonet analysis.
#' Separates DSP computation from I/O for testability.
#'
#' @param audio_samples Numeric vector of audio samples (16 kHz mono)
#' @param classes Character vector of phonological classes to extract
#' @param phonet_instance Python Phonet object
#'
#' @return List with time, phoneme, and posterior probabilities or NULL on failure
#' @keywords internal
#' @noRd
.compute_phonet_internal <- function(audio_samples, classes, phonet_instance) {

  # Phonet expects 16 kHz audio
  # Note: Caller must ensure audio is at correct sample rate

  # Convert to temporary WAV file for Phonet
  # (Phonet currently requires file input)
  temp_wav <- tempfile(fileext = ".wav")
  on.exit(unlink(temp_wav), add = TRUE)

  # Write audio to temp file
  av::av_audio_convert(
    audio = audio_samples,
    output = temp_wav,
    format = "wav",
    channels = 1,
    sample_rate = 16000
  )

  # Run Phonet extraction
  # feat_file="" returns DataFrame directly
  df_result <- phonet_instance$get_phon_wav(
    audio_file = temp_wav,
    feat_file = "",
    plot_flag = FALSE
  )

  # Convert pandas DataFrame to R list
  result <- list(
    time = as.numeric(df_result$time),
    phoneme = as.character(df_result$phoneme)
  )

  # Add phonological posteriors
  for (class_name in classes) {
    if (class_name != "all" && class_name %in% names(df_result)) {
      result[[class_name]] <- as.numeric(df_result[[class_name]])
    }
  }

  # If "all" was requested, add all available classes
  if ("all" %in% classes) {
    for (col in names(df_result)) {
      if (!col %in% c("time", "phoneme") && !col %in% names(result)) {
        result[[col]] <- as.numeric(df_result[[col]])[1]
      }
    }
  }

  if (length(result$time) == 0) {
    return(NULL)
  }

  return(result)
}

#' Summarize Phonet posteriors for JSTF (internal)
#'
#' Computes summary statistics (mean/SD) from time-series posterior probabilities
#' for storage in JSTF format.
#'
#' @param phonet_result List from .compute_phonet_internal()
#'
#' @return List with mean and SD for each posterior, or NULL on error
#' @keywords internal
#' @noRd
.summarize_phonet_posteriors <- function(phonet_result) {

  if (is.null(phonet_result) || !is.null(phonet_result$error)) {
    return(NULL)
  }

  summary_result <- list()

  for (name in names(phonet_result)) {
    # Skip non-posterior fields
    if (name %in% c("time", "phoneme", "file", "error")) {
      next
    }

    # Compute statistics for numeric posteriors
    if (is.numeric(phonet_result[[name]]) && length(phonet_result[[name]]) > 0) {
      summary_result[[paste0(name, "_mean")]] <- mean(phonet_result[[name]], na.rm = TRUE)
      summary_result[[paste0(name, "_sd")]] <- sd(phonet_result[[name]], na.rm = TRUE)
    }
  }

  return(summary_result)
}

#' Load audio for computation (internal helper)
#'
#' Loads audio file segment into memory for processing.
#' Wrapper around av_load_for_python for consistency.
#'
#' @param file_path Character string with file path
#' @param start_time Numeric start time in seconds (0 = beginning)
#' @param end_time Numeric end time in seconds (0 = end of file)
#'
#' @return List with audio_np (numpy array), sample_rate, and samples (R vector)
#' @keywords internal
#' @noRd
.load_audio_for_computation <- function(file_path, start_time = 0, end_time = 0) {

  audio_data <- av_load_for_python(
    file_path,
    start_time = if (start_time > 0) start_time else NULL,
    end_time = if (end_time > 0) end_time else NULL
  )

  # Add R samples vector for convenience
  audio_data$samples <- as.numeric(audio_data$audio_np)

  return(audio_data)
}

#' Resample audio to target rate (internal)
#'
#' Resamples audio to a target sample rate if needed.
#'
#' @param audio_samples Numeric vector of audio samples
#' @param from_rate Integer current sample rate
#' @param to_rate Integer target sample rate
#'
#' @return Numeric vector with resampled audio
#' @keywords internal
#' @noRd
.resample_audio_internal <- function(audio_samples, from_rate, to_rate) {

  if (from_rate == to_rate) {
    return(audio_samples)
  }

  # Use av package for resampling via temporary file
  # (More robust than pure R resampling)
  temp_in <- tempfile(fileext = ".wav")
  temp_out <- tempfile(fileext = ".wav")
  on.exit(unlink(c(temp_in, temp_out)), add = TRUE)

  # Write input
  av::av_audio_convert(
    audio = audio_samples,
    output = temp_in,
    format = "wav",
    sample_rate = from_rate,
    channels = 1
  )

  # Resample
  av::av_audio_convert(
    audio = temp_in,
    output = temp_out,
    format = "wav",
    sample_rate = to_rate,
    channels = 1
  )

  # Read resampled
  audio_data <- av::read_audio_bin(temp_out, channels = 1)
  samples <- as.numeric(audio_data) / 2147483647.0  # INT32_MAX normalization

  return(samples)
}

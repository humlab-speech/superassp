#' Check if pladdrr package is available
#' @keywords internal
pladdrr_available <- function() {
  requireNamespace("pladdrr", quietly = TRUE)
}

#' Load audio file as pladdrr Sound object
#'
#' @description
#' Simple wrapper to load audio files directly with pladdrr. Handles time
#' windowing via pladdrr's native \code{extract_part()} method.
#'
#' @param file_path Character path to audio file (WAV, AIFF, FLAC, MP3, etc.)
#' @param start_time Numeric start time in seconds (default: 0 = beginning)
#' @param end_time Numeric end time in seconds (default: 0 = end of file)
#' @param window_type Character window shape for part extraction (default: "Gaussian1")
#' @param relative_width Numeric relative width of window (default: 1.0)
#'
#' @return A pladdrr Sound R6 object
#'
#' @details
#' pladdrr reads audio files directly via its C++ backend (no av needed).
#' Supported formats: WAV, AIFF, FLAC, MP3, NIST (via native readers),
#' or any format via av package fallback.
#'
#' If time windowing is requested (start_time > 0 or end_time > 0),
#' the function uses pladdrr's \code{extract_part()} method with appropriate
#' windowing.
#'
#' @examples
#' \dontrun{
#' # Load entire file
#' sound <- av_load_for_pladdrr("speech.wav")
#'
#' # Load with time windowing (1.0 to 3.0 seconds)
#' sound <- av_load_for_pladdrr("speech.wav", start_time = 1.0, end_time = 3.0)
#'
#' # Use with pladdrr functions
#' pitch <- sound$to_pitch_cc(time_step = 0.01, pitch_floor = 75, pitch_ceiling = 600)
#' }
#'
av_load_for_pladdrr <- function(file_path,
                                start_time = 0.0,
                                end_time = 0.0,
                                window_type = "Gaussian1",
                                relative_width = 1.0) {

  # Check if pladdrr is available
  if (!pladdrr_available()) {
    cli::cli_abort(c(
      "x" = "pladdrr package not available",
      "i" = "Install with: install_pladdrr()"
    ))
  }

  # Check file exists
  if (!file.exists(file_path)) {
    cli::cli_abort("File not found: {file_path}")
  }

  # Load sound with pladdrr (reads files directly)
  # Fallback: if pladdrr cannot read the format, transcode to WAV via av
  sound <- tryCatch(
    pladdrr::Sound(file_path),
    error = function(e) NULL
  )

  if (is.null(sound)) {
    if (!requireNamespace("av", quietly = TRUE)) {
      stop("pladdrr could not read '", basename(file_path),
           "' and the 'av' package is not available as a fallback.")
    }
    tmp_wav <- tempfile(fileext = ".wav")
    on.exit(unlink(tmp_wav), add = TRUE)
    av::av_audio_convert(file_path, tmp_wav, verbose = FALSE)
    sound <- pladdrr::Sound(tmp_wav)
  }

  # Handle time windowing if requested
  if (start_time > 0 || end_time > 0) {
    duration <- sound$.cpp$duration
    actual_end <- if (end_time > 0) end_time else duration

    # Pin the original Sound to prevent GC — extract_part() returns a new
    # object that may share C++ memory with the original.
    original_sound <- sound
    sound <- sound$extract_part(
      from_time = start_time,
      to_time = actual_end,
      window_shape = window_type,
      relative_width = relative_width,
      preserve_times = TRUE
    )
    on.exit(invisible(original_sound), add = TRUE)
  }

  return(sound)
}


#' Get external pointer from pladdrr R6 object
#'
#' @description
#' Extract the underlying C pointer from a pladdrr R6 object.
#' Used internally for passing objects to direct API functions.
#'
#' @param obj pladdrr R6 object (Sound, Pitch, Formant, etc.)
#'
#' @return External pointer to the Praat C object
#'
#' @details
#' pladdrr uses R6 function factory pattern where all objects have a .xptr field
#' containing the C pointer.
#'
#' @examples
#' \dontrun{
#' sound <- pladdrr::Sound("audio.wav")
#' ptr <- get_pladdrr_ptr(sound)
#' }
#'
#' @keywords internal
get_pladdrr_ptr <- function(obj) {
  if (is.null(obj$.xptr)) {
    stop("Object does not have .xptr field")
  }
  return(obj$.xptr)
}


#' Convert pladdrr data frame to superassp format
#'
#' @description
#' Converts data frames from pladdrr (long format with time/frequency/value)
#' to superassp's wide format expected for AsspDataObj tracks.
#'
#' @param df Data frame from pladdrr R6 object's \code{as_data_frame()} method
#' @param type Character indicating data type: "pitch", "formant", "intensity"
#' @param n_formants Integer number of formants (default: 5, for formant data only)
#'
#' @return Data frame in wide format suitable for AsspDataObj
#'
#' @details
#' Transformations by type:
#'
#' **Pitch**: pladdrr returns (time, frequency) → superassp expects (time, F0)
#'
#' **Formant**: pladdrr returns (time, formant_number, frequency, bandwidth) →
#'   superassp expects (time, fm1, fm2, ..., fm5, bw1, bw2, ..., bw5)
#'
#' **Intensity**: pladdrr returns (time, intensity) → superassp expects same
#'
#' @examples
#' \dontrun{
#' sound <- pladdrr::Sound("audio.wav")
#' pitch <- sound$to_pitch_cc()
#' pitch_df <- pitch$as_data_frame()
#' wide_df <- pladdrr_df_to_superassp(pitch_df, type = "pitch")
#' }
#'
pladdrr_df_to_superassp <- function(df, type = c("pitch", "formant", "intensity"),
                                    n_formants = 5) {
  type <- match.arg(type)

  if (type == "pitch") {
    # pladdrr: (time, frequency) → superassp: (time, F0)
    # Replace 0 or undefined with NA
    df$frequency[df$frequency == 0 | !is.finite(df$frequency)] <- NA
    colnames(df) <- c("time", "F0")
    return(df)

  } else if (type == "formant") {
    # pladdrr long format: (time, formant_number, frequency, bandwidth)
    # superassp wide format: (time, fm1, fm2, ..., fm5, bw1, bw2, ..., bw5)

    # Get unique times
    times <- sort(unique(df$time))

    # Initialize result with time column
    result <- data.frame(time = times)

    # Add formant frequency columns (fm1, fm2, ...)
    for (i in 1:n_formants) {
      freq_col <- paste0("fm", i)
      bw_col <- paste0("bw", i)

      # Extract data for this formant
      f_data <- df[df$formant_number == i, ]

      # Match times
      time_idx <- match(result$time, f_data$time)

      result[[freq_col]] <- f_data$frequency[time_idx]
      result[[bw_col]] <- f_data$bandwidth[time_idx]
    }

    # Replace 0 with NA
    for (col in names(result)[-1]) {
      result[[col]][result[[col]] == 0 | !is.finite(result[[col]])] <- NA
    }

    return(result)

  } else if (type == "intensity") {
    # pladdrr: (time, intensity) → superassp: (time, intensity)
    colnames(df) <- c("time", "intensity")
    return(df)

  } else {
    stop("Unknown type: ", type)
  }
}

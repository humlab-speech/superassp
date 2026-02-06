#' Convert av audio data to pladdrr Sound object
#'
#' @description
#' Converts audio data loaded by \code{\link[av]{read_audio_bin}} to a
#' pladdrr Sound object for in-memory processing. This eliminates the
#' need for temporary WAV files.
#'
#' @param audio_data Integer vector of audio samples from \code{av::read_audio_bin}.
#'   Expected to be INT32 format in range \[-2147483648, 2147483647\].
#' @param sample_rate Integer sample rate in Hz. If NULL, extracted from
#'   audio_data attributes.
#' @param normalize Logical. Normalize audio to \[-1, 1\] range (default: TRUE).
#'   Required for pladdrr compatibility.
#'
#' @return A pladdrr Sound R6 object that can be used with all pladdrr
#'   functions (pitch analysis, formant extraction, etc.)
#'
#' @details
#' This function performs the conversion:
#' \enumerate{
#'   \item Extracts sample rate from av audio data attributes
#'   \item Converts INT32 samples to numeric normalized to \[-1, 1\]
#'   \item Constructs pladdrr::Sound object in memory
#' }
#'
#' The resulting Sound object can be used with any pladdrr/Praat function:
#' - \code{sound$to_pitch()} - Pitch analysis
#' - \code{sound$to_formant()} - Formant analysis
#' - \code{sound$to_intensity()} - Intensity analysis
#' - \code{sound$to_spectrogram()} - Spectrogram
#' - Direct API: \code{to_pitch_cc_direct()}, \code{to_formant_direct()}, etc.
#'
#' @examples
#' \dontrun{
#' # Load audio using av
#' audio_data <- av::read_audio_bin("speech.wav", channels = 1)
#'
#' # Convert to pladdrr Sound (in memory, no file I/O)
#' sound <- av_to_pladdrr_sound(audio_data)
#'
#' # Use with pladdrr functions
#' pitch <- sound$to_pitch(time_step = 0.01, pitch_floor = 75, pitch_ceiling = 600)
#' formants <- sound$to_formant(max_formant = 5500)
#' intensity <- sound$to_intensity()
#'
#' # Extract values
#' mean_pitch <- pitch$get_mean(0, 0, "Hertz")
#' }
#'
#' @seealso
#' \code{\link{av_load_for_pladdrr}} for complete workflow with time windowing
#'
#' @export
av_to_pladdrr_sound <- function(audio_data,
                                sample_rate = NULL,
                                normalize = TRUE) {

  # Check if pladdrr is available
  if (!pladdrr_available()) {
    cli::cli_abort(c(
      "x" = "pladdrr package not available",
      "i" = "Install with: install_pladdrr()"
    ))
  }

  # Extract sample rate from attributes if not provided
  if (is.null(sample_rate)) {
    sample_rate <- attr(audio_data, "sample_rate")
    if (is.null(sample_rate)) {
      cli::cli_abort("sample_rate must be provided or present in audio_data attributes")
    }
  }

  # Convert INT32 to float
  # av::read_audio_bin returns INT32 in range [-2147483648, 2147483647]
  INT32_MAX <- 2147483647

  if (normalize) {
    # Normalize to [-1, 1] range (required for pladdrr)
    audio_float <- as.numeric(audio_data) / INT32_MAX
  } else {
    audio_float <- as.numeric(audio_data)
  }

  # Create pladdrr Sound object
  # pladdrr works directly with R numeric vectors (no numpy needed!)
  sound <- pladdrr::Sound(
    audio_float,
    sampling_frequency = as.integer(sample_rate)
  )

  return(sound)
}


#' Load audio for pladdrr processing
#'
#' @description
#' Complete workflow to load audio using av package and convert to pladdrr
#' Sound object. Handles time windowing and format conversion automatically.
#'
#' @param file_path Character path to audio file (any format supported by av)
#' @param start_time Numeric start time in seconds (default: NULL = file start)
#' @param end_time Numeric end time in seconds (default: NULL = file end)
#' @param channels Integer number of channels (default: 1 = mono)
#' @param target_sample_rate Integer target sample rate in Hz (default: NULL = original)
#'
#' @return A pladdrr Sound R6 object ready for analysis
#'
#' @details
#' This function combines \code{av::read_audio_bin} and
#' \code{av_to_pladdrr_sound} into a single convenient call.
#'
#' Advantages over file-based approach:
#' \itemize{
#'   \item No temporary files created (pure in-memory processing)
#'   \item Time windowing handled by av (efficient)
#'   \item Supports all media formats (WAV, MP3, MP4, video, etc.)
#'   \item Automatic resampling if needed
#'   \item Cleaner code, no temp file cleanup required
#'   \item Native R integration (no Python dependency)
#' }
#'
#' Advantages of pladdrr over parselmouth:
#' \itemize{
#'   \item Pure R/C implementation (no Python/reticulate)
#'   \item R6 object-oriented interface
#'   \item Direct access to Praat C library
#'   \item No numpy conversion overhead
#'   \item Better R integration
#' }
#'
#' @examples
#' \dontrun{
#' # Load entire file
#' sound <- av_load_for_pladdrr("speech.wav")
#'
#' # Load with time windowing
#' sound <- av_load_for_pladdrr("speech.wav",
#'                              start_time = 1.0,
#'                              end_time = 3.0)
#'
#' # Load and resample
#' sound <- av_load_for_pladdrr("speech.mp3",
#'                              target_sample_rate = 16000)
#'
#' # Use with pladdrr
#' pitch <- sound$to_pitch()
#' formants <- sound$to_formant()
#' }
#'
#' @seealso
#' \code{\link{av_to_pladdrr_sound}} for the conversion step only
#'
#' @export
av_load_for_pladdrr <- function(file_path,
                                start_time = NULL,
                                end_time = NULL,
                                channels = 1,
                                target_sample_rate = NULL) {

  # Load audio using av package
  audio_data <- av::read_audio_bin(
    audio = file_path,
    start_time = start_time,
    end_time = end_time,
    channels = channels,
    sample_rate = target_sample_rate
  )

  # Convert to pladdrr Sound
  sound <- av_to_pladdrr_sound(audio_data)

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
#' containing the external pointer to the underlying Praat C object.
#'
#' @keywords internal
#' @export
get_pladdrr_ptr <- function(obj) {
  # pladdrr R6 objects have .xptr field
  if (!is.null(obj$.xptr)) {
    return(obj$.xptr)
  }
  
  # Fallback: check for .cpp field (older pladdrr versions)
  if (!is.null(obj$.cpp)) {
    if (!is.null(obj$.cpp$.xptr)) {
      return(obj$.cpp$.xptr)
    }
  }
  
  stop("Could not find .xptr in pladdrr object")
}


#' Convert pladdrr data frame to superassp format
#'
#' @description
#' Converts pladdrr's long-format data frames to the wide format expected
#' by superassp/AsspDataObj.
#'
#' @param df Data frame from pladdrr (e.g., pitch$as_data_frame())
#' @param type Type of data: "pitch", "formant", "intensity", etc.
#'
#' @return Data frame in superassp format
#'
#' @details
#' pladdrr returns data in long format:
#' \itemize{
#'   \item Pitch: time, frequency
#'   \item Formant: time, formant, frequency, bandwidth
#'   \item Intensity: time, intensity
#' }
#'
#' superassp expects wide format:
#' \itemize{
#'   \item Pitch: time, pitch
#'   \item Formant: time, F1, F2, F3, F4, F5, B1, B2, B3, B4, B5
#'   \item Intensity: time, intensity
#' }
#'
#' @keywords internal
#' @export
pladdrr_df_to_superassp <- function(df, type = "pitch") {
  
  if (type == "pitch") {
    # pladdrr: (time, frequency) → superassp: (time, pitch)
    # Replace 0 and undefined with NA
    df[[2]][df[[2]] == 0 | !is.finite(df[[2]])] <- NA
    colnames(df) <- c("time", "pitch")
    return(df)
    
  } else if (type == "formant") {
    # pladdrr: (time, formant, frequency, bandwidth) → superassp: (time, F1, F2, ...)
    # Reshape long to wide format
    times <- unique(df$time)
    formant_nums <- unique(df$formant)
    max_formant <- max(formant_nums)
    
    # Initialize result
    result <- data.frame(time = times)
    
    for (fn in formant_nums) {
      f_data <- df[df$formant == fn, ]
      freq_col <- paste0("F", fn)
      bw_col <- paste0("B", fn)
      
      # Match times
      time_idx <- match(times, f_data$time)
      
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

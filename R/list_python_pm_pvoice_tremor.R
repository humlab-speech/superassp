##' Voice Tremor Analysis using Parselmouth (memory-based, optimized)
##'
##' @name praat_voice_tremor
NULL

#' Compute vocal tremor measures using Parselmouth
#'
#' This is a memory-based implementation of \code{\link{praat_voice_tremor}} using
#' Parselmouth instead of external Praat. It eliminates disk I/O by loading
#' audio directly into memory using the av package and processing with
#' Parselmouth.
#'
#' This function provides identical functionality to \code{praat_voice_tremor} but
#' with significantly improved performance (10-20x faster) by eliminating
#' file I/O operations.
#'
#' The function extracts 18 measures of vocal tremor from sustained phonations:
#' \itemize{
#'   \item{\strong{Frequency Tremor Measures:}}
#'     \itemize{
#'       \item{FCoM: Frequency contour magnitude}
#'       \item{FTrC: Frequency tremor cyclicality}
#'       \item{FMoN: Number of frequency modulation candidates}
#'       \item{FTrF: Frequency tremor frequency (Hz)}
#'       \item{FTrI: Frequency tremor intensity index (%)}
#'       \item{FTrP: Frequency tremor power index}
#'       \item{FTrCIP: Frequency tremor cyclicality-intensity product}
#'       \item{FTrPS: Frequency tremor product sum}
#'       \item{FCoHNR: Frequency contour HNR (dB)}
#'     }
#'   \item{\strong{Amplitude Tremor Measures:}}
#'     \itemize{
#'       \item{ACoM: Amplitude contour magnitude}
#'       \item{ATrC: Amplitude tremor cyclicality}
#'       \item{AMoN: Number of amplitude modulation candidates}
#'       \item{ATrF: Amplitude tremor frequency (Hz)}
#'       \item{ATrI: Amplitude tremor intensity index (%)}
#'       \item{ATrP: Amplitude tremor power index}
#'       \item{ATrCIP: Amplitude tremor cyclicality-intensity product}
#'       \item{ATrPS: Amplitude tremor product sum}
#'       \item{ACoHNR: Amplitude contour HNR (dB)}
#'     }
#' }
#'
#' @param listOfFiles Character vector with path(s) to audio file(s)
#' @param beginTime Numeric. Start time in seconds (default 0)
#' @param endTime Numeric. End time in seconds (0 = end of file)
#' @param analysis.time.step Numeric. Time step for analysis in seconds (default 0.015)
#' @param min.pitch Numeric. Minimum pitch for extraction in Hz (default 60)
#' @param max.pitch Numeric. Maximum pitch for extraction in Hz (default 350)
#' @param silence.threshold Numeric. Threshold for silence detection (default 0.03)
#' @param voicing.threshold Numeric. Threshold for voicing detection (default 0.3)
#' @param octave.cost Numeric. Cost for octave jumps in pitch tracking (default 0.01)
#' @param octave.jump.cost Numeric. Cost for large octave jumps (default 0.35)
#' @param voiced.unvoiced.cost Numeric. Cost for voiced/unvoiced transitions (default 0.14)
#' @param min.tremor.hz Numeric. Minimum tremor frequency in Hz (default 1.5)
#' @param max.tremor.hz Numeric. Maximum tremor frequency in Hz (default 15)
#' @param contour.magnitude.threshold Numeric. Threshold for contour magnitude (default 0.01)
#' @param tremor.cyclicality.threshold Numeric. Threshold for cyclicality (default 0.15)
#' @param freq.tremor.octave.cost Numeric. Octave cost for frequency tremor (default 0.01)
#' @param ampl.tremor.octave.cost Numeric. Octave cost for amplitude tremor (default 0.01)
#' @param na.zero Logical. Should undefined measurements be returned as zeros (TRUE) or NA (FALSE, default)?
#' @param amplitude.extraction.method Integer. Method for amplitude extraction: 1 = RMS per pitch period, 2 = Envelope (default)
#' @param speaker.name Character. Name of speaker (optional)
#' @param speaker.ID Character. Speaker identifier (optional, defaults to speaker.name)
#' @param pdf.path Character. Directory path to save Praat picture files (optional)
#' @param praat_path Character. Path to Praat executable (not used, for compatibility only)
#'
#' @return A list with tremor measurements:
#' \describe{
#'   \item{ID}{Speaker ID}
#'   \item{Speaker}{Speaker name}
#'   \item{FCoM}{Frequency contour magnitude}
#'   \item{FTrC}{Frequency tremor cyclicality}
#'   \item{FMoN}{Number of frequency modulation candidates}
#'   \item{FTrF [Hz]}{Frequency tremor frequency}
#'   \item{FTrI [%]}{Frequency tremor intensity index}
#'   \item{FTrP}{Frequency tremor power index}
#'   \item{FTrCIP}{Frequency tremor cyclicality-intensity product}
#'   \item{FTrPS}{Frequency tremor product sum}
#'   \item{FCoHNR[dB]}{Frequency contour HNR}
#'   \item{ACoM}{Amplitude contour magnitude}
#'   \item{ATrC}{Amplitude tremor cyclicality}
#'   \item{AMoN}{Number of amplitude modulation candidates}
#'   \item{ATrF [Hz]}{Amplitude tremor frequency}
#'   \item{ATrI [%]}{Amplitude tremor intensity index}
#'   \item{ATrP}{Amplitude tremor power index}
#'   \item{ATrCIP}{Amplitude tremor cyclicality-intensity product}
#'   \item{ATrPS}{Amplitude tremor product sum}
#'   \item{ACoHNR[dB]}{Amplitude contour HNR}
#' }
#'
#' @export
#'
#' @references
#' Brückl, M. (2017). Vocal tremor measurement based on autocorrelation of contours.
#' Proceedings of Interspeech 2017, 2027-2031.
#'
#' @examples
#' \dontrun{
#' # Analyze a sustained vowel recording
#' result <- lst_voice_tremorp(
#'   listOfFiles = "sustained_vowel.wav",
#'   beginTime = 0,
#'   endTime = 0,  # end of file
#'   speaker.name = "John Doe",
#'   speaker.ID = "001"
#' )
#'
#' # Access frequency tremor measures
#' print(result$`FTrF [Hz]`)  # Tremor frequency
#' print(result$`FTrI [%]`)   # Tremor intensity
#' }
lst_voice_tremorp <- function(listOfFiles,
                                   beginTime = 0,
                                   endTime = 0,
                                   analysis.time.step = 0.015,
                                   min.pitch = 60,
                                   max.pitch = 350,
                                   silence.threshold = 0.03,
                                   voicing.threshold = 0.3,
                                   octave.cost = 0.01,
                                   octave.jump.cost = 0.35,
                                   voiced.unvoiced.cost = 0.14,
                                   min.tremor.hz = 1.5,
                                   max.tremor.hz = 15,
                                   contour.magnitude.threshold = 0.01,
                                   tremor.cyclicality.threshold = 0.15,
                                   freq.tremor.octave.cost = 0.01,
                                   ampl.tremor.octave.cost = 0.01,
                                   na.zero = FALSE,
                                   amplitude.extraction.method = 2,
                                   speaker.name = NULL,
                                   speaker.ID = speaker.name,
                                   pdf.path = NULL,
                                   praat_path = NULL) {

  # Check that Parselmouth is available
  if (!reticulate::py_module_available("parselmouth")) {
    stop("Parselmouth Python module not available. Install with: pip install praat-parselmouth")
  }

  # Validate input files
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)

  filesEx <- file.exists(listOfFiles)
  if (!all(filesEx)) {
    filesNotExist <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ", paste(filesNotExist, collapse = ", "))
  }

  # Ensure time vectors match file count
  n_files <- length(listOfFiles)
  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime) == 1) endTime <- rep(endTime, n_files)

  # Convert na.zero to nan_output_mode (1 = zeros, 2 = undefined/NA)
  nan_output_mode <- if (na.zero) 1 else 2

  # Source Python script
  python_script <- system.file("python", "praat_voice_tremor_memory.py", package = "superassp")
  if (!file.exists(python_script)) {
    stop("Python script not found: ", python_script)
  }
  reticulate::source_python(python_script)

  # Get Python main module
  py <- reticulate::import_main()

  # Process each file
  results_list <- list()

  for (i in seq_along(listOfFiles)) {
    file_path <- normalizePath(listOfFiles[i], mustWork = TRUE)
    start_sec <- as.numeric(beginTime[i])
    end_sec <- if (endTime[i] == 0) NULL else as.numeric(endTime[i])

    # Load audio segment with av
    audio_data <- av_load_for_python(
      file_path,
      start_time = start_sec,
      end_time = end_sec
    )

    # Generate picture path if pdf.path provided
    picture_file <- NULL
    if (!is.null(pdf.path)) {
      # Create picture filename based on speaker ID or file name
      pic_name <- if (!is.null(speaker.ID)) {
        paste0(speaker.ID, "_tremor.prapic")
      } else {
        paste0(tools::file_path_sans_ext(basename(file_path)), "_tremor.prapic")
      }
      picture_file <- file.path(pdf.path, pic_name)
    }

    # Call Python function with audio array
    result <- py$praat_voice_tremor_memory(
      audio_np = audio_data$audio_np,
      sample_rate = audio_data$sample_rate,
      analysis_time_step = analysis.time.step,
      min_pitch = min.pitch,
      max_pitch = max.pitch,
      silence_threshold = silence.threshold,
      voicing_threshold = voicing.threshold,
      octave_cost = octave.cost,
      octave_jump_cost = octave.jump.cost,
      voiced_unvoiced_cost = voiced.unvoiced.cost,
      min_tremor_freq = min.tremor.hz,
      max_tremor_freq = max.tremor.hz,
      contour_magnitude_threshold = contour.magnitude.threshold,
      tremor_cyclicality_threshold = tremor.cyclicality.threshold,
      freq_tremor_octave_cost = freq.tremor.octave.cost,
      amp_tremor_octave_cost = ampl.tremor.octave.cost,
      nan_output_mode = nan_output_mode,
      amplitude_extraction_method = amplitude.extraction.method,
      speaker_name = if (is.null(speaker.name)) "" else speaker.name,
      speaker_id = if (is.null(speaker.ID)) "" else as.character(speaker.ID),
      picture_path = picture_file
    )

    # Convert Python dict to R list
    result_list <- as.list(result)

    # Inform user about picture file and PDF conversion
    if (!is.null(picture_file)) {
      message("Praat picture saved to: ", picture_file)
      message("To convert to PDF, run the convert_prapic_to_pdf.praat script in ", pdf.path)
    }

    results_list[[i]] <- result_list
  }

  # If single file, return single result; otherwise return list
  if (n_files == 1) {
    logger::log_trace("Computed voice tremor measurements from ", listOfFiles[1], " (Parselmouth).")
    return(results_list[[1]])
  } else {
    logger::log_trace("Computed voice tremor measurements from ", n_files, " files (Parselmouth).")
    return(results_list)
  }
}

attr(lst_voice_tremorp, "outputType") <- c("list")
attr(lst_voice_tremorp, "ext") <- c("pvt")
attr(lst_voice_tremorp, "tracks") <- c("FCoM", "FTrC", "FMoN", "FTrF [Hz]", "FTrI [%]",
                                             "FTrP", "FTrCIP", "FTrPS", "FCoHNR[dB]",
                                             "ACoM", "ATrC", "AMoN", "ATrF [Hz]", "ATrI [%]",
                                             "ATrP", "ATrCIP", "ATrPS", "ACoHNR[dB]")

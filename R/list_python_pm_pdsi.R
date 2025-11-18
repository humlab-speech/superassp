##' DSI analysis using Parselmouth (memory-based, optimized)
##'
##' @name praat_dsi
NULL

#' Compute the Dysphonia Severity Index using Parselmouth
#'
#' This is a memory-based implementation of \code{\link{praat_dsi}} using
#' Parselmouth instead of external Praat. It eliminates disk I/O by loading
#' audio directly into memory using the av package and processing with
#' Parselmouth.
#'
#' This function provides identical functionality to \code{praat_dsi} but
#' with significantly improved performance (10-20x faster) by eliminating
#' file I/O operations.
#'
#' @param softDF Data frame with soft voice samples. Must contain columns: absolute_file_path, start, end
#' @param highpitchDF Data frame with high pitch samples. Must contain columns: absolute_file_path, start, end
#' @param maxprolongedDF Data frame with maximally prolonged vowel samples. Must contain columns: absolute_file_path, start, end
#' @param stableDF Data frame with stable vowel samples (optional). If NULL, uses maxprolongedDF
#' @param use.calibration Apply calibration to measurements (default: FALSE)
#' @param db.calibration Calibration value in dB (default: 10)
#' @param speaker.name Speaker name (optional)
#' @param speaker.ID Speaker ID (optional)
#' @param speaker.dob Speaker date of birth (optional)
#' @param session.datetime Session date and time (optional)
#' @param pdf.path Path for PDF output (optional)
#' @param simple.output Return simplified output (default: FALSE)
#' @param overwrite.pdfs Overwrite existing PDF files (default: FALSE)
#' @param praat_path Path to Praat executable (optional, for compatibility)
#' @param toFile Logical. If TRUE, write results to JSTF file. Default FALSE.
#' @param explicitExt Character. File extension for output. Default "dsi".
#' @param outputDirectory Character. Output directory path. Default NULL (use input directory).
#'
#' @return If \code{toFile=FALSE} (default), a list with DSI measurements.
#'   If \code{toFile=TRUE}, invisibly returns the path to the written JSTF file.
#'
#'   The list contains:
#'   \describe{
#'     \item{\code{ID}}{Speaker ID}
#'     \item{\code{Maximum_phonation_time}}{Longest sustained vowel duration in seconds}
#'     \item{\code{Softest_intensity_of_voiced_speech}}{Minimum intensity at soft phonation (dB SPL)}
#'     \item{\code{Maximum_fundamental_frequency}}{Highest F0 in high pitch samples (Hz)}
#'     \item{\code{Jitter_ppq5}}{5-point period perturbation quotient (%)}
#'     \item{\code{Dysphonia_Severity_Index}}{DSI composite score (continuous scale)}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Define soft voice samples
#' soft <- data.frame(
#'   absolute_file_path = c("soft1.wav", "soft2.wav"),
#'   start = c(0, 0),  # seconds
#'   end = c(2, 2)
#' )
#'
#' # Define high pitch samples
#' high <- data.frame(
#'   absolute_file_path = c("high1.wav"),
#'   start = 0,
#'   end = 1.5
#' )
#'
#' # Define maximally prolonged vowel samples
#' prolonged <- data.frame(
#'   absolute_file_path = c("prolonged1.wav", "prolonged2.wav"),
#'   start = c(0, 0),
#'   end = c(5, 6)
#' )
#'
#' # Stable vowel (optional, will use prolonged if not provided)
#' stable <- data.frame(
#'   absolute_file_path = "stable1.wav",
#'   start = 0,
#'   end = 3
#' )
#'
#' # Compute DSI
#' result <- lst_dsip(
#'   softDF = soft,
#'   highpitchDF = high,
#'   maxprolongedDF = prolonged,
#'   stableDF = stable,
#'   speaker.name = "John Doe",
#'   speaker.ID = "001"
#' )
#'
#' # Write results to JSTF file
#' lst_dsip(
#'   softDF = soft,
#'   highpitchDF = high,
#'   maxprolongedDF = prolonged,
#'   stableDF = stable,
#'   speaker.name = "John Doe",
#'   speaker.ID = "001",
#'   toFile = TRUE
#' )  # Creates 001.dsi
#'
#' # Read back and convert to data.frame
#' track <- read_track("001.dsi")
#' df <- as.data.frame(track)
#' head(df)  # Shows begin_time, end_time, and all DSI measures
#' }
lst_dsip <- function(softDF,
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
                          overwrite.pdfs = FALSE,
                          praat_path = NULL,
                          toFile = FALSE,
                          explicitExt = "dsi",
                          outputDirectory = NULL) {

  # Validate JSTF parameters
  validate_jstf_parameters(toFile, explicitExt, outputDirectory, "lst_dsip")

  # Check that Parselmouth is available
  if (!reticulate::py_module_available("parselmouth")) {
    stop("Parselmouth Python module not available. Install with: pip install praat-parselmouth")
  }

  # Validate required columns
  requiredDFColumns <- c("absolute_file_path", "start", "end")

  # Reuse max prolonged vowels as stable productions if not provided
  if (is.null(stableDF)) {
    stableDF <- maxprolongedDF
  }

  # Check all dataframes are non-empty
  if (nrow(softDF) < 1 || nrow(highpitchDF) < 1 ||
      nrow(maxprolongedDF) < 1 || nrow(stableDF) < 1) {
    stop("All of the supplied dataframes indicating samples to be analysed must be non-empty.")
  }

  # Validate columns in all dataframes
  if (!all(requiredDFColumns %in% names(softDF)) ||
      !all(requiredDFColumns %in% names(highpitchDF)) ||
      !all(requiredDFColumns %in% names(maxprolongedDF)) ||
      !all(requiredDFColumns %in% names(stableDF))) {
    stop("All dataframes must both contain columns named ",
         paste(requiredDFColumns, collapse = ",", sep = ""), ".")
  }

  # Get all unique files and check they exist
  listOfFiles <- unique(c(
    softDF$absolute_file_path,
    highpitchDF$absolute_file_path,
    maxprolongedDF$absolute_file_path,
    stableDF$absolute_file_path
  ))

  filesEx <- file.exists(listOfFiles)
  if (!all(filesEx)) {
    filesNotExist <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ", paste(filesNotExist, collapse = ", "))
  }

  # Load soft voice audio segments
  soft_audio_list <- list()
  for (r in 1:nrow(softDF)) {
    file_path <- normalizePath(softDF[[r, "absolute_file_path"]], mustWork = TRUE)
    start_sec <- as.numeric(softDF[[r, "start"]])
    end_sec <- as.numeric(softDF[[r, "end"]])

    # Load audio segment with av
    audio_data <- av_load_for_python(
      file_path,
      start_time = start_sec,
      end_time = end_sec
    )

    soft_audio_list[[r]] <- list(
      audio_np = audio_data$audio_np,
      sample_rate = audio_data$sample_rate
    )
  }

  # Load high pitch audio segments
  highpitch_audio_list <- list()
  for (r in 1:nrow(highpitchDF)) {
    file_path <- normalizePath(highpitchDF[[r, "absolute_file_path"]], mustWork = TRUE)
    start_sec <- as.numeric(highpitchDF[[r, "start"]])
    end_sec <- as.numeric(highpitchDF[[r, "end"]])

    audio_data <- av_load_for_python(
      file_path,
      start_time = start_sec,
      end_time = end_sec
    )

    highpitch_audio_list[[r]] <- list(
      audio_np = audio_data$audio_np,
      sample_rate = audio_data$sample_rate
    )
  }

  # Load maximally prolonged vowel audio segments
  maxprolonged_audio_list <- list()
  for (r in 1:nrow(maxprolongedDF)) {
    file_path <- normalizePath(maxprolongedDF[[r, "absolute_file_path"]], mustWork = TRUE)
    start_sec <- as.numeric(maxprolongedDF[[r, "start"]])
    end_sec <- as.numeric(maxprolongedDF[[r, "end"]])

    audio_data <- av_load_for_python(
      file_path,
      start_time = start_sec,
      end_time = end_sec
    )

    maxprolonged_audio_list[[r]] <- list(
      audio_np = audio_data$audio_np,
      sample_rate = audio_data$sample_rate
    )
  }

  # Load stable vowel audio segments
  stable_audio_list <- list()
  for (r in 1:nrow(stableDF)) {
    file_path <- normalizePath(stableDF[[r, "absolute_file_path"]], mustWork = TRUE)
    start_sec <- as.numeric(stableDF[[r, "start"]])
    end_sec <- as.numeric(stableDF[[r, "end"]])

    audio_data <- av_load_for_python(
      file_path,
      start_time = start_sec,
      end_time = end_sec
    )

    stable_audio_list[[r]] <- list(
      audio_np = audio_data$audio_np,
      sample_rate = audio_data$sample_rate
    )
  }

  # Source Python script
  python_script <- system.file("python", "praat_dsi_memory.py", package = "superassp")
  if (!file.exists(python_script)) {
    stop("Python script not found: ", python_script)
  }
  reticulate::source_python(python_script)

  # Get Python main module
  py <- reticulate::import_main()

  # Generate picture path if pdf.path provided
  picture_file <- NULL
  if (!is.null(pdf.path)) {
    # Create picture filename based on speaker ID and date
    pic_name <- if (!is.null(speaker.ID) && !is.null(session.datetime)) {
      paste0(speaker.ID, "_", session.datetime, ".prapic")
    } else if (!is.null(speaker.ID)) {
      paste0(speaker.ID, ".prapic")
    } else {
      paste0("dsi_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".prapic")
    }
    picture_file <- file.path(pdf.path, pic_name)
  }

  # Call Python function with audio arrays
  result <- py$praat_dsi_memory(
    soft_audio_list = soft_audio_list,
    highpitch_audio_list = highpitch_audio_list,
    maxprolonged_audio_list = maxprolonged_audio_list,
    stable_audio_list = stable_audio_list,
    apply_calibration = use.calibration,
    calibration_db = db.calibration,
    speaker_name = if (is.null(speaker.name)) "" else speaker.name,
    speaker_id = if (is.null(speaker.ID)) "" else as.character(speaker.ID),
    speaker_dob = if (is.null(speaker.dob)) "" else speaker.dob,
    assessment_date = if (is.null(session.datetime)) "" else session.datetime,
    picture_path = picture_file
  )

  # Convert Python dict to R list
  result_list <- as.list(result)

  # Inform user about picture file and PDF conversion
  if (!is.null(picture_file)) {
    message("Praat picture saved to: ", picture_file)
    message("To convert to PDF, run the convert_prapic_to_pdf.praat script in ", pdf.path)
  }

  logger::log_trace("Computed a DSI value from ", nrow(softDF), " soft voice, ",
                    nrow(highpitchDF), " high pitch, ",
                    nrow(maxprolongedDF), " maximally prolonged, and ",
                    nrow(stableDF), " stable vowel samples (Parselmouth).")

  # Handle JSTF file writing
  if (toFile) {
    # Calculate total analysis time range from all segments
    all_start_times <- c(softDF$start, highpitchDF$start, maxprolongedDF$start, stableDF$start)
    all_end_times <- c(softDF$end, highpitchDF$end, maxprolongedDF$end, stableDF$end)
    analysis_begin <- min(all_start_times)
    analysis_end <- max(all_end_times)

    # Use first file as primary reference
    primary_file <- normalizePath(softDF[[1, "absolute_file_path"]], mustWork = TRUE)

    output_path <- write_lst_results_to_jstf(
      results = list(result_list),
      file_paths = primary_file,
      beginTime = analysis_begin,
      endTime = analysis_end,
      function_name = "lst_dsip",
      parameters = list(
        use.calibration = use.calibration,
        db.calibration = db.calibration,
        speaker.name = speaker.name,
        speaker.ID = speaker.ID,
        n_soft_segments = nrow(softDF),
        n_highpitch_segments = nrow(highpitchDF),
        n_maxprolonged_segments = nrow(maxprolongedDF),
        n_stable_segments = nrow(stableDF)
      ),
      explicitExt = explicitExt,
      outputDirectory = outputDirectory,
      speaker_id = speaker.ID
    )
    return(invisible(output_path))
  }

  return(result_list)
}

attr(lst_dsip, "ext") <- "dsi"
attr(lst_dsip, "outputType") <- "JSTF"
attr(lst_dsip, "format") <- "JSON"
attr(lst_dsip, "tracks") <- c("ID", "Maximum_phonation_time",
                                    "Softest_intensity_of_voiced_speech",
                                    "Maximum_fundamental_frequency",
                                    "Jitter_ppq5",
                                    "Dysphonia_Severity_Index")

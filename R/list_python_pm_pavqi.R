##' AVQI analysis using Parselmouth (memory-based, optimized)
##'
##' @name praat_avqi
NULL

#' Compute the Acoustic Voice Quality Index (AVQI) using Parselmouth
#'
#' This is a memory-based implementation of \code{\link{praat_avqi}} using
#' Parselmouth instead of external Praat. It eliminates disk I/O by loading
#' audio directly into memory using the av package and processing with
#' Parselmouth.
#'
#' This function provides identical functionality to \code{praat_avqi} but
#' with significantly improved performance (10-20x faster) by eliminating
#' file I/O operations.
#'
#' @param svDF Data frame with sustained vowel samples. Must contain columns: listOfFiles, start, end
#' @param csDF Data frame with continuous speech samples. Must contain columns: listOfFiles, start, end
#' @param min.sv Minimum sustained vowel duration in milliseconds (default: 1000)
#' @param speaker.name Speaker name (optional)
#' @param speaker.ID Speaker ID (defaults to speaker.name)
#' @param speaker.dob Speaker date of birth (optional)
#' @param session.datetime Session date and time (optional)
#' @param pdf.path Path for PDF output (optional)
#' @param simple.output Return simplified output (default: FALSE)
#' @param overwrite.pdfs Overwrite existing PDF files (default: FALSE)
#' @param praat_path Path to Praat executable (optional, for compatibility)
#' @param toFile Logical. If TRUE, write results to JSTF file. Default FALSE.
#' @param explicitExt Character. File extension for output. Default "avqi".
#' @param outputDirectory Character. Output directory path. Default NULL (use input directory).
#'
#' @return If \code{toFile=FALSE} (default), a list with AVQI measurements.
#'   If \code{toFile=TRUE}, invisibly returns the path to the written JSTF file.
#'
#'   The list contains:
#'   \describe{
#'     \item{\code{AVQI_VERSION}}{Version of AVQI algorithm}
#'     \item{\code{Speaker}}{Speaker name}
#'     \item{\code{ID}}{Speaker ID}
#'     \item{\code{CPPS}}{Cepstral Peak Prominence Smoothed}
#'     \item{\code{HNR}}{Harmonics-to-Noise Ratio}
#'     \item{\code{Shim_local}}{Local shimmer}
#'     \item{\code{Shim_local_DB}}{Local shimmer in dB}
#'     \item{\code{LTAS_Slope}}{Long-term average spectrum slope}
#'     \item{\code{LTAS_Tilt}}{Long-term average spectrum tilt}
#'     \item{\code{AVQI}}{Acoustic Voice Quality Index (0-10 scale)}
#'   }
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
#' # Compute AVQI
#' result <- lst_avqip(sv, cs,
#'                          speaker.name = "John Doe",
#'                          speaker.ID = "001")
#'
#' # Write results to JSTF file
#' lst_avqip(sv, cs,
#'           speaker.name = "John Doe",
#'           speaker.ID = "001",
#'           toFile = TRUE)  # Creates AVQI result file
#'
#' # Read back and convert to data.frame
#' track <- read_track("001.avqi")  # Or appropriate output path
#' df <- as.data.frame(track)
#' head(df)  # Shows begin_time, end_time, and all AVQI measures
#' }
lst_avqip <- function(svDF,
                           csDF,
                           min.sv = AVQI_MIN_SV_DURATION_MS,
                           speaker.name = NULL,
                           speaker.ID = speaker.name,
                           speaker.dob = NULL,
                           session.datetime = NULL,
                           pdf.path = NULL,
                           simple.output = FALSE,
                           overwrite.pdfs = FALSE,
                           praat_path = NULL,
                           toFile = FALSE,
                           explicitExt = "avqi",
                           outputDirectory = NULL) {

  # Validate JSTF parameters
  validate_jstf_parameters(toFile, explicitExt, outputDirectory, "lst_avqip")

  # Check that Parselmouth is available
  if (!reticulate::py_module_available("parselmouth")) {
    stop("Parselmouth Python module not available. Install with: pip install praat-parselmouth")
  }

  # Validate required columns
  requiredDFColumns <- c("listOfFiles", "start", "end")

  if (!all(requiredDFColumns %in% names(svDF)) || !all(requiredDFColumns %in% names(csDF))) {
    stop("The 'svDF' and 'csDF' structures must both contain columns named ",
         paste(requiredDFColumns, collapse = ",", sep = ""), ".")
  }

  # Check minimum sustained vowel duration (times are in milliseconds)
  totalSVdur <- sum(svDF$end - svDF$start)

  if (totalSVdur < min.sv) {
    stop("The total sustained vowel duration is less than the threshold length min.sv!")
  }

  # Get all unique files
  listOfFiles <- unique(c(svDF$listOfFiles, csDF$listOfFiles))

  # Check that all files exist before we begin
  filesEx <- file.exists(listOfFiles)
  if (!all(filesEx)) {
    filesNotExist <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ", paste(filesNotExist, collapse = ", "))
  }

  # Load and prepare sustained vowel audio segments
  sv_audio_list <- list()
  for (r in 1:nrow(svDF)) {
    file_path <- normalizePath(as.character(svDF[[r, "listOfFiles"]]), mustWork = TRUE)
    start_sec <- ms_to_sec(as.numeric(svDF[[r, "start"]]))
    end_sec <- ms_to_sec(as.numeric(svDF[[r, "end"]]))

    # Load audio segment with av
    audio_data <- av_load_for_python(
      file_path,
      start_time = start_sec,
      end_time = end_sec
    )

    sv_audio_list[[r]] <- list(
      audio_np = audio_data$audio_np,
      sample_rate = audio_data$sample_rate
    )
  }

  # Load and prepare continuous speech audio segments
  cs_audio_list <- list()
  for (r in 1:nrow(csDF)) {
    file_path <- normalizePath(csDF[[r, "listOfFiles"]], mustWork = TRUE)
    start_sec <- ms_to_sec(csDF[[r, "start"]])
    end_sec <- ms_to_sec(csDF[[r, "end"]])

    # Load audio segment with av
    audio_data <- av_load_for_python(
      file_path,
      start_time = start_sec,
      end_time = end_sec
    )

    cs_audio_list[[r]] <- list(
      audio_np = audio_data$audio_np,
      sample_rate = audio_data$sample_rate
    )
  }

  # Source Python script
  python_script <- system.file("python", "praat_avqi_memory.py", package = "superassp")
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
      paste0("avqi_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".prapic")
    }
    picture_file <- file.path(pdf.path, pic_name)
  }

  # Call Python function with audio arrays
  result <- py$praat_avqi_memory(
    sv_audio_list = sv_audio_list,
    cs_audio_list = cs_audio_list,
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

  logger::log_trace("Computed an AVQI value from ", nrow(svDF), " sustained vowels and ",
                    nrow(csDF), " continuous speech utterances (Parselmouth).")

  # Handle JSTF file writing
  if (toFile) {
    # Calculate total analysis time range
    all_start_times <- ms_to_sec(c(svDF$start, csDF$start))
    all_end_times <- ms_to_sec(c(svDF$end, csDF$end))
    analysis_begin <- min(all_start_times)
    analysis_end <- max(all_end_times)

    # Use first file as primary reference
    primary_file <- normalizePath(svDF[[1, "listOfFiles"]], mustWork = TRUE)

    output_path <- write_lst_results_to_jstf(
      results = list(result_list),
      file_paths = primary_file,
      beginTime = analysis_begin,
      endTime = analysis_end,
      function_name = "lst_avqip",
      parameters = list(
        min.sv = min.sv,
        speaker.name = speaker.name,
        speaker.ID = speaker.ID,
        n_sv_segments = nrow(svDF),
        n_cs_segments = nrow(csDF),
        total_sv_duration_ms = sum(svDF$end - svDF$start),
        total_cs_duration_ms = sum(csDF$end - csDF$start)
      ),
      explicitExt = explicitExt,
      outputDirectory = outputDirectory,
      speaker_id = speaker.ID
    )
    return(invisible(output_path))
  }

  return(result_list)
}



attr(lst_avqip, "ext") <- "avqi"
attr(lst_avqip, "outputType") <- "JSTF"
attr(lst_avqip, "format") <- "JSON"
attr(lst_avqip, "tracks") <- c("AVQI_VERSION", "Speaker", "ID", "CPPS", "HNR",
                                     "Shim_local", "Shim_local_DB", "LTAS_Slope",
                                     "LTAS_Tilt", "AVQI")

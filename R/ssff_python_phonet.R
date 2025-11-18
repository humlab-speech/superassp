#' Extract phonological posteriors using Phonet (list format)
#'
#' Phonet \insertCite{vasquez2019phonet}{superassp} uses Bidirectional Gated Recurrent Unit (BGRU)
#' neural networks to compute posterior probabilities of phonological classes from speech audio.
#' The models analyze Mel-filterbank features to predict 18 phonological classes based on mode
#' and manner of articulation.
#'
#' This function returns phonological posteriors as a list/data.frame (not SSFF format), making it
#' suitable for feature extraction and data analysis. For SSFF track output compatible with emuR,
#' use \code{\link{trk_phonet}} instead.
#'
#' @param listOfFiles A character vector of file paths to audio files, or an AVAudio S7 object.
#' @param classes A character vector specifying which phonological classes to extract.
#'   Options: "vocalic", "consonantal", "back", "anterior", "open", "close", "nasal",
#'   "stop", "continuant", "lateral", "flap", "trill", "voice", "strident", "labial",
#'   "dental", "velar", "pause", or "all" (for all 17 classes).
#'   Default is c("vocalic", "consonantal", "nasal", "stop").
#' @param beginTime The start time of the section of the sound file that should be processed (in seconds).
#'   Default is 0.0 (start of file).
#' @param endTime The end time of the section of the sound file that should be processed (in seconds).
#'   Default is 0.0 (end of file).
#' @param conda.env The name of the conda environment in which Python and its required packages
#'   are stored. Defaults to NULL, which uses the default environment or RETICULATE_PYTHON.
#' @param verbose If TRUE (default), print progress messages.
#' @param toFile Logical. If TRUE, write results to JSTF file. Default FALSE.
#' @param explicitExt Character. File extension for output. Default "pho".
#' @param outputDirectory Character. Output directory path. Default NULL (use input directory).
#'
#' @return
#'   If \code{toFile=FALSE} (default), a list (or list of lists for multiple files) with the following elements:
#'   \itemize{
#'     \item \strong{time}: Time vector (seconds) at 10ms intervals
#'     \item \strong{phoneme}: Recognized phonemes (character vector)
#'     \item \strong{[class]}: Posterior probability (0-1) for each requested phonological class
#'     \item \strong{file}: Original file path
#'   }
#'
#'   If \code{toFile=TRUE}, invisibly returns the path(s) to the written JSTF file(s).
#'
#'   Note: The JSTF output contains time-averaged posterior probabilities (mean values across the segment).
#'   For time-series posteriors suitable for emuR, use \code{\link{trk_phonet}} instead.
#'
#' @details
#' \strong{Phonet Algorithm:}
#' \enumerate{
#'   \item \strong{Audio Loading}: File loaded and resampled to 16 kHz if needed
#'   \item \strong{Feature Extraction}: 33 Mel-filterbanks + energy (25ms window, 10ms shift)
#'   \item \strong{BGRU Processing}: 2-layer bidirectional GRU (128 units each)
#'   \item \strong{Multi-task Prediction}: 17 binary classifiers (one per phonological class)
#'   \item \strong{Post-processing}: Mask correction smooths posterior probabilities
#' }
#'
#' \strong{Specifications:}
#' \itemize{
#'   \item \strong{Sample Rate}: 16 kHz (auto-resamples if needed)
#'   \item \strong{Frame Rate}: 10ms intervals (100 Hz)
#'   \item \strong{Features}: 33 Mel-filterbanks + energy (34D)
#'   \item \strong{Window}: 25ms Hamming window
#'   \item \strong{Model}: Pre-trained BGRU on Spanish speech
#'   \item \strong{Output}: Posterior probabilities (0-1) per class
#' }
#'
#' \strong{Installation Requirements:}
#'
#' Phonet requires Python with TensorFlow, Keras, and speech processing libraries. Install with:
#'
#' \code{install_phonet()}
#'
#' Or manually:
#'
#' \code{reticulate::py_install("git+https://github.com/jcvasquezc/phonet.git")}
#'
#' \strong{Phonological Classes:}
#'
#' \itemize{
#'   \item \strong{Vowel features}: vocalic, back, anterior, open, close
#'   \item \strong{Consonant features}: consonantal, nasal, stop, continuant, lateral, flap, trill, voice, strident
#'   \item \strong{Place of articulation}: labial, dental, velar
#'   \item \strong{Other}: pause (silence/pauses)
#' }
#'
#' @examples
#' \dontrun{
#' # Install Phonet first
#' install_phonet()
#'
#' # Basic usage - extract specific classes
#' result <- lst_phonet("speech.wav",
#'                      classes = c("nasal", "stop", "vocalic"))
#'
#' # Access posteriors
#' nasal_post <- result$nasal  # Probability of nasal sounds
#' time <- result$time         # Time stamps
#' phonemes <- result$phoneme  # Recognized phonemes
#'
#' # Extract all 17 phonological classes
#' result_all <- lst_phonet("speech.wav", classes = "all")
#'
#' # Time windowing (extract from 1.0 to 3.0 seconds)
#' result_window <- lst_phonet("long_recording.wav",
#'                             classes = c("vocalic", "consonantal"),
#'                             beginTime = 1.0,
#'                             endTime = 3.0)
#'
#' # Batch processing multiple files
#' files <- c("file1.wav", "file2.wav", "file3.mp3")
#' results <- lst_phonet(files, classes = c("nasal", "stop"))
#'
#' # Convert to data frame for analysis
#' df <- as.data.frame(result)
#' library(ggplot2)
#' ggplot(df, aes(x = time, y = nasal)) +
#'   geom_line() +
#'   labs(title = "Nasal Posterior Probability")
#'
#' # Write results to JSTF file
#' lst_phonet("speech.wav",
#'            classes = c("vocalic", "consonantal", "nasal", "stop"),
#'            toFile = TRUE)  # Creates speech.pho
#'
#' # Read back and convert to data.frame
#' track <- read_track("speech.pho")
#' df <- as.data.frame(track)
#' head(df)  # Shows begin_time, end_time, and phonological features
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso
#' \code{\link{install_phonet}} to install Phonet dependencies
#' \code{\link{phonet_available}} to check availability
#' \code{\link{phonet_info}} to get configuration information
#' \code{\link{trk_phonet}} for SSFF track output (time-series format)
#'
#' @export
lst_phonet <- function(listOfFiles,
                       classes = c("vocalic", "consonantal", "nasal", "stop"),
                       beginTime = 0.0,
                       endTime = 0.0,
                       conda.env = NULL,
                       verbose = TRUE,
                       toFile = FALSE,
                       explicitExt = "pho",
                       outputDirectory = NULL) {

  # Validate JSTF parameters
  validate_jstf_parameters(toFile, explicitExt, outputDirectory, "lst_phonet")

  # Check conda environment
  if (!is.null(conda.env) && !conda.env %in% reticulate::conda_list()$name) {
    stop("The conda environment ", conda.env, " does not exist.")
  }

  # Check Phonet availability
  if (!phonet_available()) {
    stop(
      "lst_phonet() requires Python with phonet, tensorflow, keras, and speech processing libraries.\n\n",
      "Install with: install_phonet()\n",
      "Or manually: reticulate::py_install('git+https://github.com/jcvasquezc/phonet.git')\n",
      call. = FALSE
    )
  }

  # Validate phonological classes
  valid_classes <- c(
    "vocalic", "consonantal", "back", "anterior", "open", "close",
    "nasal", "stop", "continuant", "lateral", "flap", "trill", "voice",
    "strident", "labial", "dental", "velar", "pause", "all"
  )

  if (!all(classes %in% valid_classes)) {
    invalid <- classes[!classes %in% valid_classes]
    stop(
      "Invalid phonological class(es): ", paste(invalid, collapse = ", "), "\n",
      "Valid classes: ", paste(valid_classes, collapse = ", "),
      call. = FALSE
    )
  }

  # Create file/time dataframe
  tryCatch({
    fileBeginEnd <- data.frame(
      listOfFiles = listOfFiles,
      beginTime = beginTime,
      endTime = endTime,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    stop("The beginTime and endTime must either be a single value or the same length as listOfFiles",
         call. = FALSE)
  })

  # Check that all files exist
  filesEx <- file.exists(listOfFiles)
  if (!all(filesEx)) {
    filesNotExists <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s): ", paste(filesNotExists, collapse = ", "),
         call. = FALSE)
  }

  # Import phonet module
  phonet <- reticulate::import("phonet", delay_load = FALSE)

  # Results list
  results <- list()

  # Initialize Phonet instance (reuse for all files)
  if (verbose) {
    message("Initializing Phonet with classes: ", paste(classes, collapse = ", "))
  }

  phon <- phonet$Phonet(classes)

  # Process each file
  for (i in seq_len(nrow(fileBeginEnd))) {
    origSoundFile <- normalizePath(fileBeginEnd[i, "listOfFiles"], mustWork = TRUE)
    bt <- fileBeginEnd[i, "beginTime"]
    et <- fileBeginEnd[i, "endTime"]

    if (verbose) {
      message(sprintf("Processing file %d/%d: %s", i, nrow(fileBeginEnd), basename(origSoundFile)))
    }

    tryCatch({
      # Create temporary WAV file at 16 kHz (Phonet requirement)
      temp_wav <- tempfile(fileext = ".wav")
      on.exit(unlink(temp_wav), add = TRUE)

      # Use av to convert/resample to 16 kHz WAV with time windowing
      av::av_audio_convert(
        audio = origSoundFile,
        output = temp_wav,
        format = "wav",
        channels = DEFAULT_CHANNELS,
        sample_rate = PHONET_SAMPLE_RATE,
        start_time = if (bt > 0) bt else NULL,
        total_time = if (et > bt && et > 0) (et - bt) else NULL
      )

      # Run Phonet extraction
      # feat_file="" returns DataFrame directly, plot_flag=FALSE for speed
      df_result <- phon$get_phon_wav(
        audio_file = temp_wav,
        feat_file = "",
        plot_flag = FALSE
      )

      # Convert pandas DataFrame to R list
      result <- list(
        time = as.numeric(df_result$time),
        phoneme = as.character(df_result$phoneme),
        file = origSoundFile
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
            result[[col]] <- as.numeric(df_result[[col]])
          }
        }
      }

      results[[i]] <- result

      if (verbose) {
        message(sprintf("  → Extracted %d time frames with %d posteriors",
                       length(result$time), length(result) - 3))
      }

    }, error = function(e) {

      warning(format_processing_error(origSoundFile, safe_error_message(e), "Phonet phonological analysis"),
              call. = FALSE)

      results[[i]] <- list(
        time = numeric(0),
        phoneme = character(0),
        file = origSoundFile,

        error = safe_error_message(e)

      )
    })
  }


  # Handle JSTF file writing
  if (toFile) {
    # Compute summary statistics (mean/SD posteriors) for JSTF
    # JSTF stores scalar values, so we average the time-series posteriors
    summary_results <- lapply(results, function(result) {
      if (is.null(result) || !is.null(result$error)) {
        return(NULL)
      }

      summary_result <- list()
      for (name in names(result)) {
        if (name %in% c("time", "phoneme", "file", "error")) {
          next  # Skip non-posterior fields
        }
        if (is.numeric(result[[name]]) && length(result[[name]]) > 0) {
          summary_result[[paste0(name, "_mean")]] <- mean(result[[name]], na.rm = TRUE)
          summary_result[[paste0(name, "_sd")]] <- sd(result[[name]], na.rm = TRUE)
        }
      }
      summary_result
    })

    # Prepare parameters with frame count from first valid result
    first_valid <- Find(function(x) !is.null(x), results)
    params <- list(
      classes = paste(classes, collapse = ", "),
      n_frames = if (!is.null(first_valid)) length(first_valid$time) else 0,
      n_posteriors = if (!is.null(summary_results[[1]])) length(summary_results[[1]]) / 2 else 0,
      conda.env = conda.env
    )

    output_paths <- write_lst_results_to_jstf(
      results = summary_results,
      file_paths = fileBeginEnd$listOfFiles,
      beginTime = fileBeginEnd$beginTime,
      endTime = fileBeginEnd$endTime,
      function_name = "lst_phonet",
      parameters = params,
      explicitExt = explicitExt,
      outputDirectory = outputDirectory
    )
    return(invisible(output_paths))
  }

  # Return single result for single file, otherwise list
  if (length(results) == 1) {
    return(results[[1]])
  } else {
    return(results)
  }
}


# Set function attributes
attr(lst_phonet, "ext") <- "pho"
attr(lst_phonet, "outputType") <- "JSTF"
attr(lst_phonet, "format") <- "JSON"


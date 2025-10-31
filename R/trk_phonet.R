#' Track phonological posteriors using Phonet BGRU models
#'
#' Phonet \insertCite{vasquez2019phonet}{superassp} uses Bidirectional Gated Recurrent Unit (BGRU)
#' neural networks to compute posterior probabilities of phonological classes from speech audio as
#' continuous time-series data. This function returns SSFF track objects suitable for time-aligned
#' phonological analysis and visualization in the emuR framework.
#'
#' Unlike \code{\link{lst_phonet}} which returns list format, \code{trk_phonet} returns SSFF track
#' objects that can be written to disk and loaded with emuR for time-aligned annotation.
#'
#' @param listOfFiles A character vector of file paths to audio files, or an AVAudio S7 object.
#' @param classes A character vector specifying which phonological classes to track.
#'   Options: "vocalic", "consonantal", "back", "anterior", "open", "close", "nasal",
#'   "stop", "continuant", "lateral", "flap", "trill", "voice", "strident", "labial",
#'   "dental", "velar", "pause", or "all" (for all 18 classes).
#'   Default is c("vocalic", "consonantal", "nasal", "stop").
#' @param beginTime The start time of the section of the sound file that should be processed (in seconds).
#'   Default is 0.0 (start of file).
#' @param endTime The end time of the section of the sound file that should be processed (in seconds).
#'   Default is 0.0 (end of file).
#' @param explicitExt The file extension for the output SSFF file. Default is "phn" (phonet).
#' @param outputDirectory The directory where output files should be written.
#'   If NULL (default), files are written to the same directory as the input files.
#' @param toFile If TRUE (default), write results to SSFF files and return file count.
#'   If FALSE, return AsspDataObj for single file processing.
#' @param verbose If TRUE (default), print progress messages.
#' @param conda.env The name of the conda environment in which Python and its required packages
#'   are stored. Defaults to NULL, which uses the default environment or RETICULATE_PYTHON.
#'
#' @return
#'   If \code{toFile = TRUE}: Returns the number of successfully processed files.
#'   If \code{toFile = FALSE}: Returns an AsspDataObj with tracks for each phonological class:
#'   \itemize{
#'     \item \strong{[class]}: Posterior probability (0-1) for each requested phonological class (REAL32 format)
#'   }
#'
#' @details
#' \strong{Phonet Tracking Algorithm:}
#' \enumerate{
#'   \item \strong{Audio Loading}: File loaded via av package (supports all formats)
#'   \item \strong{Resampling}: Audio resampled to 16 kHz (Phonet requirement)
#'   \item \strong{Feature Extraction}: 33 Mel-filterbanks + energy (25ms window, 10ms shift)
#'   \item \strong{BGRU Processing}: 2-layer bidirectional GRU (128 units each)
#'   \item \strong{Multi-task Prediction}: Binary classifiers for each phonological class
#'   \item \strong{Post-processing}: Mask correction smooths posterior probabilities
#'   \item \strong{SSFF Conversion}: Results formatted as SSFF tracks at 100 Hz
#' }
#'
#' \strong{Specifications:}
#' \itemize{
#'   \item \strong{Sample Rate}: 16 kHz (auto-resamples if needed)
#'   \item \strong{Frame Rate}: 100 Hz (10ms intervals)
#'   \item \strong{Features}: 33 Mel-filterbanks + energy (34D)
#'   \item \strong{Window}: 25ms Hamming window
#'   \item \strong{Model}: Pre-trained BGRU on Spanish speech
#'   \item \strong{Output Format}: SSFF with REAL32 tracks
#' }
#'
#' \strong{Installation Requirements:}
#'
#' Phonet requires Python with TensorFlow, Keras 2.x API, and speech processing libraries:
#'
#' \code{install_phonet()}
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
#' \strong{Use Cases:}
#' \itemize{
#'   \item Time-aligned phonological feature analysis
#'   \item Integration with emuR databases
#'   \item Phonetic segmentation and labeling
#'   \item Articulatory feature tracking
#'   \item Speech disorder analysis (dysarthria, apraxia)
#' }
#'
#' @examples
#' \dontrun{
#' # Install Phonet first
#' install_phonet()
#'
#' # Basic usage - track specific classes, return AsspDataObj
#' result <- trk_phonet("speech.wav",
#'                      classes = c("nasal", "stop", "vocalic"),
#'                      toFile = FALSE)
#'
#' # Access tracks
#' nasal_track <- result$nasal  # Posterior probabilities over time
#' time_stamps <- seq(attr(result, "startTime"),
#'                    by = 1/attr(result, "sampleRate"),
#'                    length.out = nrow(result$nasal))
#'
#' # Track all 18 phonological classes and write to SSFF files
#' trk_phonet("speech.wav", classes = "all", toFile = TRUE)
#' # Creates: speech.phn (SSFF file with 18 tracks)
#'
#' # Time windowing (track from 1.0 to 3.0 seconds)
#' result <- trk_phonet("long_recording.wav",
#'                      classes = c("vocalic", "consonantal"),
#'                      beginTime = 1.0,
#'                      endTime = 3.0,
#'                      toFile = FALSE)
#'
#' # Batch processing multiple files
#' files <- c("file1.wav", "file2.wav", "file3.mp3")
#' num_processed <- trk_phonet(files,
#'                             classes = c("nasal", "stop"),
#'                             outputDirectory = "results/",
#'                             toFile = TRUE)
#' cat(sprintf("Processed %d files\n", num_processed))
#'
#' # Integration with emuR
#' library(emuR)
#' # After running trk_phonet with toFile=TRUE:
#' # The .phn files can be added to emuR database
#' # and used for time-aligned phonological annotation
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso
#' \code{\link{install_phonet}} to install Phonet dependencies
#' \code{\link{phonet_available}} to check availability
#' \code{\link{phonet_info}} to get configuration information
#' \code{\link{lst_phonet}} for list-based output (not SSFF format)
#'
#' @export
trk_phonet <- function(listOfFiles,
                       classes = c("vocalic", "consonantal", "nasal", "stop"),
                       beginTime = 0.0,
                       endTime = 0.0,
                       explicitExt = "phn",
                       outputDirectory = NULL,
                       toFile = TRUE,
                       verbose = TRUE,
                       conda.env = NULL) {

  # Validate inputs
  if (length(listOfFiles) > 1 && !toFile) {
    stop("length(listOfFiles) is > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.",
         call. = FALSE)
  }

  # Check conda environment
  if (!is.null(conda.env) && !conda.env %in% reticulate::conda_list()$name) {
    stop("The conda environment ", conda.env, " does not exist.")
  }

  # Check Phonet availability
  if (!phonet_available()) {
    stop(
      "trk_phonet() requires Python with phonet, tensorflow, tf-keras, and speech processing libraries.\n\n",
      "Install with: install_phonet()\n",
      "Or manually: reticulate::py_install(c('tensorflow', 'tf-keras')); ",
      "reticulate::py_install('git+https://github.com/jcvasquezc/phonet.git')\n",
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

  # Initialize Phonet instance (reuse for all files for efficiency)
  if (verbose) {
    message("Initializing Phonet with classes: ", paste(classes, collapse = ", "))
  }

  phon <- phonet$Phonet(classes)

  # Vector of successfully processed files
  outListOfFiles <- c()

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
        channels = 1,  # Mono
        sample_rate = 16000,  # Phonet requires 16 kHz
        start_time = if (bt > 0) bt else NULL,
        total_time = if (et > bt && et > 0) (et - bt) else NULL
      )

      # Get audio info for metadata
      audio_info <- av::av_media_info(temp_wav)
      sample_rate <- as.numeric(audio_info$audio$sample_rate)

      # Run Phonet extraction
      df_result <- phon$get_phon_wav(temp_wav, "", plot_flag = FALSE)

      # Create AsspDataObj
      outDataObj <- list()

      # Determine track formats and names
      if ("all" %in% classes) {
        # Get all phonological classes from result
        all_cols <- names(df_result)
        phon_cols <- all_cols[!all_cols %in% c("time", "phoneme")]
        track_names <- phon_cols
      } else {
        track_names <- classes
      }

      # Set up track formats (all REAL32 for posteriors)
      track_formats <- rep("REAL32", length(track_names))

      attr(outDataObj, "trackFormats") <- track_formats
      attr(outDataObj, "sampleRate") <- 100.0  # Phonet uses 10ms frames = 100 Hz
      attr(outDataObj, "origFreq") <- as.numeric(sample_rate)
      attr(outDataObj, "startTime") <- as.numeric(if (bt > 0) bt else 0.0)
      attr(outDataObj, "startRecord") <- as.integer(1)
      attr(outDataObj, "endRecord") <- as.integer(nrow(df_result))
      class(outDataObj) <- "AsspDataObj"

      AsspFileFormat(outDataObj) <- "SSFF"
      AsspDataFormat(outDataObj) <- as.integer(2)  # binary

      # Add phonological posterior tracks
      for (track_name in track_names) {
        if (track_name %in% names(df_result)) {
          outDataObj <- addTrack(outDataObj,
                                track_name,
                                as.matrix(as.numeric(df_result[[track_name]])),
                                "REAL32")
        }
      }

      # Validate AsspDataObj
      if (!is.AsspDataObj(outDataObj)) {
        stop("The AsspDataObj created by trk_phonet is invalid.", call. = FALSE)
      }

      # Prepare output file path
      ssff_file <- sub("\\.[^.]+$", paste0(".", explicitExt), origSoundFile)
      if (!is.null(outputDirectory)) {
        if (!dir.exists(outputDirectory)) {
          dir.create(outputDirectory, recursive = TRUE)
        }
        ssff_file <- file.path(outputDirectory, basename(ssff_file))
      }

      attr(outDataObj, "filePath") <- as.character(ssff_file)

      # Write to file or accumulate
      if (toFile) {
        write.AsspDataObj(dobj = outDataObj, file = ssff_file)
        outListOfFiles <- c(outListOfFiles, TRUE)

        if (verbose) {
          message(sprintf("  → Written to: %s (%d tracks)", ssff_file, length(track_names)))
        }
      }

    }, error = function(e) {
      warning(sprintf("Error processing %s: %s", basename(origSoundFile), e$message))
      if (toFile) {
        outListOfFiles <- c(outListOfFiles, FALSE)
      }
    })
  }

  # Return based on toFile flag
  if (toFile) {
    return(sum(outListOfFiles))  # Count of successful files
  } else {
    return(outDataObj)
  }
}

# Set function attributes (for compatibility with existing infrastructure)
attr(trk_phonet, "ext") <- "phn"
attr(trk_phonet, "tracks") <- c("vocalic", "consonantal", "nasal", "stop")  # Default tracks
attr(trk_phonet, "outputType") <- "SSFF"

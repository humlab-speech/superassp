#' Optimized pitch analysis using pladdrr
#'
#' Computes fundamental frequency (F0) using pladdrr's native R bindings to Praat's
#' C library. Extracts pitch using cross-correlation (CC) and autocorrelation (AC)
#' methods. Audio is loaded using the av package and converted directly to pladdrr
#' Sound objects, eliminating temporary file creation. Supports all media formats
#' (WAV, MP3, MP4, etc.).
#'
#' @param listOfFiles Vector of file paths to audio files
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param time_step Time step between frames in seconds
#' @param window_length Analysis window length in seconds (not used by pladdrr direct API)
#' @param minimum_f0 Minimum pitch in Hz
#' @param maximum_f0 Maximum pitch in Hz
#' @param very_accurate Use very accurate pitch tracking
#' @param number_of_candidates Number of pitch candidates to consider
#' @param silence_threshold Silence threshold (0-1)
#' @param voicing_threshold Voicing threshold (0-1)
#' @param octave_cost Cost for octave jumps
#' @param octave_jump_cost Cost for octave jumps (same as octave_cost)
#' @param voiced_voiceless_cost Cost for voiced/voiceless transitions
#' @param only_correlation_methods If TRUE, only use CC and AC methods (default: TRUE)
#' @param minimum_filter_frequency Minimum filter frequency for SPINET (not available)
#' @param maximum_filter_frequency Maximum filter frequency for SPINET (not available)
#' @param number_of_filters Number of filters for SPINET (not available)
#' @param maximum_frequency_components Maximum frequency components for SHS (not available)
#' @param maximum_number_of_subharmonics Maximum number of subharmonics for SHS (not available)
#' @param compression_factor Compression factor for SHS (not available)
#' @param number_of_points_per_octave Number of points per octave for SHS (not available)
#' @param windowShape Window shape for time windowing
#' @param relativeWidth Relative width for windowing
#' @param toFile Write output to file (TRUE) or return object (FALSE)
#' @param explicitExt File extension for output files
#' @param outputDirectory Output directory path (NULL for same as input)
#' @param verbose Show progress messages (default: TRUE)
#' @param windowShift Deprecated. Use time_step instead.
#' @param corr.only Deprecated. Use only_correlation_methods instead.
#'
#' @return If toFile=TRUE, returns number of files successfully processed.
#'   If toFile=FALSE, returns an AsspDataObj with pitch tracks (cc and ac).
#'
#' @details
#' This function uses pladdrr (R bindings to Praat's C library) instead of
#' Python's parselmouth. Advantages:
#' \itemize{
#'   \item Pure R/C implementation (no Python dependency)
#'   \item Native R6 object-oriented interface
#'   \item Direct C library access for performance
#'   \item No numpy conversion overhead
#' }
#'
#' The function extracts pitch using two methods:
#' \itemize{
#'   \item **CC (Cross-Correlation)**: Most robust for speech
#'   \item **AC (Autocorrelation)**: Fast, good for clean signals
#' }
#'
#' **Note**: SPINET and SHS methods are not available in pladdrr's direct API.
#' Parameters related to these methods are accepted for backwards compatibility
#' but are ignored.
#'
#' @examples
#' \dontrun{
#' # Analyze pitch for a single file
#' pitch <- trk_pitchp("speech.wav", toFile = FALSE)
#'
#' # Batch process multiple files
#' files <- c("speech1.wav", "speech2.wav")
#' trk_pitchp(files, toFile = TRUE, outputDirectory = "output/")
#'
#' # With custom F0 range
#' pitch <- trk_pitchp("speech.wav",
#'                     minimum_f0 = 100,
#'                     maximum_f0 = 300,
#'                     toFile = FALSE)
#' }
#'
#' @export
trk_pitchp <- function(listOfFiles,
                       beginTime = 0.0,
                       endTime = 0.0,
                       time_step = 0.005,
                       window_length = 0.040,
                       minimum_f0 = 75.0,
                       maximum_f0 = 600.0,
                       very_accurate = TRUE,
                       number_of_candidates = 15,
                       silence_threshold = 0.03,
                       voicing_threshold = 0.45,
                       octave_cost = 0.01,
                       octave_jump_cost = 0.35,
                       voiced_voiceless_cost = 0.14,
                       only_correlation_methods = TRUE,
                       minimum_filter_frequency = 70.0,
                       maximum_filter_frequency = 5000.0,
                       number_of_filters = 250,
                       maximum_frequency_components = 1250.0,
                       maximum_number_of_subharmonics = 15,
                       compression_factor = 0.84,
                       number_of_points_per_octave = 48,
                       windowShape = "Gaussian1",
                       relativeWidth = 1.0,
                       toFile = TRUE,
                       explicitExt = "pit",
                       outputDirectory = NULL,
                       verbose = TRUE,
                       # Backwards compatibility parameters
                       windowShift = NULL,
                       corr.only = NULL) {

  # Check pladdrr availability
  if (!pladdrr_available()) {
    cli::cli_abort(c(
      "x" = "pladdrr package not available",
      "i" = "Install with: install_pladdrr()"
    ))
  }

  # Handle backwards compatibility: windowShift -> time_step
  if (!is.null(windowShift)) {
    warning("Parameter 'windowShift' is deprecated. Please use 'time_step' instead.", call. = FALSE)
    time_step <- windowShift / 1000  # Convert from ms to seconds
  }

  # Handle backwards compatibility: corr.only -> only_correlation_methods
  if (!is.null(corr.only)) {
    warning("Parameter 'corr.only' is deprecated. Please use 'only_correlation_methods' instead.", call. = FALSE)
    only_correlation_methods <- corr.only
  }

  # Warn about SPINET/SHS if requested
  if (!only_correlation_methods) {
    warning("SPINET and SHS methods are not available in pladdrr direct API. Only CC and AC will be computed.", call. = FALSE)
  }

  # Check if multiple files with toFile=FALSE
  if(length(listOfFiles) > 1 & !toFile) {
    stop("length(listOfFiles) is > 1 and toFile=FALSE! toFile=FALSE only permitted for single files.")
  }

  # Create data frame for file processing
  tryCatch({
    fileBeginEnd <- data.frame(
      listOfFiles = listOfFiles,
      beginTime = beginTime,
      endTime = endTime
    )
  }, error = function(e) {
    stop("The beginTime and endTime must either be a single value or the same length as listOfFiles")
  })

  # Check that all files exist
  filesEx <- file.exists(listOfFiles)
  if(!all(filesEx)) {
    filesNotExist <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ", paste(filesNotExist, collapse = ", "))
  }

  outListOfFiles <- c()
  n_files <- nrow(fileBeginEnd)

  # Progress bar setup
  if (verbose && n_files > 1) {
    pb <- utils::txtProgressBar(min = 0, max = n_files, style = 3)
  }

  # Process each file
  for(i in 1:n_files) {
    origSoundFile <- normalizePath(fileBeginEnd[i, "listOfFiles"], mustWork = TRUE)

    bt <- fileBeginEnd[i, "beginTime"]
    et <- fileBeginEnd[i, "endTime"]

    # Load audio and create pladdrr Sound object (file-based)
    sound <- av_load_for_pladdrr(
      file_path = origSoundFile,
      start_time = bt,
      end_time = et,
      window_type = windowShape,
      relative_width = relativeWidth
    )

    # sound is already windowed by av_load_for_pladdrr
    analysis_sound <- sound

    # Get duration
    duration <- sound$.cpp$duration


    # Extract pitch using CC (cross-correlation) - Direct API
    pitch_cc_ptr <- pladdrr::to_pitch_cc_direct(
      analysis_sound,
      time_step = time_step,
      pitch_floor = minimum_f0,
      pitch_ceiling = maximum_f0,
      max_candidates = as.integer(number_of_candidates),
      very_accurate = very_accurate,
      silence_threshold = silence_threshold,
      voicing_threshold = voicing_threshold,
      octave_cost = octave_cost,
      octave_jump_cost = octave_jump_cost,
      voiced_unvoiced_cost = voiced_voiceless_cost
    )
    pitch_cc <- pladdrr::Pitch(.xptr = pitch_cc_ptr)

    # Extract pitch using AC (autocorrelation) - Direct API
    pitch_ac_ptr <- pladdrr::to_pitch_ac_direct(
      analysis_sound,
      time_step = time_step,
      pitch_floor = minimum_f0,
      pitch_ceiling = maximum_f0,
      max_candidates = as.integer(number_of_candidates),
      very_accurate = very_accurate,
      silence_threshold = silence_threshold,
      voicing_threshold = voicing_threshold,
      octave_cost = octave_cost,
      octave_jump_cost = octave_jump_cost,
      voiced_unvoiced_cost = voiced_voiceless_cost
    )
    pitch_ac <- pladdrr::Pitch(.xptr = pitch_ac_ptr)

    # Get frame data
    cc_df <- pitch_cc$as_data_frame()
    ac_df <- pitch_ac$as_data_frame()

    num_frames_cc <- nrow(cc_df)
    num_frames_ac <- nrow(ac_df)

    # Use CC pitch times as reference
    times <- cc_df[[1]]
    cc_values <- cc_df[[2]]

    # Get AC values by frame index (matching plabench behavior)
    num_frames <- num_frames_cc
    ac_values <- rep(NA_real_, num_frames)

    # For indices within AC's range, get the value at that frame index
    for (j in seq_len(min(num_frames, num_frames_ac))) {
      ac_values[j] <- ac_df[[2]][j]
    }

    # Replace 0 and undefined values with NA for unvoiced frames
    cc_values[cc_values == 0 | !is.finite(cc_values)] <- NA
    ac_values[ac_values == 0 | !is.finite(ac_values)] <- NA

    # Create AsspDataObj
    outDataObj <- list()

    # Calculate sample rate from pitch object
    actual_time_step <- pitch_cc$get_time_step()
    sampleRate <- 1 / actual_time_step

    attr(outDataObj, "sampleRate") <- sampleRate
    attr(outDataObj, "origFreq") <- as.numeric(16000)

    startTime <- times[1]
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(num_frames)

    class(outDataObj) <- "AsspDataObj"
    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2)
    attr(outDataObj, "trackFormats") <- character(0)

    # Add CC track
    cc_data <- cc_values
    cc_data[is.na(cc_data)] <- 0
    outDataObj <- wrassp::addTrack(outDataObj, "pitch_cc",
                                   as.matrix(cc_data), "REAL32")
    attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")

    # Add AC track
    ac_data <- ac_values
    ac_data[is.na(ac_data)] <- 0
    outDataObj <- wrassp::addTrack(outDataObj, "pitch_ac",
                                   as.matrix(ac_data), "REAL32")
    attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")

    assertthat::assert_that(wrassp::is.AsspDataObj(outDataObj),
                            msg = "The AsspDataObj created by pitch analysis is invalid.")

    # Determine output file path
    ssff_file <- sub("\\.[^.]*$", paste0(".", explicitExt), origSoundFile)
    if(!is.null(outputDirectory)) {
      ssff_file <- file.path(outputDirectory, basename(ssff_file))
    }

    attr(outDataObj, "filePath") <- as.character(ssff_file)

    if(toFile) {
      wrassp::write.AsspDataObj(dobj = outDataObj, file = ssff_file)
      outListOfFiles <- c(outListOfFiles, TRUE)
    }

    # Update progress bar
    if (verbose && n_files > 1) {
      utils::setTxtProgressBar(pb, i)
    }
  }

  # Close progress bar
  if (verbose && n_files > 1) {
    close(pb)
  }

  if(toFile) {
    if (verbose) {
      cli::cli_alert_success("Processed {length(outListOfFiles)} file{?s}")
    }
    return(length(outListOfFiles))
  } else {
    return(outDataObj)
  }
}

attr(trk_pitchp, "ext") <- c("pit")
attr(trk_pitchp, "tracks") <- c("pitch_cc", "pitch_ac")
attr(trk_pitchp, "outputType") <- c("SSFF")
attr(trk_pitchp, "nativeFiletypes") <- c("wav")


# ── Private helper: run one pladdrr pitch method on a single Sound object ─────

.pitchp_to_assp <- function(sound, method, time_step, minimum_f0, maximum_f0,
                             number_of_candidates, very_accurate,
                             silence_threshold, voicing_threshold,
                             octave_cost, octave_jump_cost,
                             voiced_voiceless_cost) {

  fn <- if (method == "cc") pladdrr::to_pitch_cc_direct else pladdrr::to_pitch_ac_direct
  ptr <- fn(
    sound,
    time_step             = time_step,
    pitch_floor           = minimum_f0,
    pitch_ceiling         = maximum_f0,
    max_candidates        = as.integer(number_of_candidates),
    very_accurate         = very_accurate,
    silence_threshold     = silence_threshold,
    voicing_threshold     = voicing_threshold,
    octave_cost           = octave_cost,
    octave_jump_cost      = octave_jump_cost,
    voiced_unvoiced_cost  = voiced_voiceless_cost
  )
  pitch_obj <- pladdrr::Pitch(.xptr = ptr)
  df        <- pitch_obj$as_data_frame()        # (time, frequency)
  times     <- df[[1]]
  values    <- df[[2]]
  values[values == 0 | !is.finite(values)] <- NA

  list(
    times      = times,
    values     = values,
    time_step  = pitch_obj$get_time_step(),
    n_frames   = nrow(df)
  )
}


# ── trk_pitchp_cc ─────────────────────────────────────────────────────────────

#' Pitch tracking via Praat cross-correlation (CC) method
#'
#' Computes fundamental frequency (F0) using Praat's cross-correlation pitch
#' algorithm via pladdrr. Returns a single \code{F0} track.
#'
#' @inheritParams trk_pitchp
#' @param explicitExt File extension for output files. Default \code{"pcc"}.
#'
#' @return If \code{toFile = FALSE}: an AsspDataObj with one \code{F0} track
#'   (REAL32, Hz; 0 = unvoiced). If \code{toFile = TRUE}: number of files
#'   successfully written.
#'
#' @seealso \code{\link{trk_pitchp_ac}}, \code{\link{trk_pitchp}}
#'
#' @export
trk_pitchp_cc <- function(listOfFiles,
                           beginTime              = 0.0,
                           endTime                = 0.0,
                           time_step              = 0.005,
                           minimum_f0             = 75.0,
                           maximum_f0             = 600.0,
                           very_accurate          = TRUE,
                           number_of_candidates   = 15,
                           silence_threshold      = 0.03,
                           voicing_threshold      = 0.45,
                           octave_cost            = 0.01,
                           octave_jump_cost       = 0.35,
                           voiced_voiceless_cost  = 0.14,
                           windowShape            = "Gaussian1",
                           relativeWidth          = 1.0,
                           toFile                 = TRUE,
                           explicitExt            = "pcc",
                           outputDirectory        = NULL,
                           verbose                = TRUE) {

  if (!pladdrr_available()) {
    cli::cli_abort(c("x" = "pladdrr package not available",
                     "i" = "Install with: install_pladdrr()"))
  }
  if (length(listOfFiles) > 1 && !toFile) {
    stop("toFile=FALSE only permitted for single files.", call. = FALSE)
  }

  tryCatch(
    fileBeginEnd <- data.frame(listOfFiles = listOfFiles,
                               beginTime = beginTime, endTime = endTime),
    error = function(e) stop("beginTime/endTime must be length 1 or same length as listOfFiles",
                             call. = FALSE)
  )

  filesEx <- file.exists(listOfFiles)
  if (!all(filesEx))
    stop("Unable to find: ", paste(listOfFiles[!filesEx], collapse = ", "), call. = FALSE)

  outListOfFiles <- c()
  n_files <- nrow(fileBeginEnd)
  if (verbose && n_files > 1) pb <- utils::txtProgressBar(min = 0, max = n_files, style = 3)

  for (i in seq_len(n_files)) {
    origSoundFile <- normalizePath(fileBeginEnd[i, "listOfFiles"], mustWork = TRUE)
    bt <- fileBeginEnd[i, "beginTime"]
    et <- fileBeginEnd[i, "endTime"]

    sound <- av_load_for_pladdrr(origSoundFile, start_time = bt, end_time = et,
                                  window_type = windowShape, relative_width = relativeWidth)

    res <- .pitchp_to_assp(sound, "cc", time_step, minimum_f0, maximum_f0,
                            number_of_candidates, very_accurate, silence_threshold,
                            voicing_threshold, octave_cost, octave_jump_cost,
                            voiced_voiceless_cost)

    outDataObj <- list()
    attr(outDataObj, "sampleRate")   <- 1 / res$time_step
    attr(outDataObj, "origFreq")     <- as.numeric(sound$.cpp$sampling_frequency)
    attr(outDataObj, "startTime")    <- as.numeric(res$times[1])
    attr(outDataObj, "startRecord")  <- 1L
    attr(outDataObj, "endRecord")    <- as.integer(res$n_frames)
    attr(outDataObj, "trackFormats") <- character(0)
    class(outDataObj) <- "AsspDataObj"
    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- 2L

    f0_data <- res$values
    f0_data[is.na(f0_data)] <- 0
    outDataObj <- addTrack(outDataObj, "F0", as.matrix(f0_data), "REAL32")
    attr(outDataObj, "trackFormats") <- "REAL32"

    ssff_file <- sub("\\.[^.]*$", paste0(".", explicitExt), origSoundFile)
    if (!is.null(outputDirectory))
      ssff_file <- file.path(outputDirectory, basename(ssff_file))
    attr(outDataObj, "filePath") <- as.character(ssff_file)

    if (toFile) {
      write.AsspDataObj(dobj = outDataObj, file = ssff_file)
      outListOfFiles <- c(outListOfFiles, TRUE)
    }

    if (verbose && n_files > 1) utils::setTxtProgressBar(pb, i)
  }

  if (verbose && n_files > 1) close(pb)

  if (toFile) {
    if (verbose) cli::cli_alert_success("Processed {length(outListOfFiles)} file{?s}")
    return(length(outListOfFiles))
  } else {
    return(outDataObj)
  }
}

attr(trk_pitchp_cc, "ext")             <- "pcc"
attr(trk_pitchp_cc, "tracks")          <- "F0"
attr(trk_pitchp_cc, "outputType")      <- "SSFF"
attr(trk_pitchp_cc, "nativeFiletypes") <- "wav"


# ── trk_pitchp_ac ─────────────────────────────────────────────────────────────

#' Pitch tracking via Praat autocorrelation (AC) method
#'
#' Computes fundamental frequency (F0) using Praat's autocorrelation pitch
#' algorithm via pladdrr. Returns a single \code{F0} track.
#'
#' @inheritParams trk_pitchp
#' @param explicitExt File extension for output files. Default \code{"pac"}.
#'
#' @return If \code{toFile = FALSE}: an AsspDataObj with one \code{F0} track
#'   (REAL32, Hz; 0 = unvoiced). If \code{toFile = TRUE}: number of files
#'   successfully written.
#'
#' @seealso \code{\link{trk_pitchp_cc}}, \code{\link{trk_pitchp}}
#'
#' @export
trk_pitchp_ac <- function(listOfFiles,
                           beginTime              = 0.0,
                           endTime                = 0.0,
                           time_step              = 0.005,
                           minimum_f0             = 75.0,
                           maximum_f0             = 600.0,
                           very_accurate          = TRUE,
                           number_of_candidates   = 15,
                           silence_threshold      = 0.03,
                           voicing_threshold      = 0.45,
                           octave_cost            = 0.01,
                           octave_jump_cost       = 0.35,
                           voiced_voiceless_cost  = 0.14,
                           windowShape            = "Gaussian1",
                           relativeWidth          = 1.0,
                           toFile                 = TRUE,
                           explicitExt            = "pac",
                           outputDirectory        = NULL,
                           verbose                = TRUE) {

  if (!pladdrr_available()) {
    cli::cli_abort(c("x" = "pladdrr package not available",
                     "i" = "Install with: install_pladdrr()"))
  }
  if (length(listOfFiles) > 1 && !toFile) {
    stop("toFile=FALSE only permitted for single files.", call. = FALSE)
  }

  tryCatch(
    fileBeginEnd <- data.frame(listOfFiles = listOfFiles,
                               beginTime = beginTime, endTime = endTime),
    error = function(e) stop("beginTime/endTime must be length 1 or same length as listOfFiles",
                             call. = FALSE)
  )

  filesEx <- file.exists(listOfFiles)
  if (!all(filesEx))
    stop("Unable to find: ", paste(listOfFiles[!filesEx], collapse = ", "), call. = FALSE)

  outListOfFiles <- c()
  n_files <- nrow(fileBeginEnd)
  if (verbose && n_files > 1) pb <- utils::txtProgressBar(min = 0, max = n_files, style = 3)

  for (i in seq_len(n_files)) {
    origSoundFile <- normalizePath(fileBeginEnd[i, "listOfFiles"], mustWork = TRUE)
    bt <- fileBeginEnd[i, "beginTime"]
    et <- fileBeginEnd[i, "endTime"]

    sound <- av_load_for_pladdrr(origSoundFile, start_time = bt, end_time = et,
                                  window_type = windowShape, relative_width = relativeWidth)

    res <- .pitchp_to_assp(sound, "ac", time_step, minimum_f0, maximum_f0,
                            number_of_candidates, very_accurate, silence_threshold,
                            voicing_threshold, octave_cost, octave_jump_cost,
                            voiced_voiceless_cost)

    outDataObj <- list()
    attr(outDataObj, "sampleRate")   <- 1 / res$time_step
    attr(outDataObj, "origFreq")     <- as.numeric(sound$.cpp$sampling_frequency)
    attr(outDataObj, "startTime")    <- as.numeric(res$times[1])
    attr(outDataObj, "startRecord")  <- 1L
    attr(outDataObj, "endRecord")    <- as.integer(res$n_frames)
    attr(outDataObj, "trackFormats") <- character(0)
    class(outDataObj) <- "AsspDataObj"
    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- 2L

    f0_data <- res$values
    f0_data[is.na(f0_data)] <- 0
    outDataObj <- addTrack(outDataObj, "F0", as.matrix(f0_data), "REAL32")
    attr(outDataObj, "trackFormats") <- "REAL32"

    ssff_file <- sub("\\.[^.]*$", paste0(".", explicitExt), origSoundFile)
    if (!is.null(outputDirectory))
      ssff_file <- file.path(outputDirectory, basename(ssff_file))
    attr(outDataObj, "filePath") <- as.character(ssff_file)

    if (toFile) {
      write.AsspDataObj(dobj = outDataObj, file = ssff_file)
      outListOfFiles <- c(outListOfFiles, TRUE)
    }

    if (verbose && n_files > 1) utils::setTxtProgressBar(pb, i)
  }

  if (verbose && n_files > 1) close(pb)

  if (toFile) {
    if (verbose) cli::cli_alert_success("Processed {length(outListOfFiles)} file{?s}")
    return(length(outListOfFiles))
  } else {
    return(outDataObj)
  }
}

attr(trk_pitchp_ac, "ext")             <- "pac"
attr(trk_pitchp_ac, "tracks")          <- "F0"
attr(trk_pitchp_ac, "outputType")      <- "SSFF"
attr(trk_pitchp_ac, "nativeFiletypes") <- "wav"

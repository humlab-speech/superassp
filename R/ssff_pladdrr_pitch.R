
# ── Private helper: run one pladdrr pitch method on a single Sound object ─────

.pitchp_to_assp <- function(sound, method, time_step, minimum_f0, maximum_f0,
                             number_of_candidates, very_accurate,
                             silence_threshold, voicing_threshold,
                             octave_cost, octave_jump_cost,
                             voiced_voiceless_cost) {

  on.exit(invisible(sound))
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


# ── trk_pitch_cc ─────────────────────────────────────────────────────────────

#' Pitch tracking via Praat cross-correlation (CC) method
#'
#' Computes fundamental frequency (F0) using Praat's cross-correlation pitch
#' algorithm via pladdrr. Returns a single \code{F0} track.
#'
#' @param listOfFiles Vector of file paths to audio files
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param time_step Time step between frames in seconds
#' @param minimum_f0 Minimum pitch in Hz
#' @param maximum_f0 Maximum pitch in Hz
#' @param very_accurate Use very accurate pitch tracking
#' @param number_of_candidates Number of pitch candidates to consider
#' @param silence_threshold Silence threshold (0-1)
#' @param voicing_threshold Voicing threshold (0-1)
#' @param octave_cost Cost for octave jumps
#' @param octave_jump_cost Cost for octave jumps
#' @param voiced_voiceless_cost Cost for voiced/voiceless transitions
#' @param windowShape Window shape for time windowing
#' @param relativeWidth Relative width for windowing
#' @param toFile Write output to file (TRUE) or return object (FALSE)
#' @param explicitExt File extension for output files. Default \code{"pcc"}.
#'
#' @return If \code{toFile = FALSE}: an AsspDataObj with one \code{F0} track
#'   (REAL32, Hz; 0 = unvoiced). If \code{toFile = TRUE}: number of files
#'   successfully written.
#'
#' @seealso \code{\link{trk_pitch_ac}}, \code{\link{trk_pitch_shs}}, \code{\link{trk_pitch_spinet}}
#'
#' @export
trk_pitch_cc <- function(listOfFiles,
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

attr(trk_pitch_cc, "ext")             <- "pcc"
attr(trk_pitch_cc, "tracks")          <- "F0"
attr(trk_pitch_cc, "outputType")      <- "SSFF"
attr(trk_pitch_cc, "nativeFiletypes") <- "wav"


# ── trk_pitch_ac ─────────────────────────────────────────────────────────────

#' Pitch tracking via Praat autocorrelation (AC) method
#'
#' Computes fundamental frequency (F0) using Praat's autocorrelation pitch
#' algorithm via pladdrr. Returns a single \code{F0} track.
#'
#' @inheritParams trk_pitch_cc
#' @param explicitExt File extension for output files. Default \code{"pac"}.
#'
#' @return If \code{toFile = FALSE}: an AsspDataObj with one \code{F0} track
#'   (REAL32, Hz; 0 = unvoiced). If \code{toFile = TRUE}: number of files
#'   successfully written.
#'
#' @seealso \code{\link{trk_pitch_cc}}
#'
#' @export
trk_pitch_ac <- function(listOfFiles,
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

attr(trk_pitch_ac, "ext")             <- "pac"
attr(trk_pitch_ac, "tracks")          <- "F0"
attr(trk_pitch_ac, "outputType")      <- "SSFF"
attr(trk_pitch_ac, "nativeFiletypes") <- "wav"


# ── Private helpers: SHS and SPINET ──────────────────────────────────────────

.pitchp_shs_to_assp <- function(sound, time_step, minimum_f0, maximum_f0,
                                 maximum_frequency_components,
                                 maximum_number_of_subharmonics,
                                 number_of_candidates,
                                 compression_factor,
                                 number_of_points_per_octave) {
  on.exit(invisible(sound))
  ptr <- pladdrr::to_pitch_shs_direct(
    sound,
    time_step           = time_step,
    pitch_floor         = minimum_f0,
    max_frequency       = maximum_frequency_components,
    pitch_ceiling       = maximum_f0,
    max_subharmonics    = as.integer(maximum_number_of_subharmonics),
    max_candidates      = as.integer(number_of_candidates),
    compression_factor  = compression_factor,
    n_points_per_octave = as.integer(number_of_points_per_octave)
  )
  pitch_obj <- pladdrr::Pitch(.xptr = ptr)
  df     <- pitch_obj$as_data_frame()
  values <- df[[2]]
  values[values == 0 | !is.finite(values)] <- NA
  list(times = df[[1]], values = values,
       time_step = pitch_obj$get_time_step(), n_frames = nrow(df))
}

.pitchp_spinet_to_assp <- function(sound, time_step, window_length,
                                    minimum_filter_frequency,
                                    maximum_filter_frequency,
                                    number_of_filters,
                                    maximum_f0,
                                    number_of_candidates) {
  on.exit(invisible(sound))
  ptr <- pladdrr::to_pitch_spinet_direct(
    sound,
    time_step       = time_step,
    window_duration = window_length,
    min_frequency   = minimum_filter_frequency,
    max_frequency   = maximum_filter_frequency,
    n_filters       = as.integer(number_of_filters),
    pitch_ceiling   = maximum_f0,
    max_candidates  = as.integer(number_of_candidates)
  )
  pitch_obj <- pladdrr::Pitch(.xptr = ptr)
  df     <- pitch_obj$as_data_frame()
  values <- df[[2]]
  values[values == 0 | !is.finite(values)] <- NA
  list(times = df[[1]], values = values,
       time_step = pitch_obj$get_time_step(), n_frames = nrow(df))
}


# ── trk_pitch_shs ──────────────────────────────────────────────────────────

#' Pitch tracking via Praat subharmonic summation (SHS) method
#'
#' Computes fundamental frequency (F0) using Praat's subharmonic summation
#' algorithm via pladdrr. Returns a single \code{F0} track.
#'
#' @param listOfFiles Vector of file paths to audio files
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param time_step Time step between frames in seconds
#' @param minimum_f0 Minimum pitch in Hz
#' @param maximum_f0 Maximum pitch (ceiling) in Hz
#' @param maximum_frequency_components Maximum frequency for spectral analysis in Hz
#' @param maximum_number_of_subharmonics Maximum number of subharmonics to consider
#' @param number_of_candidates Number of pitch candidates
#' @param compression_factor Spectral compression factor
#' @param number_of_points_per_octave Spectral resolution in points per octave
#' @param windowShape Window shape for time windowing
#' @param relativeWidth Relative width for windowing
#' @param toFile Write output to file (TRUE) or return object (FALSE)
#' @param explicitExt File extension for output files. Default \code{"psh"}.
#' @param outputDirectory Output directory path (NULL for same as input)
#' @param verbose Show progress messages (default: TRUE)
#'
#' @return If \code{toFile = FALSE}: an AsspDataObj with one \code{F0} track
#'   (REAL32, Hz; 0 = unvoiced). If \code{toFile = TRUE}: number of files
#'   successfully written.
#'
#' @seealso \code{\link{trk_pitch_cc}}, \code{\link{trk_pitch_ac}},
#'   \code{\link{trk_pitch_spinet}}
#'
#' @export
trk_pitch_shs <- function(listOfFiles,
                            beginTime                       = 0.0,
                            endTime                         = 0.0,
                            time_step                       = 0.01,
                            minimum_f0                      = 50.0,
                            maximum_f0                      = 500.0,
                            maximum_frequency_components    = 1250.0,
                            maximum_number_of_subharmonics  = 15,
                            number_of_candidates            = 15,
                            compression_factor              = 0.84,
                            number_of_points_per_octave     = 48,
                            windowShape                     = "Gaussian1",
                            relativeWidth                   = 1.0,
                            toFile                          = TRUE,
                            explicitExt                     = "psh",
                            outputDirectory                 = NULL,
                            verbose                         = TRUE) {

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

    res <- .pitchp_shs_to_assp(sound, time_step, minimum_f0, maximum_f0,
                                maximum_frequency_components,
                                maximum_number_of_subharmonics,
                                number_of_candidates,
                                compression_factor,
                                number_of_points_per_octave)

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

attr(trk_pitch_shs, "ext")             <- "psh"
attr(trk_pitch_shs, "tracks")          <- "F0"
attr(trk_pitch_shs, "outputType")      <- "SSFF"
attr(trk_pitch_shs, "nativeFiletypes") <- "wav"


# ── trk_pitch_spinet ───────────────────────────────────────────────────────

#' Pitch tracking via Praat SPINET method
#'
#' Computes fundamental frequency (F0) using Praat's SPINET (SPectral
#' INtegration and Evaluation of Temporal patterns) algorithm via pladdrr.
#' Returns a single \code{F0} track.
#'
#' @param listOfFiles Vector of file paths to audio files
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param time_step Time step between frames in seconds
#' @param window_length Analysis window duration in seconds
#' @param minimum_filter_frequency Minimum filter frequency in Hz
#' @param maximum_filter_frequency Maximum filter frequency in Hz
#' @param number_of_filters Number of gammatone filters
#' @param maximum_f0 Maximum pitch (ceiling) in Hz
#' @param number_of_candidates Number of pitch candidates
#' @param windowShape Window shape for time windowing
#' @param relativeWidth Relative width for windowing
#' @param toFile Write output to file (TRUE) or return object (FALSE)
#' @param explicitExt File extension for output files. Default \code{"psp"}.
#' @param outputDirectory Output directory path (NULL for same as input)
#' @param verbose Show progress messages (default: TRUE)
#'
#' @return If \code{toFile = FALSE}: an AsspDataObj with one \code{F0} track
#'   (REAL32, Hz; 0 = unvoiced). If \code{toFile = TRUE}: number of files
#'   successfully written.
#'
#' @seealso \code{\link{trk_pitch_cc}}, \code{\link{trk_pitch_ac}},
#'   \code{\link{trk_pitch_shs}}
#'
#' @export
trk_pitch_spinet <- function(listOfFiles,
                               beginTime                = 0.0,
                               endTime                  = 0.0,
                               time_step                = 0.005,
                               window_length            = 0.04,
                               minimum_filter_frequency = 70.0,
                               maximum_filter_frequency = 5000.0,
                               number_of_filters        = 250,
                               maximum_f0               = 500.0,
                               number_of_candidates     = 15,
                               windowShape              = "Gaussian1",
                               relativeWidth            = 1.0,
                               toFile                   = TRUE,
                               explicitExt              = "psp",
                               outputDirectory          = NULL,
                               verbose                  = TRUE) {

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

    res <- .pitchp_spinet_to_assp(sound, time_step, window_length,
                                   minimum_filter_frequency,
                                   maximum_filter_frequency,
                                   number_of_filters,
                                   maximum_f0,
                                   number_of_candidates)

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

attr(trk_pitch_spinet, "ext")             <- "psp"
attr(trk_pitch_spinet, "tracks")          <- "F0"
attr(trk_pitch_spinet, "outputType")      <- "SSFF"
attr(trk_pitch_spinet, "nativeFiletypes") <- "wav"

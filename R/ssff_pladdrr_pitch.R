
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

#' Pitch tracking via Praat's cross-correlation method
#'
#' Tracks fundamental frequency (F0) using Praat's cross-correlation (CC) pitch
#' algorithm via pladdrr. The CC method correlates short-term waveforms across
#' time and is more robust to noise than autocorrelation at the cost of slightly
#' higher compute time.
#'
#' @param listOfFiles Character vector of audio file paths. Any format supported by
#'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
#' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
#' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
#' @param time_step Numeric. Frame shift in seconds; sets output frame rate
#'   (1 / time_step Hz). Default 0.005 s (200 Hz).
#' @param minimum_f0 Numeric. Lower F0 bound in Hz. Default 75 Hz.
#' @param maximum_f0 Numeric. Upper F0 bound (ceiling) in Hz. Default 600 Hz.
#' @param very_accurate Logical. Use slower, higher-accuracy candidate search.
#'   Default \code{TRUE}.
#' @param number_of_candidates Integer. Maximum pitch candidates per frame.
#'   Default 15.
#' @param silence_threshold Numeric. Frames with amplitude below this fraction of
#'   the global maximum are treated as silent (0–1). Default 0.03.
#' @param voicing_threshold Numeric. Minimum strength for a frame to be voiced (0–1).
#'   Default 0.45.
#' @param octave_cost Numeric. Penalty per octave above \code{minimum_f0} to
#'   discourage high-frequency candidates. Default 0.01.
#' @param octave_jump_cost Numeric. Penalty for octave jumps between adjacent frames.
#'   Default 0.35.
#' @param voiced_voiceless_cost Numeric. Penalty for voiced/unvoiced transitions.
#'   Default 0.14.
#' @param windowShape Character. Window shape applied to audio before loading.
#'   Default \code{"Gaussian1"}.
#' @param relativeWidth Numeric. Relative width of the window. Default 1.0.
#' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
#'   count written (invisibly). If \code{FALSE}, return an \code{AsspDataObj}.
#'   Default \code{TRUE}.
#' @param explicitExt Character. Output file extension. Default \code{"pcc"}.
#' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
#'   writes alongside the input file.
#' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
#'
#' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
#'   \describe{
#'     \item{\code{F0}}{REAL32, Hz, n_frames x 1. Fundamental frequency.
#'       0 encodes unvoiced frames.}
#'   }
#'   Frame rate: \code{1 / time_step} Hz (default 200 Hz).
#'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
#'
#' @seealso \code{\link{trk_pitch_ac}}, \code{\link{trk_pitch_shs}}, \code{\link{trk_pitch_spinet}}
#'
#' @examples
#' \dontrun{
#' trk_pitch_cc(
#'   system.file("samples", "sustained", "a1.wav", package = "superassp"),
#'   toFile = FALSE
#' )
#' }
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

#' Pitch tracking via Praat's autocorrelation method
#'
#' Tracks fundamental frequency (F0) using Praat's autocorrelation (AC) pitch
#' algorithm via pladdrr. The AC method is Praat's standard pitch tracker and
#' is well-suited to clean speech; prefer \code{trk_pitch_cc} for noisier signals.
#'
#' @inheritParams trk_pitch_cc
#' @param explicitExt Character. Output file extension. Default \code{"pac"}.
#'
#' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
#'   \describe{
#'     \item{\code{F0}}{REAL32, Hz, n_frames x 1. Fundamental frequency.
#'       0 encodes unvoiced frames.}
#'   }
#'   Frame rate: \code{1 / time_step} Hz (default 200 Hz).
#'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
#'
#' @seealso \code{\link{trk_pitch_cc}}, \code{\link{trk_pitch_shs}},
#'   \code{\link{trk_pitch_spinet}}
#'
#' @examples
#' \dontrun{
#' trk_pitch_ac(
#'   system.file("samples", "sustained", "a1.wav", package = "superassp"),
#'   toFile = FALSE
#' )
#' }
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

#' Pitch tracking via Praat's subharmonic summation (SHS) method
#'
#' Tracks fundamental frequency (F0) using Praat's subharmonic summation (SHS)
#' algorithm via pladdrr. SHS works in the frequency domain by summing
#' subharmonically compressed spectral representations and is robust to
#' spectral-domain interference; prefer it over AC/CC for speech with strong
#' vocal fry or very low F0.
#'
#' @param listOfFiles Character vector of audio file paths. Any format supported by
#'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
#' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
#' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
#' @param time_step Numeric. Frame shift in seconds; sets output frame rate
#'   (1 / time_step Hz). Default 0.01 s (100 Hz).
#' @param minimum_f0 Numeric. Lower F0 bound in Hz. Default 50 Hz.
#' @param maximum_f0 Numeric. Upper F0 bound (ceiling) in Hz. Default 500 Hz.
#' @param maximum_frequency_components Numeric. Highest frequency included in the
#'   spectral analysis in Hz. Default 1250 Hz.
#' @param maximum_number_of_subharmonics Integer. Maximum number of subharmonic
#'   levels to sum. Default 15.
#' @param number_of_candidates Integer. Maximum pitch candidates per frame.
#'   Default 15.
#' @param compression_factor Numeric. Spectral compression factor applied during
#'   subharmonic summation. Default 0.84.
#' @param number_of_points_per_octave Integer. Spectral resolution in the
#'   log-frequency representation (points per octave). Default 48.
#' @param windowShape Character. Window shape applied to audio before loading.
#'   Default \code{"Gaussian1"}.
#' @param relativeWidth Numeric. Relative width of the window. Default 1.0.
#' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
#'   count written (invisibly). If \code{FALSE}, return an \code{AsspDataObj}.
#'   Default \code{TRUE}.
#' @param explicitExt Character. Output file extension. Default \code{"psh"}.
#' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
#'   writes alongside the input file.
#' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
#'
#' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
#'   \describe{
#'     \item{\code{F0}}{REAL32, Hz, n_frames x 1. Fundamental frequency.
#'       0 encodes unvoiced frames.}
#'   }
#'   Frame rate: \code{1 / time_step} Hz (default 100 Hz).
#'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
#'
#' @seealso \code{\link{trk_pitch_cc}}, \code{\link{trk_pitch_ac}},
#'   \code{\link{trk_pitch_spinet}}
#'
#' @examples
#' \dontrun{
#' trk_pitch_shs(
#'   system.file("samples", "sustained", "a1.wav", package = "superassp"),
#'   toFile = FALSE
#' )
#' }
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

#' Pitch tracking via Praat's SPINET method
#'
#' Tracks fundamental frequency (F0) using Praat's SPINET (SPectral INtegration
#' and Evaluation of Temporal patterns) algorithm via pladdrr. SPINET uses a
#' bank of gammatone filters followed by temporal integration, mimicking auditory
#' processing; it can detect F0 in conditions where time-domain methods fail.
#'
#' @param listOfFiles Character vector of audio file paths. Any format supported by
#'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
#' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
#' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
#' @param time_step Numeric. Frame shift in seconds; sets output frame rate
#'   (1 / time_step Hz). Default 0.005 s (200 Hz).
#' @param window_length Numeric. Duration of the integration window in seconds.
#'   Default 0.04 s.
#' @param minimum_filter_frequency Numeric. Center frequency of the lowest gammatone
#'   filter in Hz. Default 70 Hz.
#' @param maximum_filter_frequency Numeric. Center frequency of the highest gammatone
#'   filter in Hz. Default 5000 Hz.
#' @param number_of_filters Integer. Number of gammatone filters in the bank.
#'   Default 250.
#' @param maximum_f0 Numeric. Upper F0 bound (ceiling) in Hz. Default 500 Hz.
#' @param number_of_candidates Integer. Maximum pitch candidates per frame.
#'   Default 15.
#' @param windowShape Character. Window shape applied to audio before loading.
#'   Default \code{"Gaussian1"}.
#' @param relativeWidth Numeric. Relative width of the window. Default 1.0.
#' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
#'   count written (invisibly). If \code{FALSE}, return an \code{AsspDataObj}.
#'   Default \code{TRUE}.
#' @param explicitExt Character. Output file extension. Default \code{"psp"}.
#' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
#'   writes alongside the input file.
#' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
#'
#' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
#'   \describe{
#'     \item{\code{F0}}{REAL32, Hz, n_frames x 1. Fundamental frequency.
#'       0 encodes unvoiced frames.}
#'   }
#'   Frame rate: \code{1 / time_step} Hz (default 200 Hz).
#'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
#'
#' @details
#' SPINET is slower than AC/CC/SHS due to the gammatone filterbank but can
#' outperform them on atypical voices (very low F0, strong noise).
#' No lower F0 bound is exposed — the minimum detectable F0 is determined
#' by \code{minimum_filter_frequency}.
#'
#' @seealso \code{\link{trk_pitch_cc}}, \code{\link{trk_pitch_ac}},
#'   \code{\link{trk_pitch_shs}}
#'
#' @examples
#' \dontrun{
#' trk_pitch_spinet(
#'   system.file("samples", "sustained", "a1.wav", package = "superassp"),
#'   toFile = FALSE
#' )
#' }
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

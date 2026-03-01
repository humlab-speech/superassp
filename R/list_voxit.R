#' Extract Voxit prosodic complexity features from audio files
#'
#' @description
#' Computes voice and articulation complexity measures from audio files using
#' word-level alignments and pitch contours. Features include speaking rate,
#' pause statistics, rhythmic complexity, and pitch dynamics.
#'
#' This is a pure R/C++ implementation requiring no Python dependencies.
#' Pitch tracking uses \code{\link{trk_rapt}} (SPTK C++).
#'
#' @details
#' \strong{Audio Loading:}
#'
#' Audio files are loaded via the \code{\link[av]{read_audio_bin}} function,
#' supporting all media formats (WAV, MP3, MP4, video files, etc.).
#'
#' \strong{Features Computed (11 total):}
#'
#' \enumerate{
#'   \item \strong{WPM}: Words per minute (speaking rate)
#'   \item \strong{pause_count}: Number of pauses (100-3000ms)
#'   \item \strong{long_pause_count}: Number of pauses > 3s
#'   \item \strong{average_pause_length}: Mean pause duration (seconds)
#'   \item \strong{average_pause_rate}: Pauses per second
#'   \item \strong{rhythmic_complexity_of_pauses}: Normalized Lempel-Ziv complexity
#'   \item \strong{average_pitch}: Mean F0 (Hz)
#'   \item \strong{pitch_range}: F0 range (octaves)
#'   \item \strong{pitch_speed}: F0 velocity (octaves/second)
#'   \item \strong{pitch_acceleration}: F0 acceleration (octaves/second^2)
#'   \item \strong{pitch_entropy}: F0 distribution entropy (bits)
#' }
#'
#' @param listOfFiles Character vector of file paths to audio files
#' @param alignmentFiles Character vector of alignment CSV file paths.
#'   CSV must have columns: word, start, end (and optionally case).
#' @param beginTime Numeric. Start time in seconds (default: 0 = file start)
#' @param endTime Numeric. End time in seconds (default: 0 = file end)
#' @param minF Numeric. Minimum F0 for pitch tracking in Hz (default: 60)
#' @param maxF Numeric. Maximum F0 for pitch tracking in Hz (default: 600)
#' @param verbose Logical. Show progress messages (default: TRUE)
#' @param parallel Logical. Use parallel processing for multiple files (default: TRUE)
#' @param n_cores Integer. Number of cores for parallel processing (default: NULL = auto-detect)
#' @param toFile Logical. If TRUE, write results to JSTF file. Default FALSE.
#' @param explicitExt Character. File extension for output. Default "vxt".
#' @param outputDirectory Character. Output directory path. Default NULL (use input directory).
#'
#' @return
#' If \code{toFile=FALSE} (default): For single file, a named list with 11 prosodic features.
#'   For multiple files, named list where each element is a file's feature list.
#'   If \code{toFile=TRUE}, invisibly returns the path(s) to the written JSTF file(s).
#'
#' @section Word Alignments:
#' Word alignments must be provided as CSV files with columns:
#' \describe{
#'   \item{word}{Word text}
#'   \item{start}{Start time in seconds}
#'   \item{end}{End time in seconds}
#' }
#'
#' @examples
#' \dontrun{
#' # Single file analysis with alignments
#' features <- lst_voxit("speech.wav",
#'                       alignmentFiles = "speech_align.csv")
#' print(features$WPM)
#' print(features$average_pitch)
#'
#' # Write results to JSTF file
#' lst_voxit("speech.wav", alignmentFiles = "speech_align.csv", toFile = TRUE)
#' }
#'
#' @references
#' Voxit toolbox: Voice and articulation complexity measures
#'
#' @seealso \code{\link{trk_rapt}}
#'
#' @export
lst_voxit <- function(listOfFiles,
                      alignmentFiles = NULL,
                      beginTime = 0.0,
                      endTime = 0.0,
                      minF = 60,
                      maxF = 600,
                      verbose = TRUE,
                      parallel = TRUE,
                      n_cores = NULL,
                      toFile = FALSE,
                      explicitExt = "vxt",
                      outputDirectory = NULL) {

  # Input validation
  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No files provided")
  }

  listOfFiles <- as.vector(listOfFiles)
  n_files <- length(listOfFiles)

  missing <- !file.exists(listOfFiles)
  if (any(missing)) {
    cli::cli_abort(c(
      "x" = "{sum(missing)} file{?s} not found",
      "i" = "First missing: {.file {listOfFiles[which(missing)[1]]}}"
    ))
  }

  if (is.null(alignmentFiles)) {
    cli::cli_abort(c(
      "x" = "Alignment files required",
      "i" = "Provide alignmentFiles parameter (CSV with word, start, end columns)"
    ))
  }

  if (length(alignmentFiles) != n_files) {
    cli::cli_abort(c(
      "x" = "Number of alignment files ({length(alignmentFiles)}) must match audio files ({n_files})"
    ))
  }

  missing_align <- !file.exists(alignmentFiles)
  if (any(missing_align)) {
    cli::cli_abort(c(
      "x" = "{sum(missing_align)} alignment file{?s} not found",
      "i" = "First missing: {.file {alignmentFiles[which(missing_align)[1]]}}"
    ))
  }

  # Single-file processor
  process_single_file <- function(i) {
    file_path <- listOfFiles[i]
    align_path <- alignmentFiles[i]

    if (verbose && !parallel) {
      cli::cli_progress_step("Processing {.file {basename(file_path)}}")
    }

    tryCatch({
      # Read alignments
      align_df <- readr::read_csv(align_path, show_col_types = FALSE)

      # Filter valid words
      valid <- !is.na(align_df$start) & !is.na(align_df$end)
      if ("word" %in% names(align_df)) {
        valid <- valid & (align_df$word != "[noise]")
      }

      # Apply time windowing
      if (beginTime > 0) valid <- valid & (align_df$start >= beginTime)
      if (endTime > 0) valid <- valid & (align_df$end <= endTime)

      align_df <- align_df[valid, , drop = FALSE]

      if (nrow(align_df) == 0) {
        return(.voxit_nan_result())
      }

      starts <- align_df$start
      ends <- align_df$end
      word_count <- nrow(align_df)
      gentle_length <- ends[word_count]  # end time of last word

      if (gentle_length <= 0) return(.voxit_nan_result())

      # --- WPM ---
      WPM <- as.integer(word_count / (gentle_length / 60))

      # --- Pause analysis ---
      pauses <- .voxit_pause_analysis(starts, ends, gentle_length)

      # --- Rhythmic complexity (Lempel-Ziv) ---
      rhythmic_complexity <- .voxit_rhythmic_complexity(starts, ends, gentle_length)

      # --- Pitch features ---
      pitch_features <- .voxit_pitch_features(
        file_path, beginTime, endTime, minF, maxF, starts, ends
      )

      result <- c(
        list(WPM = WPM),
        pauses,
        list(rhythmic_complexity_of_pauses = rhythmic_complexity),
        pitch_features
      )

      return(result)

    }, error = function(e) {
      if (verbose) {
        cli::cli_alert_warning("Failed {.file {basename(file_path)}}: {e$message}")
      }
      return(NULL)
    })
  }

  # Process files
  if (n_files == 1) {
    result <- process_single_file(1)

    if (toFile && !is.null(result)) {
      return(invisible(.voxit_write_jstf(
        result, listOfFiles[1], beginTime, endTime, minF, maxF,
        explicitExt, outputDirectory
      )))
    }

  } else if (parallel && n_files > 1) {
    if (is.null(n_cores)) n_cores <- parallel::detectCores() - 1
    n_cores <- min(n_cores, n_files)

    if (verbose) {
      cli::cli_alert_info("Processing {n_files} files using {n_cores} cores")
    }

    if (.Platform$OS.type == "windows") {
      cl <- parallel::makeCluster(n_cores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      result <- parallel::parLapply(cl, 1:n_files, process_single_file)
    } else {
      result <- parallel::mclapply(1:n_files, process_single_file,
                                   mc.cores = n_cores)
    }
    names(result) <- basename(listOfFiles)
  } else {
    result <- lapply(1:n_files, process_single_file)
    names(result) <- basename(listOfFiles)
  }

  # JSTF writing for multi-file
  if (toFile && n_files > 1) {
    output_paths <- vapply(seq_len(n_files), function(i) {
      r <- result[[basename(listOfFiles[i])]]
      if (!is.null(r)) {
        .voxit_write_jstf(r, listOfFiles[i], beginTime, endTime, minF, maxF,
                          explicitExt, outputDirectory)
      } else NA_character_
    }, character(1))
    return(invisible(output_paths))
  }

  if (n_files == 1) return(result)
  return(result)
}


# --- Internal helpers ---

#' @noRd
.voxit_nan_result <- function() {
  list(
    WPM = NaN, pause_count = NaN, long_pause_count = NaN,
    average_pause_length = NaN, average_pause_rate = NaN,
    rhythmic_complexity_of_pauses = NaN, average_pitch = NaN,
    pitch_range = NaN, pitch_speed = NaN, pitch_acceleration = NaN,
    pitch_entropy = NaN
  )
}

#' @noRd
.voxit_pause_analysis <- function(starts, ends, gentle_length) {
  min_pause <- 0.1
  max_pause <- 3.0

  if (length(starts) < 2) {
    return(list(
      pause_count = 0L, long_pause_count = 0L,
      average_pause_length = 0.0, average_pause_rate = 0.0
    ))
  }

  gaps <- starts[-1] - ends[-length(ends)]
  valid_pauses <- gaps[gaps >= min_pause & gaps <= max_pause]
  long_pauses <- gaps[gaps > max_pause]

  pause_count <- length(valid_pauses)
  long_pause_count <- length(long_pauses)

  list(
    pause_count = pause_count,
    long_pause_count = long_pause_count,
    average_pause_length = if (pause_count > 0) mean(valid_pauses) else 0.0,
    average_pause_rate = if (pause_count > 0) pause_count / gentle_length else 0.0
  )
}

#' @noRd
.voxit_rhythmic_complexity <- function(starts, ends, gentle_length) {
  min_pause <- 0.1
  max_pause <- 3.0
  sampling_interval <- 0.01

  n_samples <- as.integer(gentle_length / sampling_interval) + 1L
  if (n_samples < 2) return(0.0)

  s <- rep(1L, n_samples)  # 1 = voiced/non-pause

  if (length(starts) >= 2) {
    gaps <- starts[-1] - ends[-length(ends)]
    for (k in seq_along(gaps)) {
      if (gaps[k] >= min_pause && gaps[k] <= max_pause) {
        start_idx <- as.integer(ends[k] / sampling_interval) + 1L
        end_idx <- as.integer(starts[k + 1] / sampling_interval)
        if (end_idx >= start_idx && start_idx >= 1 && end_idx <= n_samples) {
          s[start_idx:end_idx] <- 0L
        }
      }
    }
  }

  lz <- .lempel_ziv_complexity(s)
  normalized <- lz / (n_samples / log2(n_samples))
  return(normalized * 100)
}

#' Lempel-Ziv complexity of a binary sequence
#' Matches Python lempel_ziv_complexity fallback implementation
#' @noRd
.lempel_ziv_complexity <- function(s) {
  # Works on character string (matching Python API)
  s_str <- paste0(s, collapse = "")
  n <- nchar(s_str)
  if (n == 0) return(0L)

  i <- 1L  # R 1-indexed (Python i=0)
  k <- 1L
  l <- 2L  # R 1-indexed (Python l=1)
  c_val <- 1L
  k_max <- 1L

  while (TRUE) {
    if (i + k - 1L > n) return(c_val)
    sub_ik <- substr(s_str, i, i + k - 1L)
    sub_lk <- substr(s_str, l, l + k - 1L)
    if (sub_ik != sub_lk) {
      c_val <- c_val + 1L
      i <- l + k
      l <- l + k
      k <- 1L
      k_max <- 1L
    } else {
      k <- k + 1L
      if (k > k_max) k_max <- k
      if (i + k - 1L > n) return(c_val + 1L)
    }
  }
}

#' @noRd
.voxit_pitch_features <- function(file_path, beginTime, endTime, minF, maxF,
                                   word_starts, word_ends) {
  nan_result <- list(
    average_pitch = NaN, pitch_range = NaN,
    pitch_speed = NaN, pitch_acceleration = NaN,
    pitch_entropy = NaN
  )

  # Get F0 via trk_rapt (C++)
  f0_obj <- tryCatch(
    trk_rapt(file_path,
             beginTime = beginTime, endTime = endTime,
             minF = minF, maxF = maxF,
             toFile = FALSE, verbose = FALSE),
    error = function(e) NULL
  )

  if (is.null(f0_obj)) return(nan_result)

  f0_vals <- as.numeric(f0_obj[["pitch[Hz]"]])
  sr <- attr(f0_obj, "sampleRate")
  start_time_attr <- attr(f0_obj, "startTime")
  if (is.null(start_time_attr)) start_time_attr <- 0

  # Compute times for each frame
  n_frames <- length(f0_vals)
  frame_shift <- 1.0 / sr  # sr here is the track sample rate (frames/sec)
  times <- start_time_attr + (seq_len(n_frames) - 1) * frame_shift

  # Filter voiced frames (f0 > 0)
  voiced <- f0_vals > 0 & !is.na(f0_vals)
  if (sum(voiced) < 3) return(nan_result)

  drift_time <- times[voiced]
  drift_pitch <- f0_vals[voiced]

  # average_pitch
  average_pitch <- mean(drift_pitch)

  # pitch_range (in octaves)
  f0log <- log2(drift_pitch)
  f0mean_geom <- 2^mean(f0log)  # geometric mean
  diffoctf0 <- log2(drift_pitch) - log2(f0mean_geom)
  pitch_range <- max(diffoctf0) - min(diffoctf0)

  # pitch_entropy (Shannon, 25 bins over +-1 octave)
  h <- hist(diffoctf0, breaks = seq(-1, 1, length.out = 26), plot = FALSE)
  f0prob <- h$counts / sum(h$counts)
  f0prob_nz <- f0prob[f0prob > 0]
  pitch_entropy <- -sum(f0prob_nz * log2(f0prob_nz))

  # pitch_speed and pitch_acceleration via Savitzky-Golay
  pitch_speed <- 0.0
  pitch_acceleration <- 0.0

  ts <- if (length(drift_time) > 1) drift_time[2] - drift_time[1] else 0.01

  # Segment voiced regions (break at gaps > 3*ts)
  dminvoice <- 0.100  # 100 ms minimum voiced segment
  vdurthresh <- as.integer(dminvoice / ts)
  gap_threshold <- 3 * ts

  if (length(drift_time) > 1) {
    time_diff <- diff(drift_time)
    gap_indices <- which(time_diff > gap_threshold)

    seg_starts <- c(1L, gap_indices + 1L)
    seg_ends <- c(gap_indices, length(drift_pitch))

    all_velocity <- numeric(0)
    all_accel <- numeric(0)

    for (seg_i in seq_along(seg_starts)) {
      si <- seg_starts[seg_i]
      ei <- seg_ends[seg_i]
      seg_len <- ei - si + 1

      if (seg_len <= vdurthresh) next

      seg_oct <- diffoctf0[si:ei]

      # Savitzky-Golay smoothing (window=7, order=2)
      if (seg_len >= 7 && requireNamespace("gsignal", quietly = TRUE)) {
        seg_smooth <- gsignal::sgolayfilt(seg_oct, p = 2, n = 7)
      } else {
        seg_smooth <- seg_oct
      }

      vel <- diff(seg_smooth) / ts
      acc <- diff(diff(seg_smooth)) / ts  # matches Python: single division

      all_velocity <- c(all_velocity, vel)
      all_accel <- c(all_accel, acc)
    }

    if (length(all_velocity) > 0) {
      pitch_speed <- mean(abs(all_velocity)) * sign(mean(all_velocity))
    }
    if (length(all_accel) > 0) {
      pitch_acceleration <- mean(abs(all_accel)) * sign(mean(all_accel))
    }
  }

  list(
    average_pitch = average_pitch,
    pitch_range = pitch_range,
    pitch_speed = pitch_speed,
    pitch_acceleration = pitch_acceleration,
    pitch_entropy = pitch_entropy
  )
}

#' @noRd
.voxit_write_jstf <- function(result, file_path, beginTime, endTime, minF, maxF,
                               explicitExt, outputDirectory) {
  audio_info <- av::av_media_info(file_path)
  sample_rate <- audio_info$audio$sample_rate
  audio_duration <- audio_info$duration

  json_obj <- create_json_track_obj(
    results = result,
    function_name = "lst_voxit",
    file_path = file_path,
    sample_rate = sample_rate,
    audio_duration = audio_duration,
    beginTime = beginTime,
    endTime = if (endTime > 0) endTime else audio_duration,
    parameters = list(minF = minF, maxF = maxF)
  )

  base_name <- tools::file_path_sans_ext(basename(file_path))
  out_dir <- if (is.null(outputDirectory)) dirname(file_path) else outputDirectory
  output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))

  write_json_track(json_obj, output_path)
  return(output_path)
}

# Set function attributes
attr(lst_voxit, "module") <- "voxit"
attr(lst_voxit, "type") <- "summary"
attr(lst_voxit, "features") <- 11
attr(lst_voxit, "ext") <- "vxt"
attr(lst_voxit, "outputType") <- "JSTF"
attr(lst_voxit, "format") <- "JSON"

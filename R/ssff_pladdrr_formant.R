#' Formant frequencies and bandwidths via Praat's Burg method
#'
#' Tracks formant frequencies (F1–F5), bandwidths (B1–B5), and optionally
#' spectral intensities (L1–L5) using Praat's Burg LPC algorithm via pladdrr.
#' Optional HMM-based formant tracking smooths trajectories across time. Prefer
#' this over \code{trk_formant_forest} when Praat-compatible Burg estimates are needed.
#'
#' @param listOfFiles Character vector of audio file paths. Any format supported by
#'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
#' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
#' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
#' @param timeStep Numeric. Frame shift in seconds; sets output frame rate
#'   (1 / timeStep Hz). Default 0.005 s (200 Hz).
#' @param number_of_formants Integer. Number of formant candidates to find (up to 5).
#'   Default 5.
#' @param maxHzFormant Numeric. Formant ceiling in Hz; typically 5500 for female,
#'   5000 for male speakers. Default 5500 Hz.
#' @param windowLength Numeric. LPC analysis window length in seconds. Default 0.025 s.
#' @param pre_emphasis Numeric. Pre-emphasis onset frequency in Hz. Default 50 Hz.
#' @param track_formants Logical. Apply HMM-based formant tracking to smooth
#'   trajectories. Default \code{FALSE}.
#' @param number_of_tracks Integer. Number of tracks to retain when
#'   \code{track_formants = TRUE}. Default 3.
#' @param reference_F1 Numeric. Reference F1 frequency for HMM tracking in Hz.
#'   Default 550 Hz.
#' @param reference_F2 Numeric. Reference F2 frequency for HMM tracking in Hz.
#'   Default 1650 Hz.
#' @param reference_F3 Numeric. Reference F3 frequency for HMM tracking in Hz.
#'   Default 2750 Hz.
#' @param reference_F4 Numeric. Reference F4 frequency for HMM tracking in Hz.
#'   Default 3850 Hz.
#' @param reference_F5 Numeric. Reference F5 frequency for HMM tracking in Hz.
#'   Default 4950 Hz.
#' @param frequency_cost Numeric. HMM cost weight for deviation from reference
#'   frequencies. Default 1.0.
#' @param bandwidth_cost Numeric. HMM cost weight for formant bandwidth. Default 1.0.
#' @param transition_cost Numeric. HMM cost weight for frame-to-frame transitions.
#'   Default 1.0.
#' @param windowShape Character. Window shape applied before analysis. Default
#'   \code{"Gaussian1"}.
#' @param relativeWidth Numeric. Relative width of the analysis window. Default 1.0.
#' @param include_intensity Logical. Extract spectral intensity at each formant
#'   frequency (L1–L5 tracks via spectrogram lookup). Default \code{TRUE}.
#' @param spectrogram_resolution Numeric. Frequency resolution of the spectrogram
#'   used for intensity extraction in Hz. Default 40 Hz.
#' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
#'   count written (invisibly). If \code{FALSE}, return an \code{AsspDataObj}.
#'   Default \code{TRUE}.
#' @param explicitExt Character. Output file extension. Default \code{"pfm"}.
#' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
#'   writes alongside the input file.
#' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
#'
#' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
#'   \describe{
#'     \item{\code{fm1}–\code{fm5}}{REAL32, Hz, n_frames x 1. Formant frequencies F1–F5.
#'       0 encodes an undefined (unvoiced) frame.}
#'     \item{\code{bw1}–\code{bw5}}{REAL32, Hz, n_frames x 1. Formant bandwidths B1–B5.
#'       0 encodes an undefined frame.}
#'     \item{\code{L1}–\code{L5}}{REAL32, Pa^2/Hz, n_frames x 1. Spectral power at each
#'       formant frequency (only when \code{include_intensity = TRUE}).}
#'   }
#'   Frame rate: \code{1 / timeStep} Hz (default 200 Hz).
#'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
#'
#' @details
#' Uses Praat's Burg all-pole LPC estimator. Spectral intensity tracks (L1–L5)
#' are extracted from a separate spectrogram and require pladdrr >= 4.8.20;
#' earlier versions may crash with \code{include_intensity = TRUE}.
#' HMM tracking (\code{track_formants = TRUE}) may not be available in all
#' pladdrr builds; a warning is issued and untracked formants are returned.
#'
#' @examples
#' \dontrun{
#' # Analyze formants for a single file
#' formants <- trk_formant_burg("speech.wav", toFile = FALSE)
#'
#' # Batch process multiple files with tracking
#' files <- c("speech1.wav", "speech2.wav")
#' trk_formant_burg(files, track_formants = TRUE, toFile = TRUE)
#'
#' # With custom formant ceiling (e.g., male speaker)
#' formants <- trk_formant_burg("speech.wav",
#'                          maxHzFormant = 5000,
#'                          number_of_formants = 4,
#'                          toFile = FALSE)
#' }
#'
#' @export
trk_formant_burg <- function(listOfFiles,
                         beginTime = 0.0,
                         endTime = 0.0,
                         timeStep = 0.005,
                         number_of_formants = 5,
                         maxHzFormant = 5500.0,
                         windowLength = 0.025,
                         pre_emphasis = 50.0,
                         track_formants = FALSE,
                         number_of_tracks = 3,
                         reference_F1 = 550,
                         reference_F2 = 1650,
                         reference_F3 = 2750,
                         reference_F4 = 3850,
                         reference_F5 = 4950,
                         frequency_cost = 1.0,
                         bandwidth_cost = 1.0,
                         transition_cost = 1.0,
                         windowShape = "Gaussian1",
                         relativeWidth = 1.0,
                         include_intensity = TRUE,   # Now enabled by default (fixed in pladdrr 4.8.20+)
                         spectrogram_resolution = 40.0,
                         toFile = TRUE,
                         explicitExt = "pfm",
                         outputDirectory = NULL,
                         verbose = TRUE) {

  # Check pladdrr availability
  if (!pladdrr_available()) {
    cli::cli_abort(c(
      "x" = "pladdrr package not available",
      "i" = "Install with: install_pladdrr()"
    ))
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

    # Extract formants using Burg's method - Direct API
    formant_ptr <- pladdrr::to_formant_direct(
      analysis_sound,
      time_step = timeStep,
      max_formants = as.numeric(number_of_formants),
      max_formant = maxHzFormant,
      window_length = windowLength,
      pre_emphasis = pre_emphasis
    )
    formant <- pladdrr::Formant(.xptr = formant_ptr)

    actual_num_formants <- as.integer(number_of_formants)

    # Apply tracking if requested
    if (track_formants) {
      track_outcome <- tryCatch({
        formant <- formant$track(
          as.integer(number_of_tracks),
          reference_F1, reference_F2, reference_F3, reference_F4, reference_F5,
          frequency_cost, bandwidth_cost, transition_cost
        )
        actual_num_formants <- min(actual_num_formants, as.integer(number_of_tracks))
        list(ok = TRUE, formant = formant, actual = actual_num_formants)
      }, error = function(e) {
        if (verbose) {
          warning("Formant tracking not available or failed, using untracked formants")
        }
        list(ok = FALSE)
      })
      if (track_outcome$ok) {
        formant             <- track_outcome$formant
        actual_num_formants <- track_outcome$actual
      } else {
        track_formants <- FALSE
      }
    }

    # Create spectrogram for intensity extraction if needed
    spectrogram <- NULL
    if (include_intensity) {
      spec_outcome <- tryCatch({
        max_hz <- maxHzFormant + 2000.0
        spec <- analysis_sound$to_spectrogram(
          window_length  = windowLength,
          max_frequency  = max_hz,
          time_step      = timeStep,
          frequency_step = spectrogram_resolution,
          window_shape   = "Gaussian"
        )
        list(ok = TRUE, spec = spec)
      }, error = function(e) {
        if (verbose) {
          warning("Spectrogram creation failed, skipping intensity extraction")
        }
        list(ok = FALSE)
      })
      if (spec_outcome$ok) {
        spectrogram <- spec_outcome$spec
      } else {
        include_intensity <- FALSE
      }
    }

    # Get frame data using as_data_frame
    # pladdrr returns long format: time, formant, frequency, bandwidth
    formant_df_long <- formant$as_data_frame()

    # Get unique times
    times <- sort(unique(formant_df_long$time))
    num_frames <- length(times)

    # Initialize result data frame
    result_list <- list(time = times)

    # Extract formant frequencies and bandwidths for each formant
    for (f_num in seq_len(actual_num_formants)) {
      f_col <- paste0("fm", f_num)
      b_col <- paste0("bw", f_num)

      # Extract values for this formant number from long format
      f_data <- formant_df_long[formant_df_long$formant == f_num, ]

      # Initialize vectors
      f_vals <- rep(NA_real_, num_frames)
      b_vals <- rep(NA_real_, num_frames)

      # Match times
      time_idx <- match(times, f_data$time)
      valid_idx <- !is.na(time_idx)

      if (any(valid_idx)) {
        f_vals[valid_idx] <- f_data$frequency[time_idx[valid_idx]]
        b_vals[valid_idx] <- f_data$bandwidth[time_idx[valid_idx]]
      }

      # Replace 0 with NA (Praat convention for undefined formants)
      f_vals[f_vals == 0 | !is.finite(f_vals)] <- NA
      b_vals[b_vals == 0 | !is.finite(b_vals)] <- NA

      result_list[[f_col]] <- f_vals
      result_list[[b_col]] <- b_vals

      # Get intensity at formant frequency
      if (include_intensity && !is.null(spectrogram)) {
        l_col <- paste0("L", f_num)
        l_vals <- rep(NA_real_, num_frames)

        for (j in seq_len(num_frames)) {
          if (!is.na(f_vals[j])) {
            val <- tryCatch(
              spectrogram$get_power_at(times[j], f_vals[j]),
              error = function(e) NA
            )
            l_vals[j] <- if (is.null(val) || is.na(val)) NA else val
          }
        }
        result_list[[l_col]] <- l_vals
      }
    }

    # Convert to data frame
    result_df <- as.data.frame(result_list)

    # Create AsspDataObj
    outDataObj <- list()

    # Calculate sample rate from formant object
    actual_time_step <- formant$get_time_step()
    sampleRate <- 1 / actual_time_step

    attr(outDataObj, "sampleRate") <- sampleRate
    attr(outDataObj, "origFreq") <- as.numeric(sound$.cpp$sampling_frequency)

    startTime <- times[1]
    attr(outDataObj, "startTime") <- as.numeric(startTime)
    attr(outDataObj, "startRecord") <- as.integer(1)
    attr(outDataObj, "endRecord") <- as.integer(num_frames)

    class(outDataObj) <- "AsspDataObj"
    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2)
    attr(outDataObj, "trackFormats") <- character(0)

    # Add formant frequency tracks (fm1-fm5)
    for (f_num in seq_len(actual_num_formants)) {
      f_col <- paste0("fm", f_num)
      f_data <- result_df[[f_col]]
      f_data[is.na(f_data)] <- 0  # SSFF uses 0 for missing values
      outDataObj <- wrassp::addTrack(outDataObj, f_col,
                                     as.matrix(f_data), "REAL32")
      attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
    }

    # Add bandwidth tracks (bw1-bw5)
    for (f_num in seq_len(actual_num_formants)) {
      b_col <- paste0("bw", f_num)
      b_data <- result_df[[b_col]]
      b_data[is.na(b_data)] <- 0
      outDataObj <- wrassp::addTrack(outDataObj, b_col,
                                     as.matrix(b_data), "REAL32")
      attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
    }

    # Add intensity tracks (L1-L5) if available
    if (include_intensity) {
      for (f_num in seq_len(actual_num_formants)) {
        l_col <- paste0("L", f_num)
        if (l_col %in% names(result_df)) {
          l_data <- result_df[[l_col]]
          l_data[is.na(l_data)] <- 0
          outDataObj <- wrassp::addTrack(outDataObj, l_col,
                                         as.matrix(l_data), "REAL32")
          attr(outDataObj, "trackFormats") <- c(attr(outDataObj, "trackFormats"), "REAL32")
        }
      }
    }

    # Write file or store result
    if(toFile) {
      base_name <- tools::file_path_sans_ext(basename(origSoundFile))
      out_dir <- if (is.null(outputDirectory)) dirname(origSoundFile) else outputDirectory

      # Create output directory if needed
      if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
      }

      output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))
      wrassp::write.AsspDataObj(outDataObj, output_path)
      outListOfFiles <- c(outListOfFiles, output_path)
    } else {
      outListOfFiles <- outDataObj
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

  # Return results
  if(toFile) {
    return(invisible(n_files))
  } else {
    return(outListOfFiles)
  }
}

# Set function attributes
attr(trk_formant_burg, "ext") <- "pfm"
attr(trk_formant_burg, "tracks") <- c("fm1", "fm2", "fm3", "fm4", "fm5",
                                   "bw1", "bw2", "bw3", "bw4", "bw5",
                                   "L1", "L2", "L3", "L4", "L5")
attr(trk_formant_burg, "outputType") <- "SSFF"
attr(trk_formant_burg, "nativeFiletypes") <- c("wav")

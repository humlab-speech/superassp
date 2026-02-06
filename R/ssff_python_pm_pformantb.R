#' Optimized formant analysis using pladdrr (Burg method)
#'
#' Computes formant frequencies, bandwidths, and spectral intensities using pladdrr's
#' native R bindings to Praat's C library. Uses Burg's algorithm for formant extraction
#' with optional HMM tracking. Audio is loaded directly by pladdrr, eliminating temporary
#' file creation. Supports all media formats via pladdrr's native readers.
#'
#' @param listOfFiles Vector of file paths to audio files
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param timeStep Time step between frames in seconds
#' @param number_of_formants Maximum number of formants to extract
#' @param maxHzFormant Maximum formant frequency in Hz
#' @param windowLength Analysis window length in seconds
#' @param pre_emphasis Pre-emphasis frequency in Hz
#' @param track_formants Whether to apply HMM formant tracking
#' @param number_of_tracks Number of formant tracks (used if track_formants=TRUE)
#' @param reference_F1 Reference frequency for F1 in Hz (for tracking)
#' @param reference_F2 Reference frequency for F2 in Hz (for tracking)
#' @param reference_F3 Reference frequency for F3 in Hz (for tracking)
#' @param reference_F4 Reference frequency for F4 in Hz (for tracking)
#' @param reference_F5 Reference frequency for F5 in Hz (for tracking)
#' @param frequency_cost Weight for frequency deviation in tracking
#' @param bandwidth_cost Weight for bandwidth in tracking
#' @param transition_cost Weight for formant transitions in tracking
#' @param windowShape Window shape for time windowing
#' @param relativeWidth Relative width for windowing
#' @param include_intensity Whether to extract spectral intensity at formant frequencies
#' @param spectrogram_resolution Frequency resolution for spectrogram in Hz
#' @param toFile Write output to file (TRUE) or return object (FALSE)
#' @param explicitExt File extension for output files
#' @param outputDirectory Output directory path (NULL for same as input)
#' @param verbose Show progress messages (default: TRUE)
#'
#' @return If toFile=TRUE, returns number of files successfully processed.
#'   If toFile=FALSE, returns an AsspDataObj with formant tracks.
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
#' The function extracts:
#' \itemize{
#'   \item **Formant frequencies**: F1-F5 (fm1-fm5 tracks)
#'   \item **Bandwidths**: B1-B5 (bw1-bw5 tracks)
#'   \item **Spectral intensities**: L1-L5 (optional, if include_intensity=TRUE)
#' }
#'
#' **Formant tracking**: If track_formants=TRUE, applies HMM-based tracking
#' to smooth formant trajectories across time. May not be fully available in
#' all pladdrr versions.
#'
#' **Note on formant bug**: pladdrr v4.6.4 had a known bug where F1/F2/F3 values
#' were systematically 35-55% too low. This appears to be FIXED in v4.8.16+ based
#' on testing with sustained vowels. Formant values now match expected ranges.
#'
#' **Known limitation**: Spectral intensity extraction (include_intensity=TRUE) may
#' cause crashes in some pladdrr versions due to spectrogram issues. Set
#' include_intensity=FALSE if you encounter segfaults.
#'
#' @examples
#' \dontrun{
#' # Analyze formants for a single file
#' formants <- trk_formantp("speech.wav", toFile = FALSE)
#'
#' # Batch process multiple files with tracking
#' files <- c("speech1.wav", "speech2.wav")
#' trk_formantp(files, track_formants = TRUE, toFile = TRUE)
#'
#' # With custom formant ceiling (e.g., male speaker)
#' formants <- trk_formantp("speech.wav",
#'                          maxHzFormant = 5000,
#'                          number_of_formants = 4,
#'                          toFile = FALSE)
#' }
#'
#' @export
trk_formantp <- function(listOfFiles,
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
                         include_intensity = FALSE,  # Disabled by default due to pladdrr spectrogram issues
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
      tryCatch({
        formant <- formant$track(
          as.integer(number_of_tracks),
          reference_F1, reference_F2, reference_F3, reference_F4, reference_F5,
          frequency_cost, bandwidth_cost, transition_cost
        )
        actual_num_formants <- min(actual_num_formants, as.integer(number_of_tracks))
      }, error = function(e) {
        if (verbose) {
          warning("Formant tracking not available or failed, using untracked formants")
        }
        track_formants <<- FALSE
      })
    }

    # Create spectrogram for intensity extraction if needed
    spectrogram <- NULL
    if (include_intensity) {
      tryCatch({
        max_hz <- maxHzFormant + 2000.0  # Avoid edge effects
        spectrogram <- analysis_sound$to_spectrogram(
          window_length = windowLength,
          max_frequency = max_hz,
          time_step = timeStep,
          frequency_step = spectrogram_resolution,
          window_shape = "Gaussian"
        )
      }, error = function(e) {
        if (verbose) {
          warning("Spectrogram creation failed, skipping intensity extraction")
        }
        include_intensity <<- FALSE
      })
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
attr(trk_formantp, "ext") <- "pfm"
attr(trk_formantp, "tracks") <- c("fm1", "fm2", "fm3", "fm4", "fm5",
                                   "bw1", "bw2", "bw3", "bw4", "bw5",
                                   "L1", "L2", "L3", "L4", "L5")
attr(trk_formantp, "outputType") <- "SSFF"
attr(trk_formantp, "nativeFiletypes") <- c("wav")

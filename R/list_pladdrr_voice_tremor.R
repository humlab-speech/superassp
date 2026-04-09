#' Vocal Tremor Analysis Using pladdrr
#'
#' Analyzes vocal tremor from sustained vowel recordings using pladdrr's Praat
#' bindings. Extracts 18 measures of frequency and amplitude tremor based on
#' Brückl (2012) autocorrelation algorithm.
#'
#' This function processes sustained phonations to detect tremor characteristics
#' in both pitch (frequency) and intensity (amplitude) contours. It applies
#' Gaussian1 windowing and uses autocorrelation-based analysis to identify
#' tremor frequency, intensity, and cyclicality.
#'
#' @param listOfFiles Character vector with path(s) to audio file(s)
#' @param beginTime Numeric. Start time in seconds (default 0)
#' @param endTime Numeric. End time in seconds (0 = end of file)
#' @param analysisTimeStep Numeric. Time step for analysis in seconds (default 0.015)
#' @param minPitch Numeric. Minimum pitch for extraction in Hz (default 60)
#' @param maxPitch Numeric. Maximum pitch for extraction in Hz (default 350)
#' @param silenceThreshold Numeric. Threshold for silence detection (default 0.03)
#' @param voicingThreshold Numeric. Threshold for voicing detection (default 0.3)
#' @param octaveCost Numeric. Cost for octave jumps in pitch tracking (default 0.01)
#' @param octaveJumpCost Numeric. Cost for large octave jumps (default 0.35)
#' @param voicedUnvoicedCost Numeric. Cost for voiced/unvoiced transitions (default 0.14)
#' @param minTremorFreq Numeric. Minimum tremor frequency in Hz (default 1.5)
#' @param maxTremorFreq Numeric. Maximum tremor frequency in Hz (default 15)
#' @param tremorMagThresh Numeric. Threshold for contour magnitude (default 0.01)
#' @param tremorCyclicalThresh Numeric. Threshold for cyclicality (default 0.15)
#' @param freqTremorOctaveCost Numeric. Octave cost for frequency tremor (default 0.01)
#' @param ampTremorOctaveCost Numeric. Octave cost for amplitude tremor (default 0.01)
#' @param nanAsZero Logical. Convert undefined measurements to zeros (default FALSE)
#' @param toFile Logical. If TRUE, write results to JSTF file. Default FALSE.
#' @param explicitExt Character. File extension for output. Default "pvt".
#' @param outputDirectory Character. Output directory path. Default NULL (use input directory).
#' @param verbose Logical. Print progress messages (default TRUE)
#'
#' @return If \code{toFile=FALSE} (default), a data.frame (or list of data.frames for multiple files)
#'   with 18 tremor measurements. If \code{toFile=TRUE}, invisibly returns the path(s) to the
#'   written JSTF file(s).
#'
#'   Each result contains 18 columns:
#' \describe{
#'   \item{FCoM}{Frequency contour magnitude}
#'   \item{FTrC}{Frequency tremor cyclicality (0-1)}
#'   \item{FMoN}{Number of frequency modulation candidates}
#'   \item{FTrF}{Frequency tremor frequency (Hz)}
#'   \item{FTrI}{Frequency tremor intensity index (percent)}
#'   \item{FTrP}{Frequency tremor power index}
#'   \item{FTrCIP}{Frequency tremor cyclicality-intensity product}
#'   \item{FTrPS}{Frequency tremor product sum}
#'   \item{FCoHNR}{Frequency contour HNR (dB)}
#'   \item{ACoM}{Amplitude contour magnitude}
#'   \item{ATrC}{Amplitude tremor cyclicality (0-1)}
#'   \item{AMoN}{Number of amplitude modulation candidates}
#'   \item{ATrF}{Amplitude tremor frequency (Hz)}
#'   \item{ATrI}{Amplitude tremor intensity index (percent)}
#'   \item{ATrP}{Amplitude tremor power index}
#'   \item{ATrCIP}{Amplitude tremor cyclicality-intensity product}
#'   \item{ATrPS}{Amplitude tremor product sum}
#'   \item{ACoHNR}{Amplitude contour HNR (dB)}
#' }
#'
#' @references
#' Brückl, M. (2012). Vocal tremor measurement based on autocorrelation of contours.
#' Proceedings of Interspeech 2012, 2027-2031.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Analyze sustained vowel
#' result <- lst_voice_tremor("sustained_vowel.wav")
#' print(result$FTrF)  # Frequency tremor frequency
#' print(result$FTrI)  # Frequency tremor intensity
#'
#' # Write to JSTF file
#' lst_voice_tremor("sustained_vowel.wav", toFile = TRUE)
#' track <- read_track("sustained_vowel.pvt")
#' df <- as.data.frame(track)
#' }
lst_voice_tremor <- function(listOfFiles,
                              beginTime = 0.0,
                              endTime = 0.0,
                              analysisTimeStep = 0.015,
                              minPitch = 60,
                              maxPitch = 350,
                              silenceThreshold = 0.03,
                              voicingThreshold = 0.3,
                              octaveCost = 0.01,
                              octaveJumpCost = 0.35,
                              voicedUnvoicedCost = 0.14,
                              minTremorFreq = 1.5,
                              maxTremorFreq = 15.0,
                              tremorMagThresh = 0.01,
                              tremorCyclicalThresh = 0.15,
                              freqTremorOctaveCost = 0.01,
                              ampTremorOctaveCost = 0.01,
                              nanAsZero = FALSE,
                              toFile = FALSE,
                              explicitExt = "pvt",
                              outputDirectory = NULL,
                              verbose = TRUE) {
  
  # Check pladdrr availability
  if (!pladdrr_available()) {
    stop("pladdrr package not available. Install with: install.packages('pladdrr')")
  }
  
  # Validate files
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)
  
  filesEx <- file.exists(listOfFiles)
  if (!all(filesEx)) {
    filesNotExist <- listOfFiles[!filesEx]
    stop("Unable to find the sound file(s) ", paste(filesNotExist, collapse = ", "))
  }
  
  # Ensure time vectors match file count
  n_files <- length(listOfFiles)
  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime) == 1) endTime <- rep(endTime, n_files)
  
  # Progress bar
  if (verbose && n_files > 1) {
    pb <- txtProgressBar(min = 0, max = n_files, style = 3)
  }
  
  # Process each file
  results_list <- list()
  
  for (i in seq_along(listOfFiles)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[i]
    
    tryCatch({
      # Load sound with pladdrr
      sound <- pladdrr::Sound(file_path)
      
      # Apply time windowing if requested
      if (bt > 0 || et > 0) {
        duration <- sound$.cpp$duration
        start_time <- if (bt > 0) bt else 0
        end_time <- if (et > 0) et else duration
        
        # Extract time window
        sound <- sound$extract_part(
          from_time = start_time,
          to_time = end_time,
          window_shape = "rectangular",
          relative_width = 1.0,
          preserve_times = FALSE
        )
      }
      
      # CRITICAL: Apply Gaussian1 windowing (Praat standard for tremor analysis)
      original_duration <- sound$.cpp$duration
      sound <- sound$extract_part(
        from_time = 0,
        to_time = original_duration,
        window_shape = "Gaussian1",
        relative_width = 1.0,
        preserve_times = FALSE
      )
      
      slength <- sound$.cpp$duration
      
      # Analyze frequency tremor
      ftrem <- analyze_frequency_tremor_pladdrr(
        sound = sound,
        slength = slength,
        analysisTimeStep = analysisTimeStep,
        minPitch = minPitch,
        maxPitch = maxPitch,
        silenceThreshold = silenceThreshold,
        voicingThreshold = voicingThreshold,
        octaveCost = octaveCost,
        octaveJumpCost = octaveJumpCost,
        voicedUnvoicedCost = voicedUnvoicedCost,
        minTremorFreq = minTremorFreq,
        maxTremorFreq = maxTremorFreq,
        tremorMagThresh = tremorMagThresh,
        tremorCyclicalThresh = tremorCyclicalThresh,
        freqTremorOctaveCost = freqTremorOctaveCost,
        nanAsZero = nanAsZero
      )
      
      # Analyze amplitude tremor
      atrem <- analyze_amplitude_tremor_pladdrr(
        sound = sound,
        slength = slength,
        analysisTimeStep = analysisTimeStep,
        minPitch = minPitch,
        silenceThreshold = silenceThreshold,
        voicingThreshold = voicingThreshold,
        octaveCost = octaveCost,
        octaveJumpCost = octaveJumpCost,
        voicedUnvoicedCost = voicedUnvoicedCost,
        minTremorFreq = minTremorFreq,
        maxTremorFreq = maxTremorFreq,
        tremorMagThresh = tremorMagThresh,
        tremorCyclicalThresh = tremorCyclicalThresh,
        ampTremorOctaveCost = ampTremorOctaveCost,
        nanAsZero = nanAsZero
      )
      
      # Combine results
      result_df <- data.frame(
        file = basename(file_path),
        FCoM = ftrem$FCoM,
        FTrC = ftrem$FTrC,
        FMoN = ftrem$FMoN,
        FTrF = ftrem$FTrF,
        FTrI = ftrem$FTrI,
        FTrP = ftrem$FTrP,
        FTrCIP = ftrem$FTrCIP,
        FTrPS = ftrem$FTrPS,
        FCoHNR = ftrem$FCoHNR,
        ACoM = atrem$ACoM,
        ATrC = atrem$ATrC,
        AMoN = atrem$AMoN,
        ATrF = atrem$ATrF,
        ATrI = atrem$ATrI,
        ATrP = atrem$ATrP,
        ATrCIP = atrem$ATrCIP,
        ATrPS = atrem$ATrPS,
        ACoHNR = atrem$ACoHNR,
        stringsAsFactors = FALSE
      )
      
      results_list[[i]] <- result_df
      
    }, error = function(e) {
      warning("Error processing file ", file_path, ": ", conditionMessage(e))
      # Return NA results
      result_df <- data.frame(
        file = basename(file_path),
        FCoM = NA_real_, FTrC = NA_real_, FMoN = NA_integer_,
        FTrF = NA_real_, FTrI = NA_real_, FTrP = NA_real_,
        FTrCIP = NA_real_, FTrPS = NA_real_, FCoHNR = NA_real_,
        ACoM = NA_real_, ATrC = NA_real_, AMoN = NA_integer_,
        ATrF = NA_real_, ATrI = NA_real_, ATrP = NA_real_,
        ATrCIP = NA_real_, ATrPS = NA_real_, ACoHNR = NA_real_,
        stringsAsFactors = FALSE
      )
      results_list[[i]] <- result_df
    })
    
    if (verbose && n_files > 1) {
      setTxtProgressBar(pb, i)
    }
  }
  
  if (verbose && n_files > 1) {
    close(pb)
  }
  
  # Handle JSTF file writing
  if (toFile) {
    # Convert each data.frame row to named list
    results_for_jstf <- lapply(results_list, function(df) {
      df$file <- NULL  # Remove file column
      as.list(df[1, ])  # Convert to named list
    })
    
    output_paths <- write_lst_results_to_jstf(
      results = results_for_jstf,
      file_paths = listOfFiles,
      beginTime = beginTime,
      endTime = endTime,
      function_name = "lst_voice_tremor",
      parameters = list(
        analysisTimeStep = analysisTimeStep,
        minPitch = minPitch,
        maxPitch = maxPitch,
        silenceThreshold = silenceThreshold,
        voicingThreshold = voicingThreshold,
        octaveCost = octaveCost,
        octaveJumpCost = octaveJumpCost,
        voicedUnvoicedCost = voicedUnvoicedCost,
        minTremorFreq = minTremorFreq,
        maxTremorFreq = maxTremorFreq,
        tremorMagThresh = tremorMagThresh,
        tremorCyclicalThresh = tremorCyclicalThresh,
        freqTremorOctaveCost = freqTremorOctaveCost,
        ampTremorOctaveCost = ampTremorOctaveCost,
        nanAsZero = nanAsZero
      ),
      explicitExt = explicitExt,
      outputDirectory = outputDirectory
    )
    
    return(invisible(output_paths))
  }
  
  # Combine all results into single data.frame
  results_df <- do.call(rbind, results_list)
  rownames(results_df) <- NULL
  
  return(results_df)
}


#' Analyze frequency tremor using pladdrr
#' 
#' @keywords internal
analyze_frequency_tremor_pladdrr <- function(sound, slength, analysisTimeStep,
                                             minPitch, maxPitch,
                                             silenceThreshold, voicingThreshold,
                                             octaveCost, octaveJumpCost, voicedUnvoicedCost,
                                             minTremorFreq, maxTremorFreq,
                                             tremorMagThresh, tremorCyclicalThresh,
                                             freqTremorOctaveCost,
                                             nanAsZero) {
  
  tryCatch({
    # Step 1: Extract pitch with Direct API
    pitch_ptr <- pladdrr::to_pitch_cc_direct(
      sound,
      time_step = analysisTimeStep,
      pitch_floor = minPitch,
      max_candidates = 15,
      very_accurate = FALSE,
      silence_threshold = silenceThreshold,
      voicing_threshold = voicingThreshold,
      octave_cost = octaveCost,
      octave_jump_cost = octaveJumpCost,
      voiced_unvoiced_cost = voicedUnvoicedCost,
      pitch_ceiling = maxPitch
    )
    pitch <- pladdrr::Pitch(.xptr = pitch_ptr)
    
    # Step 2: Extract F0 values using vectorized methods
    f0_raw <- pitch$get_values_vector()
    voiced_mask <- pitch$get_voiced_mask()
    num_frames <- length(f0_raw)
    number_voiced <- sum(voiced_mask)
    
    if (number_voiced == 0) {
      return(.undefined_ftrem_results(nanAsZero))
    }
    
    # Replace NA/<=0 with 0
    f0_list <- ifelse(is.na(f0_raw) | f0_raw <= 0, 0, f0_raw)
    
    if (sum(f0_list > 0) < 10) {
      return(.undefined_ftrem_results(nanAsZero))
    }
    
    # Calculate mean F0 from voiced frames
    mean_f0 <- mean(f0_list[f0_list > 0])
    
    # Step 3: Detrend using native Praat detrending (v4.0.14+) or R fallback
    if (.has_pitch_detrend()) {
      # FAST: Native Praat detrending
      f0_detrended <- pitch$get_values_detrended(unit = "hertz")
      f0_detrended <- ifelse(is.na(f0_detrended), 0, f0_detrended)
    } else {
      # FALLBACK: R linear regression
      voiced_indices <- which(f0_list > 0)
      lm_fit <- lm(f0_list[voiced_indices] ~ voiced_indices)
      all_indices <- seq_len(num_frames)
      predicted <- predict(lm_fit, newdata = data.frame(voiced_indices = all_indices))
      f0_detrended <- ifelse(f0_list > 0, f0_list - predicted, 0)
    }
    
    # Step 4: Normalize by mean F0
    f0_normalized <- f0_detrended / mean_f0
    
    # Step 5: Create Sound with exact sampling rate
    sampling_rate <- 1.0 / analysisTimeStep
    f0_sound <- pladdrr::Sound$from_values(f0_normalized, sampling_rate = sampling_rate)
    f0_duration <- f0_sound$.cpp$duration
    
    # Step 6: Apply Pitch to detect tremor (single-frame analysis)
    tremor_pitch_ptr <- pladdrr::to_pitch_cc_direct(
      f0_sound,
      time_step = f0_duration,
      pitch_floor = minTremorFreq,
      max_candidates = 15,
      very_accurate = TRUE,
      silence_threshold = tremorMagThresh,
      voicing_threshold = tremorCyclicalThresh,
      octave_cost = freqTremorOctaveCost,
      octave_jump_cost = 0.35,
      voiced_unvoiced_cost = 0.14,
      pitch_ceiling = maxTremorFreq
    )
    tremor_pitch <- pladdrr::Pitch(.xptr = tremor_pitch_ptr)
    
    # Step 7: Extract tremor frequency at center
    ftrf <- tremor_pitch$get_value_at_time(f0_duration / 2)
    
    if (is.na(ftrf) || ftrf == 0) {
      return(.undefined_ftrem_results(nanAsZero))
    }
    
    # Step 8: Extract intensity and strength
    pitch_info <- .extract_pitch_intensity_and_strength(tremor_pitch, maxTremorFreq, slength)
    ftrm <- pitch_info$intensity
    ftrc <- pitch_info$strength
    
    # Step 9: Calculate tremor intensity index
    ftri <- .calculate_tremor_intensity(f0_sound, tremor_pitch, f0_duration, apply_correction = TRUE)
    
    # Step 10: Calculate HNR
    ftr_hnr <- .calculate_contour_hnr(f0_sound, minTremorFreq, maxTremorFreq)
    
    # Calculate derived measures
    ftrp <- if (ftrf > 0) ftri * ftrf / (ftrf + 1) else 0.0
    ftrcip <- ftri * ftrc
    ftrps <- ftrcip
    fmon <- as.integer(slength * ftrf)
    
    return(list(
      FCoM = ftrm,
      FTrC = ftrc,
      FMoN = fmon,
      FTrF = ftrf,
      FTrI = ftri,
      FTrP = ftrp,
      FTrCIP = ftrcip,
      FTrPS = ftrps,
      FCoHNR = ftr_hnr
    ))
    
  }, error = function(e) {
    warning("Error in frequency tremor analysis: ", conditionMessage(e))
    return(.undefined_ftrem_results(nanAsZero))
  })
}


#' Analyze amplitude tremor using pladdrr
#' 
#' @keywords internal
analyze_amplitude_tremor_pladdrr <- function(sound, slength, analysisTimeStep,
                                             minPitch,
                                             silenceThreshold, voicingThreshold,
                                             octaveCost, octaveJumpCost, voicedUnvoicedCost,
                                             minTremorFreq, maxTremorFreq,
                                             tremorMagThresh, tremorCyclicalThresh,
                                             ampTremorOctaveCost,
                                             nanAsZero) {
  
  tryCatch({
    # Step 1: Extract amplitude contour using Intensity
    intensity <- sound$to_intensity(minimum_pitch = minPitch)
    
    # Step 2: Extract amplitude values for all frames
    intensity_df <- as.data.frame(intensity)
    amp_db_raw <- intensity_df$intensity_db
    
    # Convert dB to linear, NA becomes 0
    amp_list <- ifelse(is.na(amp_db_raw), 0, 10^(amp_db_raw / 20.0))
    
    if (sum(amp_list > 0) < 10) {
      return(.undefined_atrem_results(nanAsZero))
    }
    
    # Step 3: Normalize amplitude
    mean_amp <- mean(amp_list[amp_list > 0])
    amp_normalized <- (amp_list - mean_amp) / mean_amp
    
    # Step 4: Create Sound with exact sampling rate
    sampling_rate <- 1.0 / analysisTimeStep
    amp_sound <- pladdrr::Sound$from_values(amp_normalized, sampling_rate = sampling_rate)
    amp_duration <- amp_sound$.cpp$duration
    
    # Step 5: Apply Pitch to detect tremor
    tremor_pitch_ptr <- pladdrr::to_pitch_cc_direct(
      amp_sound,
      time_step = amp_duration,
      pitch_floor = minTremorFreq,
      max_candidates = 15,
      very_accurate = TRUE,
      silence_threshold = tremorMagThresh,
      voicing_threshold = tremorCyclicalThresh,
      octave_cost = ampTremorOctaveCost,
      octave_jump_cost = 0.35,
      voiced_unvoiced_cost = 0.14,
      pitch_ceiling = maxTremorFreq
    )
    tremor_pitch <- pladdrr::Pitch(.xptr = tremor_pitch_ptr)
    
    # Step 6: Extract tremor frequency at center
    atrf <- tremor_pitch$get_value_at_time(amp_duration / 2)
    
    # Step 7: Extract intensity and strength (BEFORE checking atrf)
    pitch_info <- .extract_pitch_intensity_and_strength(tremor_pitch, maxTremorFreq, amp_duration)
    atrm <- pitch_info$intensity
    atrc <- pitch_info$strength
    
    if (is.na(atrf) || atrf == 0) {
      # No tremor detected, but return ACoM
      val <- if (nanAsZero) 0.0 else NA_real_
      atr_hnr <- .calculate_contour_hnr(amp_sound, minTremorFreq, maxTremorFreq)
      return(list(
        ACoM = atrm,
        ATrC = val,
        AMoN = 0L,
        ATrF = val,
        ATrI = val,
        ATrP = val,
        ATrCIP = val,
        ATrPS = val,
        ACoHNR = atr_hnr
      ))
    }
    
    # Step 8: Calculate tremor intensity index
    atri <- .calculate_tremor_intensity(amp_sound, tremor_pitch, amp_duration, apply_correction = FALSE)
    
    # Step 9: Calculate HNR
    atr_hnr <- .calculate_contour_hnr(amp_sound, minTremorFreq, maxTremorFreq)
    
    # Calculate derived measures
    atrp <- if (atrf > 0) atri * atrf / (atrf + 1) else 0.0
    atrcip <- atri * atrc
    atrps <- atrcip
    amon <- as.integer(slength * atrf)
    
    return(list(
      ACoM = atrm,
      ATrC = atrc,
      AMoN = amon,
      ATrF = atrf,
      ATrI = atri,
      ATrP = atrp,
      ATrCIP = atrcip,
      ATrPS = atrps,
      ACoHNR = atr_hnr
    ))
    
  }, error = function(e) {
    warning("Error in amplitude tremor analysis: ", conditionMessage(e))
    return(.undefined_atrem_results(nanAsZero))
  })
}


#' Extract pitch intensity and strength from Pitch object
#' 
#' @keywords internal
.extract_pitch_intensity_and_strength <- function(pitch_object, max_freq, slength = NULL) {
  tryCatch({
    # Get strength at center time or use mean
    if (!is.null(slength)) {
      strength <- pitch_object$get_strength_at_time(slength / 2)
    } else {
      strength <- pitch_object$get_mean_strength()
    }
    
    # If strength is NA or invalid, use mean as fallback
    if (is.na(strength) || strength < 0) {
      strength <- pitch_object$get_mean_strength()
      if (is.na(strength) || strength < 0) {
        strength <- 0.0
      }
    }
    
    # For intensity, use mean strength as proxy
    intensity <- strength
    
    return(list(intensity = intensity, strength = strength))
    
  }, error = function(e) {
    warning("Error extracting pitch strength: ", conditionMessage(e))
    return(list(intensity = 0.0, strength = 0.0))
  })
}


#' Calculate tremor intensity following Praat's tremIntIndex procedure
#' 
#' @keywords internal
.calculate_tremor_intensity <- function(contour_sound, tremor_pitch, slength, apply_correction = FALSE) {
  tryCatch({
    # Create PointProcess of peaks (maxima)
    pp_max <- contour_sound$pitch_to_pointprocess_peaks(
      tremor_pitch,
      include_maxima = TRUE,
      include_minima = FALSE
    )
    
    num_max_points <- pp_max$get_number_of_points()
    
    # Calculate peak contribution
    tri_max <- 0.0
    no_f_max <- 0
    
    if (num_max_points > 0) {
      if (.has_pp_values_from_sound()) {
        # FAST: Single call (v4.0.14+)
        peak_values <- pp_max$get_values_from_sound(
          contour_sound,
          channel = 1L,
          interpolation = "sinc70"
        )
      } else {
        # FALLBACK: Two calls
        peak_times <- pp_max$get_pointprocess_times()
        peak_values <- contour_sound$get_values_at_times(
          peak_times,
          channel = 1L,
          interpolation = "sinc70"
        )
      }
      
      # Vectorized sum and count
      valid_mask <- !is.na(peak_values)
      no_f_max <- sum(!valid_mask)
      tri_max <- sum(abs(peak_values[valid_mask]))
    }
    
    # Calculate percentage, excluding undefined points
    number_of_maxima <- num_max_points - no_f_max
    if (number_of_maxima > 0) {
      tri_max <- 100 * tri_max / number_of_maxima
    } else {
      tri_max <- 0.0
    }
    
    # Create PointProcess of troughs (minima)
    pp_min <- contour_sound$pitch_to_pointprocess_peaks(
      tremor_pitch,
      include_maxima = FALSE,
      include_minima = TRUE
    )
    
    num_min_points <- pp_min$get_number_of_points()
    
    tri_min <- 0.0
    no_f_min <- 0
    
    if (num_min_points > 0) {
      if (.has_pp_values_from_sound()) {
        # FAST: Single call
        trough_values <- pp_min$get_values_from_sound(
          contour_sound,
          channel = 1L,
          interpolation = "sinc70"
        )
      } else {
        # FALLBACK: Two calls
        trough_times <- pp_min$get_pointprocess_times()
        trough_values <- contour_sound$get_values_at_times(
          trough_times,
          channel = 1L,
          interpolation = "sinc70"
        )
      }
      
      # Vectorized sum and count
      valid_mask <- !is.na(trough_values)
      no_f_min <- sum(!valid_mask)
      tri_min <- sum(abs(trough_values[valid_mask]))
    }
    
    number_of_minima <- num_min_points - no_f_min
    if (number_of_minima > 0) {
      tri_min <- 100 * tri_min / number_of_minima
    } else {
      tri_min <- 0.0
    }
    
    # Average of peaks and troughs
    tri <- (tri_max + tri_min) / 2
    
    # Apply correction factor for frequency tremor
    if (apply_correction) {
      tri <- tri * 0.1
    }
    
    return(tri)
    
  }, error = function(e) {
    warning("Error calculating tremor intensity: ", conditionMessage(e))
    return(0.0)
  })
}


#' Calculate harmonicity (HNR) for contour sound
#' 
#' @keywords internal
.calculate_contour_hnr <- function(sound, min_tremor_freq, max_tremor_freq) {
  tryCatch({
    # Time step for harmonicity analysis
    time_step <- 1.0 / max_tremor_freq / 4.0
    
    # Calculate harmonicity
    harmonicity <- sound$to_harmonicity_ac(
      time_step = time_step,
      min_pitch = min_tremor_freq,
      silence_threshold = 0.1,
      periods_per_window = 1.0
    )
    
    # Get mean HNR
    hnr <- harmonicity$get_mean(from_time = 0, to_time = 0)
    
    if (is.na(hnr) || is.infinite(hnr)) {
      return(0.0)
    }
    
    return(hnr)
    
  }, error = function(e) {
    # Known limitation: HNR fails on short tremor contours
    return(0.0)
  })
}


#' Check if pladdrr v4.0.14+ detrending is available
#' 
#' @keywords internal
.has_pitch_detrend <- function() {
  tryCatch({
    pkg_version <- packageVersion("pladdrr")
    pkg_version >= "4.0.14"
  }, error = function(e) FALSE)
}


#' Check if pladdrr v4.0.14+ pp$get_values_from_sound is available
#' 
#' @keywords internal
.has_pp_values_from_sound <- function() {
  tryCatch({
    pkg_version <- packageVersion("pladdrr")
    pkg_version >= "4.0.14"
  }, error = function(e) FALSE)
}


#' Return undefined frequency tremor results
#' 
#' @keywords internal
.undefined_ftrem_results <- function(nanAsZero) {
  val <- if (nanAsZero) 0.0 else NA_real_
  return(list(
    FCoM = val, FTrC = val, FMoN = 0L, FTrF = val,
    FTrI = val, FTrP = val, FTrCIP = val, FTrPS = val, FCoHNR = val
  ))
}


#' Return undefined amplitude tremor results
#' 
#' @keywords internal
.undefined_atrem_results <- function(nanAsZero) {
  val <- if (nanAsZero) 0.0 else NA_real_
  return(list(
    ACoM = val, ATrC = val, AMoN = 0L, ATrF = val,
    ATrI = val, ATrP = val, ATrCIP = val, ATrPS = val, ACoHNR = val
  ))
}


# Set function attributes
attr(lst_voice_tremor, "ext") <- "pvt"
attr(lst_voice_tremor, "outputType") <- "JSTF"
attr(lst_voice_tremor, "format") <- "JSON"
attr(lst_voice_tremor, "tracks") <- c(
  "FCoM", "FTrC", "FMoN", "FTrF", "FTrI",
  "FTrP", "FTrCIP", "FTrPS", "FCoHNR",
  "ACoM", "ATrC", "AMoN", "ATrF", "ATrI",
  "ATrP", "ATrCIP", "ATrPS", "ACoHNR"
)

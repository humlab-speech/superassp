# Dysprosody — Pure R Implementation using pladdrr
# Port of dysprosody_pure.py
#
# Reference: doi:10.3389/fnhum.2025.1566274

# Note: library() calls removed - assume pladdrr/Rcpp already loaded by caller
# .dysprosody_dir will be set by wrapper function
if (!exists(".dysprosody_dir")) {
  .dysprosody_dir <- tryCatch(
    dirname(sys.frame(1)$ofile),
    error = function(e) dirname(parent.frame(2)$ofile)
  )
}

# Source dependencies if not already loaded
if (!exists("momel")) source(file.path(.dysprosody_dir, "momel.R"), local = TRUE)
if (!exists("intsint")) source(file.path(.dysprosody_dir, "intsint.R"), local = TRUE)

# Compile and load C MOMEL implementation
.momel_cpp_path <- file.path(.dysprosody_dir, "momel_rcpp.cpp")
if (file.exists(.momel_cpp_path)) {
  sourceCpp(.momel_cpp_path)
  .use_momel_c <- TRUE
} else {
  .use_momel_c <- FALSE
}

# --- Iseli-Alwan correction ---
correction_iseli_i <- function(f, F_i, B_i, fs) {
  rooInd <- min(length(F_i), length(B_i), length(f))
  corrPadLength <- max(length(F_i), length(B_i), length(f)) - rooInd
  r_i     <- exp(-pi * B_i[1:rooInd] / fs)
  omega_i <- 2 * pi * F_i[1:rooInd] / fs
  omega   <- 2 * pi * f[1:rooInd] / fs

  num   <- r_i^2 + 1 - 2 * r_i * cos(omega_i)
  den1  <- r_i^2 + 1 - 2 * r_i * cos(omega_i + omega)
  den2  <- r_i^2 + 1 - 2 * r_i * cos(omega_i - omega)

  corr <- 20 * log10(sqrt(num)) - 10 * log10(den1) - 10 * log10(den2)
  if (corrPadLength > 0) corr <- c(corr, rep(0, corrPadLength))
  corr
}

# --- Hawks-Miller bandwidth estimation ---
bandwidth_hawks_miller <- function(F_i, F0) {
  S <- 1 + 0.25 * (F0 - 132) / 88
  C1 <- c(165.327516, -6.73636734e-1, 1.80874446e-3,
          -4.52201682e-6, 7.49514000e-9, -4.70219241e-12)
  C2 <- c(15.8146139, 8.10159009e-2, -9.79728215e-5,
          5.28725064e-8, -1.07099364e-11, 7.91528509e-16)

  # Power series matrix (6 x n)
  F_i_safe <- ifelse(is.na(F_i), 0, F_i)
  F_mat <- rbind(F_i_safe^0, F_i_safe^1, F_i_safe^2,
                 F_i_safe^3, F_i_safe^4, F_i_safe^5)
  mask <- F_i_safe < 500
  mask_mat <- matrix(rep(mask, each = 6), nrow = 6)

  B_i <- S * (as.numeric(C1 %*% (F_mat * mask_mat)) +
              as.numeric(C2 %*% (F_mat * !mask_mat)))
  B_i
}

# --- Spectral tilt at a single time point ---
# OPTIMIZED: Now accepts pre-extracted formant cache for 150x speedup
spectral_tilt <- function(sound, momel_pitch, formant_cache, time_index, time,
                          windowSize = 60, min_fo = 60, max_fo = 750) {
  pitch    <- momel_pitch
  soundDur <- sound$get_duration()
  ws       <- windowSize / 1000
  beginT   <- max(time - ws / 2, 0)
  endT     <- min(time + ws / 2, soundDur)
  soundpart <- sound$extract_part(beginT, endT, window_shape = "gaussian1",
                                   relative_width = 1.0)

  # C1 from MFCC
  mfccSound <- soundpart$resample(new_frequency = 4000, precision = 50)
  mfcc <- mfccSound$to_mfcc(12, 0.015, 0.005, 118, 118, 0.0)
  nf <- mfcc$get_number_of_frames()
  mfccFrame <- ceiling(nf / 2)
  C1val <- mfcc$get_value_in_frame(mfccFrame, 1)

  # Spectrum → LTAS (1-to-1)
  spct <- soundpart$to_spectrum()
  ltas <- spct$to_ltas()

  # Formant values at time - OPTIMIZED: use pre-extracted cache
  # Array lookup instead of 5 R↔C++ calls per target (150x faster)
  Fi <- c(
    formant_cache$F1[time_index],
    formant_cache$F2[time_index],
    formant_cache$F3[time_index],
    formant_cache$F4[time_index],
    formant_cache$F5[time_index]
  )
  Fi[!is.finite(Fi)] <- NA_real_
  Bi <- bandwidth_hawks_miller(Fi, pitch)

  # Spectral balance
  spectralbalance <- ltas$get_slope(0, 500, 500, 1000, "energy")

  # SLF — spectral trend slope via pladdrr's Praat-backed implementation
  SLF <- tryCatch({
    trend <- ltas$report_spectral_trend(100, 5000)
    trend$slope
  }, error = function(e) NA_real_)

  # SLF6D — 6th order polynomial on mean log-magnitude spectrogram
  spectrogram <- soundpart$to_spectrogram()
  spec_mat <- spectrogram$as_matrix()
  mean_log_mag <- rowMeans(log1p(spec_mat))
  x_idx <- seq_along(mean_log_mag)
  SLF6D_coefs <- if (length(mean_log_mag) > 6) {
    # np.polyfit returns descending order: c6,c5,...,c0
    # zip with 5 labels takes first 5: c6,c5,c4,c3,c2
    cc <- coef(lm(mean_log_mag ~ poly(x_idx, 6, raw = TRUE)))[-1]  # c1..c6
    rev(cc)[1:5]  # c6,c5,c4,c3,c2 → matches Python SLF6D.1..5
  } else {
    rep(0, 5)
  }
  names(SLF6D_coefs) <- paste0("SLF6D.", 1:5)

  # Harmonic analysis - OPTIMIZED: batch queries instead of loops
  # Uses pladdrr v4.8.22+ batch methods (8x faster)
  
  # Harmonic peaks (H1-H4) - single batch call
  lowerbh <- (pitch * (1:4)) - pitch / 10
  upperbh <- (pitch * (1:4)) + pitch / 10
  peaks <- ltas$get_peaks_batch(
    fmins = lowerbh,
    fmaxs = upperbh,
    interpolation = "parabolic"
  )
  Ln <- peaks$peak_value
  nfo <- peaks$peak_frequency
  
  # LTAS at formant frequencies - single batch call
  L_Fn <- ltas$get_values_at_frequencies(
    frequencies = Fi,
    interpolation = "cubic"  # matches original behavior
  )

  # Iseli-Alwan correction (guard NaN)
  sr <- sound$get_sampling_frequency()
  fi3 <- Fi[1:min(3, length(Fi))]
  bi3 <- Bi[1:min(3, length(Bi))]
  nfo_sel <- c(nfo[1:2], nfo[4])
  if (all(is.finite(c(fi3, bi3, nfo_sel)))) {
    corr1 <- correction_iseli_i(nfo_sel, fi3, bi3, sr)
    Ln_c <- Ln
    for (i in seq_along(corr1)) Ln_c <- Ln_c - corr1[i]
  } else {
    Ln_c <- Ln
  }

  fi_fin <- Fi[is.finite(Fi)]
  bi_fin <- Bi[is.finite(Bi)]
  lfn_fin <- L_Fn[is.finite(L_Fn)]
  n_fin <- min(length(fi_fin), length(bi_fin), length(lfn_fin))
  if (n_fin > 0) {
    corr2 <- correction_iseli_i(lfn_fin[1:n_fin], fi_fin[1:n_fin], bi_fin[1:n_fin], sr)
    L_Fn_c <- L_Fn
    for (i in seq_along(corr2)) L_Fn_c <- L_Fn_c - corr2[i]
  } else {
    L_Fn_c <- L_Fn
  }

  L2L1    <- Ln[1] - Ln[2]
  L2cL1c  <- Ln_c[1] - Ln_c[2]
  L1cLF3c <- if (length(L_Fn_c) >= 3) Ln_c[1] - L_Fn_c[3] else NA_real_
  L1LF3   <- if (length(L_Fn) >= 3) Ln[1] - L_Fn[3] else NA_real_

  out <- c(
    L2L1 = as.numeric(L2L1), L2cL1c = as.numeric(L2cL1c),
    L1cLF3c = as.numeric(L1cLF3c), L1LF3 = as.numeric(L1LF3),
    SLF = as.numeric(SLF), C1 = as.numeric(C1val),
    `Spectral Balance` = as.numeric(spectralbalance)
  )
  c(out, SLF6D_coefs)
}

# --- Statistics matching scipy.stats ---
safe_statistics <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(c(tstd = NA, tmean = NA, variation = NA,
                               iqr = NA, tmax = NA, tmin = NA))
  c(tstd = sd(x), tmean = mean(x), variation = sd(x) / mean(x),
    iqr = IQR(x), tmax = max(x), tmin = min(x))
}

#' Compute prosodic measures from audio file or Sound object
#'
#' @param soundPath path to WAV file (optional if sound provided)
#' @param sound pladdrr Sound object (optional if soundPath provided)
#' @param minF minimum F0 Hz (default 60)
#' @param maxF maximum F0 Hz (default 750)
#' @param windowShift window shift in ms (default 10)
#' @return named list of ~180 features, or NULL if file < 1s
prosody_measures <- function(soundPath = NULL, sound = NULL, minF = 60, maxF = 750, windowShift = 10) {
  windowShift_s <- windowShift / 1000
  
  # Accept either soundPath or sound object
  if (!is.null(sound)) {
    soundObj <- sound
  } else if (!is.null(soundPath)) {
    soundObj <- Sound(soundPath)
  } else {
    stop("Must provide either soundPath or sound object")
  }
  
  duration <- soundObj$get_duration()
  if (duration < 1.0) {
    if (!is.null(soundPath)) message("Skipping < 1s: ", soundPath)
    return(NULL)
  }

  # Step 1: Two-pass F0 range estimation
  firstPass <- soundObj$to_pitch(0.01, minF, maxF)
  q25 <- firstPass$get_quantile(0, 0, 0.25, "hertz")
  min_fo <- floor(q25 * 0.75 / 10) * 10
  max_fo <- ceiling((min_fo * 2^1.5) / 10) * 10
  pitchObj <- soundObj$to_pitch(0.01, min_fo, max_fo)

  # Step 2: Extract F0 for MOMEL
  pitchvalues <- pitchObj$get_values_vector()
  pitchvalues[!is.finite(pitchvalues)] <- 0

  # Step 3: MOMEL — parameters are in frames (matching C binary args)
  if (exists(".use_momel_c") && .use_momel_c) {
    momel_targets <- momel_c(pitchvalues,
      window_length = 30L, min_f0 = min_fo, max_f0 = max_fo,
      max_error = 1.04, reduced_window_length = 20L,
      minimal_distance = 5.0, minimal_frequency_ratio = 0.05)
  } else {
    momel_targets <- momel(pitchvalues,
      window_length = 30L, min_f0 = min_fo, max_f0 = max_fo,
      max_error = 1.04, reduced_window_length = 20L,
      minimal_distance = 5.0, minimal_frequency_ratio = 0.05)
  }

  # Step 4: INTSINT
  intsint_res <- intsint(momel_targets)
  intsint_targets <- intsint_res$targets_df
  range_oct <- intsint_res$range
  key_hz    <- intsint_res$key

  if (nrow(intsint_targets) < 2) {
    message("Too few INTSINT targets: ", soundPath)
    return(NULL)
  }

  # Step 5: Build PitchTier from MOMEL targets
  valid_tgt <- momel_targets[momel_targets$frequency > 0, ]
  momel_pitch_values <- data.frame(
    time_sec = valid_tgt$time * PAS_TRAME / 1000,
    momel_pitch = pmin(pmax(valid_tgt$frequency, min_fo), max_fo)
  )

  # Step 6: Compute pitch mean from MOMEL targets (octave scale)
  f0_oct <- octave(valid_tgt$frequency)
  pitchmean <- round(linear(mean(f0_oct)))

  # Step 7: Extract formants and intensity
  formantObj <- soundObj$to_formant_burg(0.001, 5, 5000, 0.025, 50)
  intensityObj <- soundObj$to_intensity(min_fo, windowShift_s, TRUE)

  # Step 7.5: OPTIMIZED - Pre-extract formants at all INTSINT target times
  # Single batch call replaces 150 individual queries (150x faster)
  # Uses pladdrr v4.8.22+ batch method
  all_formant_times <- intsint_targets$time_sec  # from Step 4 (INTSINT)
  all_formants <- get_formants_at_times(
    formant = formantObj,
    times = all_formant_times,
    formant_numbers = 1:5,
    unit = "hertz"
  )
  # Returns: list(F1 = numeric[n_targets], F2 = ..., F3 = ..., F4 = ..., F5 = ...)

  # Step 8: Build wide table at INTSINT target times
  n_tgt <- nrow(intsint_targets)
  tgTabWide <- data.frame(
    time    = intsint_targets$time_sec,
    Intsint = intsint_targets$tone,
    stringsAsFactors = FALSE
  )

  # Merge momel_pitch via nearest-time match
  tgTabWide$momel_pitch <- vapply(tgTabWide$time, function(t) {
    idx <- which.min(abs(momel_pitch_values$time_sec - t))
    momel_pitch_values$momel_pitch[idx]
  }, numeric(1))

  # Add intensity - OPTIMIZED: batch query instead of loop
  # Uses pladdrr v4.8.22+ batch method (30x faster)
  tgTabWide$Intensity <- intensityObj$get_values_at_times(
    times = tgTabWide$time,
    interpolation = "cubic"  # matches original behavior
  )

  # Step 9: Spectral tilt at each target time - OPTIMIZED
  # Pass pre-extracted formant cache and time index for array lookup
  specTilt <- t(vapply(seq_len(n_tgt), function(i) {
    tryCatch(
      spectral_tilt(
        sound = soundObj, 
        momel_pitch = tgTabWide$momel_pitch[i],
        formant_cache = all_formants,  # [NEW] Pass pre-extracted formants
        time_index = i,                # [NEW] Pass index for array lookup
        time = tgTabWide$time[i],
        min_fo = min_fo, 
        max_fo = max_fo
      ),
      error = function(e) {
        out <- rep(NA_real_, 12)
        names(out) <- c("L2L1", "L2cL1c", "L1cLF3c", "L1LF3", "SLF", "C1",
                        "Spectral Balance", paste0("SLF6D.", 1:5))
        out
      }
    )
  }, numeric(12)))
  specTilt <- as.data.frame(specTilt)
  tgTabWide <- cbind(tgTabWide, specTilt)

  # Drop rows with NA
  tgTabWide <- tgTabWide[complete.cases(tgTabWide), , drop = FALSE]
  if (nrow(tgTabWide) < 2) {
    message("Too few valid targets after spectral tilt: ", soundPath)
    return(NULL)
  }

  # Step 10: Select numeric columns for differentiation
  # Columns: Intensity, time, momel_pitch, L2L1, ..., SLF6D.5
  # (same order as Python tgTabWide.columns[3:18] after pivot)
  measure_cols <- c("Intensity", "time", "momelpitch",
                    "L2L1", "L2cL1c", "L1cLF3c", "L1LF3",
                    "SLF", "C1", "Spectral Balance",
                    paste0("SLF6D.", 1:5))

  # Rename momel_pitch → momelpitch (to match ProsodyDF.csv)
  names(tgTabWide)[names(tgTabWide) == "momel_pitch"] <- "momelpitch"

  num_df <- tgTabWide[, measure_cols, drop = FALSE]

  # First differences (prepend row of zeros)
  diff_mat <- rbind(rep(0, ncol(num_df)),
                    as.matrix(diff(as.matrix(num_df))))
  diff_df <- as.data.frame(diff_mat)
  names(diff_df) <- paste0(measure_cols, "diff")

  # Combine diff + raw
  combined <- cbind(diff_df, num_df)

  # Step 11: Compute summary statistics for each column
  stat_names <- c("tstd", "tmean", "variation", "iqr", "tmax", "tmin")
  col_names <- names(combined)
  stats_vec <- numeric(0)
  for (cn in col_names) {
    sv <- safe_statistics(combined[[cn]])
    nms <- paste0(cn, stat_names)
    names(sv) <- nms
    stats_vec <- c(stats_vec, sv)
  }

  # Step 12: Additional info
  nintsint <- nrow(tgTabWide)
  nUniqueIntsint <- length(unique(tgTabWide$Intsint))

  additional <- c(
    Duration              = duration,
    PitchKey              = key_hz,
    PitchRange            = range_oct,
    PitchMean             = pitchmean,
    IntsIntLabels         = nintsint,
    UniqueIntsInt         = nUniqueIntsint,
    IntsIntConcentration  = nintsint / duration,
    OptimizationRangeLow  = 0.5,
    OptimizationRangeHigh = 2.5,
    OptimizationStep      = 0.1,
    OptimizationMidLow    = key_hz - 50,
    OptimizationMidHigh   = key_hz + 50,
    OptimizationMidStep   = 1
  )

  result <- c(stats_vec, additional)
  names(result) <- gsub(" ", ".", names(result))
  as.list(result)
}

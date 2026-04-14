#' SEDREAMS Glottal Closure Instant Detection
#'
#' Detects glottal closure instants (GCIs) using SEDREAMS algorithm.
#' Returns event times (GCI instants), not a regular frame grid.
#'
#' @param listOfFiles Vector of file paths (WAV, MP3, MP4, etc.) to analyze
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param f0mean Estimated mean F0 in Hz. If NULL, auto-estimated from signal.
#' @param polarity Signal polarity (1 or -1). If NULL, auto-detected.
#' @param verbose Show progress messages (default: TRUE)
#'
#' @return
#' Data frame with columns:
#' - `file`: Input file path
#' - `n_gcis`: Number of detected GCIs
#' - `gci_times`: List column of numeric vectors (GCI times in seconds)
#'
#' @details
#' **SEDREAMS Algorithm**:
#' 1. Compute LPC residual (25ms frames, 5ms shift, order ≈ fs/1000 + 2)
#' 2. Bandpass filter signal around estimated F0 (mean-based signal)
#' 3. Find maxima/minima pairs in mean-based signal
#' 4. Locate GCI positions in LP residual peaks within windows
#'
#' **Interpretation**:
#' - GCI times: precise locations of glottal closure events (in seconds)
#' - Can be used as anchors for GCI-based analysis (e.g., trk_covarep_vq_gci)
#' - Prerequisite for voice quality measures that require GCI timing
#'
#' @examples
#' \dontrun{
#' # Single file
#' gcis <- lst_covarep_gci_sedreams("speech.wav", f0mean = 100)
#'
#' # Batch process
#' files <- c("file1.wav", "file2.wav")
#' results <- lst_covarep_gci_sedreams(files, f0mean = 110)
#'
#' # View results
#' results$gci_times[[1]]  # GCI times for first file (in seconds)
#' }
#'
#' @export
lst_covarep_gci_sedreams <- function(listOfFiles,
                                     beginTime = 0.0,
                                     endTime = 0.0,
                                     f0mean = NULL,
                                     polarity = NULL,
                                     verbose = TRUE) {

  nFiles <- length(listOfFiles)

  # Handle time parameters
  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime <- if (is.null(endTime)) 0.0 else endTime
  beginTime <- fast_recycle_times(beginTime, nFiles)
  endTime <- fast_recycle_times(endTime, nFiles)

  # Handle F0mean
  f0mean <- if (is.null(f0mean)) 100.0 else f0mean
  f0mean <- fast_recycle_times(f0mean, nFiles)

  # Handle polarity
  polarity <- if (is.null(polarity)) 1L else as.integer(polarity)
  polarity <- if (length(polarity) == 1) rep(polarity, nFiles) else polarity[1:nFiles]

  # Process each file
  results <- vector("list", nFiles)
  file_vec <- character(nFiles)
  n_gcis_vec <- integer(nFiles)
  gci_list <- vector("list", nFiles)

  for (idx in seq_len(nFiles)) {
    if (verbose) {
      cli::cli_progress_step("Processing {.file {basename(listOfFiles[idx])}} ({idx}/{nFiles})")
    }

    file_vec[idx] <- listOfFiles[idx]

    # Load audio
    audio_data <- av::read_audio_bin(
      listOfFiles[[idx]],
      channels = 1,
      start_time = if (beginTime[idx] > 0) beginTime[idx] else NULL,
      end_time = if (endTime[idx] > 0) endTime[idx] else NULL
    )

    fs <- attr(audio_data, "sample_rate")
    wave <- as.numeric(audio_data)

    # Auto-detect polarity if needed
    if (is.null(polarity) || is.na(polarity[idx])) {
      pol <- .gci_detect_polarity(wave, fs)
    } else {
      pol <- polarity[idx]
    }

    # Run GCI SEDREAMS algorithm
    gci_result <- .gci_sedreams(wave, fs, f0mean = f0mean[idx], polarity = pol)

    gci_times <- gci_result$gci
    n_gcis_vec[idx] <- length(gci_times)
    gci_list[[idx]] <- gci_times
  }

  if (verbose) {
    cli::cli_progress_done()
  }

  # Return as data frame
  data.frame(
    file = file_vec,
    n_gcis = n_gcis_vec,
    gci_times = I(gci_list)  # I() keeps list column
  )
}

# Internal: SEDREAMS GCI detection
.gci_sedreams <- function(wave, fs, f0mean = 100, polarity = 1) {
  wave <- polarity * as.numeric(wave)

  # Compute LPC residual
  frame_length <- round(25 / 1000 * fs)
  frame_shift <- round(5 / 1000 * fs)
  order <- round(fs / 1000) + 2L
  res <- get_lpc_residual_cpp(wave, frame_length, frame_shift, order)
  res[is.na(res)] <- 0

  # Mean-based signal extraction
  t0mean <- round(fs / f0mean)
  half_l <- round((1.7 * t0mean) / 2)
  blackwin <- as.numeric(av::blackman(2 * half_l + 1L))

  # Simple convolution for mean-based signal
  mean_based <- .gci_filter_signal_simple(wave, blackwin)

  if (length(mean_based) > 2 * half_l) {
    mean_based[(half_l + 1):(length(mean_based) - half_l)] <-
      mean_based[(1 + 2 * half_l):length(mean_based)]
  }
  mean_based[c(seq_len(half_l), seq.int(length(mean_based) - half_l + 1L, length(mean_based)))] <- 0

  # High-pass filter
  hp_filter <- .gci_fir_highpass(fs = fs, cutoff_hz = 50, order = max(32L, 2L * round(fs / 100)))
  mean_based <- .gci_filtfilt_fir(mean_based, hp_filter)

  # Normalize
  max_abs_mbs <- max(abs(mean_based))
  if (max_abs_mbs > 0) {
    mean_based <- mean_based / max_abs_mbs
  }

  # Find peaks
  pot_maxis <- .gci_find_peaks(mean_based)$index
  pot_minis <- .gci_find_peaks(-mean_based)$index

  if (!length(pot_maxis) || !length(pot_minis)) {
    return(list(gci = numeric()))
  }

  # Ensure alternation: max then min
  while (length(pot_maxis) && length(pot_minis) && pot_maxis[1] < pot_minis[1]) {
    pot_maxis <- pot_maxis[-1]
  }
  while (length(pot_minis) && length(pot_maxis) && pot_minis[length(pot_minis)] > pot_maxis[length(pot_maxis)]) {
    pot_minis <- pot_minis[-length(pot_minis)]
  }

  if (!length(pot_maxis) || !length(pot_minis)) {
    return(list(gci = numeric()))
  }

  minis <- pot_minis
  maxis <- pot_maxis

  # Normalize residual
  max_abs_res <- max(abs(res))
  if (max_abs_res > 0) {
    res <- res / max_abs_res
  }

  # Find median relative position of peaks
  posis <- which(res > 0.4)
  rel_posis <- numeric(length(posis))
  for (k in seq_along(posis)) {
    dists <- abs(minis - posis[k])
    pos <- which.min(dists)
    if (length(pos) && pos <= length(maxis)) {
      interv <- maxis[pos] - minis[pos]
      if (interv > 0) {
        rel_posis[k] <- (posis[k] - minis[pos]) / interv
      }
    }
  }
  ratio_gci <- if (length(rel_posis)) stats::median(rel_posis) else 0

  # Extract GCI positions
  gci <- numeric(length(minis))
  ind <- 1L
  for (k in seq_along(minis)) {
    interv <- maxis[k] - minis[k]
    start <- minis[k] + round((ratio_gci - 0.35) * interv)
    stop <- minis[k] + round((ratio_gci + 0.35) * interv)

    if (start < 1L) {
      start <- 1L
    } else if (start > length(res)) {
      break
    }
    if (stop > length(res)) {
      stop <- length(res)
    }
    if (stop > 1L) {
      vec <- res[start:stop]
      posi <- which.max(vec)
      gci[ind] <- start + posi[1] - 1L
      ind <- ind + 1L
    }
  }

  gci <- gci[gci > 0]
  list(gci = (gci - 1) / fs)
}

# Helper: Simple filter signal (convolution)
.gci_filter_signal_simple <- function(x, b) {
  # Use stats::filter for efficiency
  y <- stats::filter(x, b, method = "convolution", sides = 1)
  y[is.na(y)] <- 0
  y
}

# Helper: FIR high-pass filter design
.gci_fir_highpass <- function(fs, cutoff_hz = 490, order = NULL) {
  if (is.null(order)) {
    order <- max(32L, 2L * round(fs / 250))
  }
  if (order %% 2L == 1L) {
    order <- order + 1L
  }

  n <- 0:order
  mid <- order / 2
  fc <- cutoff_hz / fs
  h <- ifelse(
    n == mid,
    2 * fc,
    sin(2 * pi * fc * (n - mid)) / (pi * (n - mid))
  )
  w <- 0.54 - 0.46 * cos(2 * pi * n / order)
  lp <- h * w
  hp <- -lp
  hp[mid + 1L] <- hp[mid + 1L] + 1
  hp / sum(abs(hp))
}

# Helper: FIR filter application
.gci_fir_filter <- function(x, b) {
  y <- stats::filter(x, b, method = "convolution", sides = 1)
  y[is.na(y)] <- 0
  y
}

# Helper: Zero-phase (forward-backward) FIR filtering
.gci_filtfilt_fir <- function(x, b) {
  forward <- .gci_fir_filter(x, b)
  rev(.gci_fir_filter(rev(forward), b))
}

# Helper: Find peaks in signal
.gci_find_peaks <- function(x, min_distance = 1) {
  if (length(x) < 3) {
    return(list(index = integer(), value = numeric()))
  }

  # Find zero crossings in diff(sign(diff(x)))
  candidates <- which(diff(sign(diff(x))) < 0) + 1L
  if (!length(candidates)) {
    return(list(index = integer(), value = numeric()))
  }

  if (min_distance <= 1) {
    return(list(index = candidates, value = x[candidates]))
  }

  # Apply minimum distance constraint
  selected <- integer()
  for (idx in candidates[order(x[candidates], decreasing = TRUE)]) {
    if (!length(selected) || all(abs(idx - selected) >= min_distance)) {
      selected <- c(selected, idx)
    }
  }
  selected <- sort(selected)
  list(index = selected, value = x[selected])
}

# Helper: Detect signal polarity
.gci_detect_polarity <- function(wave, fs) {
  # Simple skewness-based polarity detection
  # Compute residual for both polarities and check which gives sparser result

  # Use residual skewness as proxy
  l <- round(25 / 1000 * fs)
  shift <- round(5 / 1000 * fs)
  order <- round(fs / 1000) + 2L

  # Try positive polarity
  res1 <- get_lpc_residual_cpp(wave, l, shift, order)
  skew1 <- .gci_skewness(res1)

  # Try negative polarity
  res2 <- get_lpc_residual_cpp(-wave, l, shift, order)
  skew2 <- .gci_skewness(res2)

  # Return polarity that maximizes skewness (more sparse)
  if (skew2 > skew1) -1L else 1L
}

# Helper: Skewness calculation
.gci_skewness <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(0)
  }
  centered <- x - mean(x)
  s <- stats::sd(centered)
  if (is.na(s) || s == 0) {
    return(0)
  }
  mean(centered^3) / (s^3)
}

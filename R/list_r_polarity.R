#' Signal polarity detection (RESKEW algorithm)
#'
#' Detects microphone/recording polarity from LP residual skewness.
#' The RESKEW algorithm compares residual skewness characteristics with
#' and without high-pass filtering to determine signal polarity.
#'
#' @param listOfFiles Vector of file paths (WAV, MP3, MP4, etc.) to analyze
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param toFile Write output to file (default: FALSE, not supported for scalar output)
#' @param verbose Show progress messages (default: TRUE)
#'
#' @return
#' Data frame with columns:
#' - `file`: Input file basename
#' - `polarity`: Polarity sign (+1 or -1)
#'
#' @details
#' **Algorithm** (RESKEW, Drugman et al.):
#' 1. High-pass filter signal at 490 Hz cutoff (removes low-frequency noise/drift)
#' 2. Compute LP residual with order = fs/1000 + 2 samples, 25ms frames, 5ms shift
#' 3. Compare residual skewness: unfiltered vs. high-pass filtered
#' 4. Polarity = sign(skew_filtered - skew_unfiltered)
#'
#' **Interpretation**:
#' - `+1`: Normal polarity (positive peaks are vocal pulses)
#' - `-1`: Inverted polarity (flip signal before further processing)
#'
#' **Note**: If skewness values are very close, result may be unstable.
#' Recommend confidence threshold: |result| > 0.1 for reliable polarity detection.
#'
#' @examples
#' \dontrun{
#' # Single file
#' pol <- lst_polarity("speech.wav")
#' cat("Polarity:", pol$polarity, "\n")
#'
#' # Batch process
#' files <- c("file1.wav", "file2.wav")
#' polarities <- lst_polarity(files)
#' }
#'
#' @export
lst_polarity <- function(listOfFiles,
                         beginTime = 0.0,
                         endTime = 0.0,
                         toFile = FALSE,
                         verbose = TRUE) {

  if (toFile) {
    cli::cli_abort("toFile is not supported for lst_polarity (scalar output)")
  }

  nFiles <- length(listOfFiles)

  # Handle time parameters
  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime <- if (is.null(endTime)) 0.0 else endTime
  beginTime <- fast_recycle_times(beginTime, nFiles)
  endTime <- fast_recycle_times(endTime, nFiles)

  # Process each file
  results <- data.frame(
    file = character(nFiles),
    polarity = integer(nFiles),
    stringsAsFactors = FALSE
  )

  for (idx in seq_len(nFiles)) {
    if (verbose) {
      cli::cli_progress_step("Processing {.file {basename(listOfFiles[idx])}} ({idx}/{nFiles})")
    }

    # Load audio
    audio_data <- av::read_audio_bin(
      listOfFiles[[idx]],
      channels = 1,
      start_time = if (beginTime[idx] > 0) beginTime[idx] else NULL,
      end_time = if (endTime[idx] > 0) endTime[idx] else NULL
    )

    fs <- attr(audio_data, "sample_rate")
    wave <- as.numeric(audio_data)

    # Compute polarity
    pol <- .polarity_reskew(wave, fs)

    results$file[idx] <- basename(listOfFiles[[idx]])
    results$polarity[idx] <- as.integer(pol)
  }

  if (verbose) {
    cli::cli_progress_done()
  }

  results
}

# Internal: RESKEW polarity detection algorithm
.polarity_reskew <- function(s, fs) {
  s <- as.numeric(s)

  # High-pass filter at 490 Hz
  hpf <- .polarity_fir_highpass(fs, cutoff_hz = 490)
  s_h <- .polarity_filtfilt(s, hpf)

  # LPC parameters
  frame_length <- round(25 / 1000 * fs)
  frame_shift <- round(5 / 1000 * fs)
  order <- round(fs / 1000) + 2L

  # LP residual on original signal (filtered with highpass coefficients applied to analysis)
  res1 <- .polarity_lpc_residual_two_signals(
    filter_signal = s,
    analysis_signal = s_h,
    frame_length = frame_length,
    frame_shift = frame_shift,
    order = order
  )

  # LP residual on original signal (analysis signal = original)
  res2 <- get_lpc_residual_cpp(
    wave = s,
    l = frame_length,
    shift = frame_shift,
    order = order
  )

  # Compare skewness
  skew1 <- .polarity_skewness(res1)
  skew2 <- .polarity_skewness(res2)

  # Polarity = sign of skewness difference
  sign(skew2 - skew1)
}

# FIR highpass filter design (Parks-McClellan-like Hamming window approach)
.polarity_fir_highpass <- function(fs, cutoff_hz = 490, order = NULL) {
  if (is.null(order)) {
    order <- max(32L, 2L * round(fs / 250))
  }
  if (order %% 2L == 1L) {
    order <- order + 1L
  }

  n <- 0:order
  mid <- order / 2
  fc <- cutoff_hz / fs

  # Ideal lowpass impulse response
  h <- ifelse(
    n == mid,
    2 * fc,
    sin(2 * pi * fc * (n - mid)) / (pi * (n - mid))
  )

  # Hamming window
  w <- 0.54 - 0.46 * cos(2 * pi * n / order)
  lp <- h * w

  # Transform to highpass by subtracting from impulse
  hp <- -lp
  hp[mid + 1L] <- hp[mid + 1L] + 1

  hp / sum(abs(hp))
}

# Forward-backward filtering (zero-phase FIR filtering)
.polarity_filtfilt <- function(x, b) {
  # Forward filter
  forward <- .polarity_fir_filter(x, b)
  # Reverse, filter again, reverse back (zero-phase)
  backward <- .polarity_fir_filter(rev(forward), b)
  rev(backward)
}

# Apply FIR filter via difference equation (causal)
# y[n] = sum(b[k] * x[n-k] for k = 0 to length(b)-1)
.polarity_fir_filter <- function(x, b) {
  x <- as.numeric(x)
  b <- as.numeric(b) / sum(abs(b))
  # Use stats::filter for causal FIR filtering
  # filter applies: y = filter(b, a, x) -> a[1]*y[n] = b[1]*x[n] + ... + b[m]*x[n-m+1]
  # Since a=1, this is: y[n] = b[1]*x[n] + b[2]*x[n-1] + ...
  y <- stats::filter(x, b, sides = 1)
  # Replace NAs with 0 (occurs at the beginning)
  y[is.na(y)] <- 0
  y
}

# LPC residual with two signals (filter_signal for LPC coefficients, analysis_signal for residual)
.polarity_lpc_residual_two_signals <- function(filter_signal, analysis_signal,
                                              frame_length, frame_shift, order) {
  filter_signal <- as.numeric(filter_signal)
  analysis_signal <- as.numeric(analysis_signal)

  n_frames <- floor((length(filter_signal) - frame_length) / frame_shift) + 1L
  residuals <- numeric()

  for (i in seq_len(n_frames)) {
    start_idx <- (i - 1L) * frame_shift + 1L
    end_idx <- start_idx + frame_length - 1L

    if (end_idx > length(filter_signal) || end_idx > length(analysis_signal)) {
      break
    }

    # Extract frame from filter_signal for LPC computation
    frame_filt <- filter_signal[start_idx:end_idx]

    # Apply Hann window
    w <- 0.5 * (1 - cos(2 * pi * (0:(frame_length - 1L)) / (frame_length - 1L)))
    frame_filt_windowed <- frame_filt * w

    # Compute LPC coefficients via autocorrelation method
    a <- .polarity_lpc(frame_filt_windowed, order)

    # Apply LPC filter to analysis_signal frame
    frame_ana <- analysis_signal[start_idx:end_idx]

    # Filter: a[1]*y[n] = b[0]*x[n] + b[1]*x[n-1] + ...
    # where b = c(1, -a[-1]) and a[1] = 1
    if (length(a) > 1) {
      filter_b <- c(1, -a[-1])
    } else {
      filter_b <- 1
    }

    res_frame <- stats::filter(filter_b, 1, frame_ana, method = "recursive")

    # Collect residuals (skip first few due to filter transient)
    residuals <- c(residuals, res_frame[!is.na(res_frame)])
  }

  residuals
}

# LPC via autocorrelation + Levinson-Durbin
.polarity_lpc <- function(x, order) {
  x <- as.numeric(x)
  n <- length(x)

  if (n <= order) {
    a <- rep(0, order + 1L)
    a[1L] <- 1
    return(a)
  }

  # Autocorrelation
  r <- numeric(order + 1L)
  for (k in 0:order) {
    idx <- 1:(n - k)
    if (length(idx) > 0) {
      r[k + 1L] <- sum(x[idx] * x[(k + 1):n])
    }
  }

  # Levinson-Durbin recursion
  a <- rep(0, order + 1L)
  a[1L] <- 1

  if (r[1L] <= 0) {
    return(a)
  }

  e <- r[1L]
  for (m in 1:order) {
    # Compute reflection coefficient k_m
    numerator <- r[m + 1L]
    for (j in 1:m) {
      numerator <- numerator - a[j] * r[m + 2L - j]
    }
    k_m <- numerator / e

    # Update coefficients
    a_old <- a
    for (j in 1:m) {
      a[j] <- a_old[j] - k_m * a_old[m + 2L - j]
    }
    a[m + 1L] <- -k_m

    # Update prediction error
    e <- e * (1 - k_m^2)
    if (e <= 1e-10) break
  }

  a
}

# Compute skewness (third moment / SD^3)
.polarity_skewness <- function(x) {
  x <- as.numeric(x[!is.na(x)])
  if (length(x) < 3) {
    return(0)
  }
  centered <- x - mean(x)
  s <- stats::sd(centered)
  if (is.na(s) || s == 0) {
    return(0)
  }
  mean(centered^3) / (s^3)
}

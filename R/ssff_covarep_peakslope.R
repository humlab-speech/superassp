#' D'Alessandro spectral tilt correlate (PeakSlope)
#'
#' Computes frame-level spectral tilt using D'Alessandro Morlet wavelet
#' decomposition. PeakSlope is a voice quality correlate that tracks
#' spectral energy distribution across frequency scales.
#'
#' @param listOfFiles Vector of file paths (WAV, MP3, MP4, etc.) to analyze
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param toFile Write output to file (TRUE) or return object (FALSE). Default: FALSE
#' @param explicitExt Output file extension (default: "psl")
#' @param outputDirectory Output directory (NULL for same as input file)
#' @param verbose Show progress messages (default: TRUE)
#'
#' @return
#' If `toFile=FALSE` (default): AsspDataObj with columns `peakslope` (spectral tilt slope)
#'
#' If `toFile=TRUE`: invisibly returns vector of output file paths
#'
#' @details
#' **Algorithm** \insertCite{Henrich2001}{superassp}:
#' 1. Apply Morlet wavelets at 7 scales (2^0 to 2^6)
#' 2. For each frame (40ms, 10ms shift): extract maximum magnitude across scales
#' 3. Convert magnitudes to log10 scale
#' 4. Fit linear regression across scales: log(mag) ~ scale
#' 5. Slope coefficient = peakslope value (spectral tilt)
#'
#' **Interpretation**:
#' - Negative slope: energy concentrated at low scales (bright voice, high spectral tilt)
#' - Positive slope: energy concentrated at high scales (breathy voice, low spectral tilt)
#' - Magnitude indicates overall spectral energy distribution steepness
#'
#' **Related measures**: CPP (cepstral peak prominence), HNR (harmonic-to-noise ratio)
#' PeakSlope captures similar voice quality information but via wavelet decomposition
#' instead of cepstral or harmonic analysis.
#'
#' @examples
#' \dontrun{
#' # Single file, return object
#' peakslope <- trk_peakslope("speech.wav", toFile = FALSE)
#'
#' # Batch process multiple files
#' files <- c("file1.wav", "file2.wav")
#' trk_peakslope(files, toFile = TRUE, outputDirectory = "output/")
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @export
trk_peakslope <- function(listOfFiles,
                          beginTime = 0.0,
                          endTime = 0.0,
                          toFile = FALSE,
                          explicitExt = "psl",
                          outputDirectory = NULL,
                          verbose = TRUE) {

  nFiles <- length(listOfFiles)

  # Handle time parameters
  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime <- if (is.null(endTime)) 0.0 else endTime
  beginTime <- fast_recycle_times(beginTime, nFiles)
  endTime <- fast_recycle_times(endTime, nFiles)

  # Process each file
  externalRes <- vector("list", nFiles)

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

    # Run peakslope algorithm
    ps_res <- .peakslope_analysis(wave, fs)

    # Convert to AsspDataObj
    assp_obj <- list(
      peakslope = ps_res[, 2]
    )

    class(assp_obj) <- "AsspDataObj"
    attr(assp_obj, "sampleRate") <- fs
    attr(assp_obj, "startTime") <- as.numeric(beginTime[idx])
    attr(assp_obj, "startRecord") <- 1L
    attr(assp_obj, "endRecord") <- length(ps_res[, 2])
    attr(assp_obj, "trackFormats") <- "FLOAT"

    if (toFile) {
      # Write to SSFF file
      out_dir <- if (is.null(outputDirectory)) dirname(listOfFiles[[idx]]) else outputDirectory
      out_file <- file.path(
        out_dir,
        paste0(tools::file_path_sans_ext(basename(listOfFiles[[idx]])), ".", explicitExt)
      )
      write.AsspDataObj(assp_obj, out_file)
      externalRes[[idx]] <- out_file
    } else {
      externalRes[[idx]] <- assp_obj
    }
  }

  if (verbose) {
    cli::cli_progress_done()
  }

  # Return results
  if (toFile) {
    invisible(unlist(externalRes))
  } else {
    if (nFiles == 1) externalRes[[1]] else externalRes
  }
}

# Internal: PeakSlope analysis
.peakslope_analysis <- function(s, fs) {
  frame_len_ms <- 40
  frame_shift_ms <- 10
  frame_len <- round(frame_len_ms / 1000 * fs)
  frame_shift <- round(frame_shift_ms / 1000 * fs)

  times <- numeric()
  slopes <- numeric()
  i <- 0:6
  y <- matrix(0, nrow = length(i), ncol = length(s))

  # Compute D'Alessandro wavelets at 7 scales
  for (n in seq_along(i)) {
    h_i <- .peakslope_daless_mw(i, n, fs)
    y[n, ] <- .peakslope_do_daless_filt(s, h_i)
  }

  # Frame-wise spectral slope extraction
  start <- 1L
  finish <- start + frame_len - 1L

  while (finish <= length(s)) {
    # Maximum magnitude across scales within frame
    maxima <- apply(abs(y[, start:finish, drop = FALSE]), 1, max)

    # Convert to log10 scale (add epsilon for numerical stability)
    maxima <- log10(pmax(rev(maxima), .Machine$double.eps))

    # Fit linear regression: log_mag ~ scale
    # seq_along gives scale indices (1-7)
    p <- stats::lm(maxima ~ seq_along(maxima))
    slopes <- c(slopes, unname(stats::coef(p)[2]))

    # Frame center time
    times <- c(times, (round((start + finish) / 2) - 1) / fs)

    start <- start + frame_shift
    finish <- start + frame_len - 1L
  }

  # Return [times, slopes] matrix
  ps <- cbind(times, slopes)
  ps[is.nan(ps[, 2]), 2] <- 0
  unname(ps)
}

# D'Alessandro Morlet wavelet
.peakslope_daless_mw <- function(i, n, fs) {
  s <- 2^i
  f_o <- fs / 2
  tau <- 1 / (2 * f_o)
  t <- (-1000:1000) / fs
  -cos(2 * pi * f_o * (t / s[n])) * exp(-((t / s[n])^2) / (2 * (tau^2)))
}

# Filter signal with Morlet wavelet (manual convolution for clarity)
.peakslope_do_daless_filt <- function(s, h_i) {
  # Pad signal with zeros
  s_padded <- c(s, rep(0, length(h_i) - 1L))

  # Causal FIR filter: y[n] = sum(h[k] * x[n-k] for k=0 to len(h)-1)
  # Using stats::filter with only b and a parameters (no deprecated sides)
  y <- numeric(length(s_padded))
  for (i in seq_along(s_padded)) {
    start_idx <- max(1L, i - length(h_i) + 1L)
    h_idx_start <- length(h_i) - (i - start_idx)
    if (h_idx_start > 0) {
      y[i] <- sum(h_i[h_idx_start:length(h_i)] * s_padded[start_idx:i])
    }
  }

  # Extract valid part (skip filter transient)
  half_len <- ceiling(length(h_i) / 2)
  y[half_len:(half_len + length(s) - 1L)]
}

# Set function attributes
attr(trk_peakslope, "ext") <- "psl"
attr(trk_peakslope, "tracks") <- c("peakslope")
attr(trk_peakslope, "outputType") <- "SSFF"
attr(trk_peakslope, "nativeFiletypes") <- character()

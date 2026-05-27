#' Track spectral tilt using D'Alessandro PeakSlope (Morlet wavelet)
#'
#' Returns a per-frame spectral tilt estimate related to breathiness and voice
#' quality \insertCite{Henrich2001}{superassp}. Captures similar information to
#' CPP and H1-H2 but via multi-scale wavelet analysis.
#'
#' @details
#' Spectral tilt is estimated by fitting a linear slope across peak magnitudes
#' at seven Morlet wavelet scales. Lower (more negative) values indicate a steeper
#' spectral slope, associated with creakier phonation; higher values indicate
#' breathier or modal phonation.
#'
#' @inheritParams trk_acf
#' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
#'   paths written invisibly. If \code{FALSE}, return an \code{AsspDataObj}.
#'   Default \code{FALSE}.
#' @param explicitExt Character. Output file extension. Default \code{"psl"}.
#'
#' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
#'   \describe{
#'     \item{\code{peakslope}}{FLOAT, spectral tilt slope (log10(magnitude) per scale
#'       unit), n_frames × 1. Negative = energy at low scales (bright/tense voice);
#'       positive = energy at high scales (breathy/relaxed voice).}
#'   }
#'   Frame rate: 100 Hz (10 ms hop; 40 ms analysis window).
#'   If \code{toFile = TRUE}: character vector of output file paths, returned invisibly.
#'
#' @details
#' Seven D'Alessandro Morlet wavelets at octave-spaced scales (2^0 to 2^6) are
#' convolved with the signal. The maximum magnitude in each 40 ms frame is taken per
#' scale, log10-transformed, and a linear regression slope is computed across the
#' seven scale-magnitude pairs.
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
#' @references \insertAllCited{}
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

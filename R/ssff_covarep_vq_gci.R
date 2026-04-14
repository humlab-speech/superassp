#' GCI-Based Voice Quality Parameters
#'
#' Computes voice quality parameters per glottal closure instant (GCI),
#' interpolated to regular 10ms frame grid.
#' Complements `lst_covarep_vq()` (scalar summaries) with full time series.
#'
#' @param listOfFiles Vector of file paths (WAV, MP3, MP4, etc.) to analyze
#' @param gci_times Optional list of GCI times (seconds). If NULL, computed via SEDREAMS.
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param toFile Write output to file (TRUE) or return object (FALSE). Default: FALSE
#' @param explicitExt Output file extension (default: "vqg")
#' @param outputDirectory Output directory (NULL for same as input file)
#' @param verbose Show progress messages (default: TRUE)
#'
#' @return
#' If `toFile=FALSE` (default): AsspDataObj with 5 tracks (naq, qoq, h1h2, hrf, psp),
#' each interpolated to 10ms frame grid.
#'
#' If `toFile=TRUE`: invisibly returns vector of output file paths
#'
#' @details
#' **Voice Quality Parameters** (computed per GCI, interpolated to 10ms grid):
#' - **NAQ** (Normalized Amplitude Quotient): Amplitude shape, related to open quotient
#' - **QOQ** (Quasi-Open Quotient): Duration of open phase relative to pitch period
#' - **H1H2** (Fundamental to 2nd Harmonic): Spectral tilt measure
#' - **HRF** (Harmonic Richness Factor): Overall harmonic content
#' - **PSP** (Parabolic Spectral Peak): Spectral peakedness
#'
#' **Method**:
#' 1. Detect GCIs (SEDREAMS) if not provided
#' 2. Compute glottal flow via IAIF (Inverse Filtering)
#' 3. Extract voice quality measures per GCI
#' 4. Interpolate to regular 10ms frame grid using linear interpolation
#'
#' **Interpretation**:
#' - NAQ/QOQ: Phonation mode indicators (creaky, normal, breathy)
#' - H1H2/HRF: Spectral characteristics (bright, rich, sparse voice)
#' - PSP: Voice quality correlate (higher = smoother spectral envelope)
#'
#' @examples
#' \dontrun{
#' # Single file, auto GCI detection
#' vq <- trk_covarep_vq_gci("speech.wav", toFile = FALSE)
#'
#' # With pre-computed GCIs
#' gcis <- lst_covarep_gci_sedreams("speech.wav", f0mean = 100)
#' gci_times <- gcis$gci_times[[1]]
#' vq <- trk_covarep_vq_gci("speech.wav", gci_times = gci_times, toFile = FALSE)
#' }
#'
#' @export
trk_covarep_vq_gci <- function(listOfFiles,
                               gci_times = NULL,
                               beginTime = 0.0,
                               endTime = 0.0,
                               toFile = FALSE,
                               explicitExt = "vqg",
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

    # Get GCIs if not provided
    if (is.null(gci_times)) {
      gci_result <- lst_covarep_gci_sedreams(
        listOfFiles[[idx]],
        beginTime = if (beginTime[idx] > 0) beginTime[idx] else 0,
        endTime = if (endTime[idx] > 0) endTime[idx] else 0,
        f0mean = 100,
        verbose = FALSE
      )
      gci_sec <- gci_result$gci_times[[1]]
    } else if (is.list(gci_times)) {
      gci_sec <- gci_times[[idx]]
    } else {
      gci_sec <- gci_times
    }

    # Compute IAIF (glottal flow + derivative)
    iaif_result <- trk_covarep_iaif(
      listOfFiles[[idx]],
      beginTime = if (beginTime[idx] > 0) beginTime[idx] else 0,
      endTime = if (endTime[idx] > 0) endTime[idx] else 0,
      toFile = FALSE,
      verbose = FALSE
    )

    if (!is.null(iaif_result) && length(iaif_result) > 0) {
      gf <- as.numeric(iaif_result$g)  # glottal flow
      gfd <- as.numeric(iaif_result$dg)  # glottal flow derivative
    } else {
      # If IAIF fails, skip this file
      if (verbose) {
        cli::cli_warn("IAIF computation failed for {.file {basename(listOfFiles[idx])}}")
      }
      externalRes[[idx]] <- NULL
      next
    }

    # Compute VQ parameters at GCIs
    if (length(gci_sec) > 0) {
      vq_result <- .vqg_get_vq_params(gf, gfd, fs, gci_sec)

      # Interpolate to 10ms frame grid
      frame_shift <- round(10 / 1000 * fs)
      n_frames <- ceiling(length(wave) / frame_shift)

      frame_times <- (0:(n_frames - 1)) * 10 / 1000

      # Interpolate each VQ parameter
      naq_interp <- .vqg_interpolate_gci_values(
        gci_sec, vq_result$naq, frame_times, fs
      )
      qoq_interp <- .vqg_interpolate_gci_values(
        gci_sec, vq_result$qoq, frame_times, fs
      )
      h1h2_interp <- .vqg_interpolate_gci_values(
        gci_sec, vq_result$h1h2, frame_times, fs
      )
      hrf_interp <- .vqg_interpolate_gci_values(
        gci_sec, vq_result$hrf, frame_times, fs
      )
      psp_interp <- .vqg_interpolate_gci_values(
        gci_sec, vq_result$psp, frame_times, fs
      )

      # Convert to AsspDataObj
      assp_obj <- list(
        naq = naq_interp,
        qoq = qoq_interp,
        h1h2 = h1h2_interp,
        hrf = hrf_interp,
        psp = psp_interp
      )

      class(assp_obj) <- "AsspDataObj"
      attr(assp_obj, "sampleRate") <- fs
      attr(assp_obj, "startTime") <- as.numeric(beginTime[idx])
      attr(assp_obj, "startRecord") <- 1L
      attr(assp_obj, "endRecord") <- n_frames
      attr(assp_obj, "trackFormats") <- rep("FLOAT", 5)
    } else {
      # No GCIs detected, return empty AsspDataObj
      assp_obj <- list(
        naq = numeric(),
        qoq = numeric(),
        h1h2 = numeric(),
        hrf = numeric(),
        psp = numeric()
      )
      class(assp_obj) <- "AsspDataObj"
      attr(assp_obj, "sampleRate") <- fs
      attr(assp_obj, "startTime") <- as.numeric(beginTime[idx])
      attr(assp_obj, "startRecord") <- 1L
      attr(assp_obj, "endRecord") <- 0L
      attr(assp_obj, "trackFormats") <- rep("FLOAT", 5)
    }

    if (toFile) {
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
    invisible(unlist(externalRes[!sapply(externalRes, is.null)]))
  } else {
    if (nFiles == 1) externalRes[[1]] else externalRes
  }
}

# Helper: Extract VQ parameters from glottal flow at GCI locations
.vqg_get_vq_params <- function(gf, gfd, fs, gci) {
  gci_samp <- round(as.numeric(gci) * fs) + 1L

  f0min <- 20
  f0max <- 500
  naq <- numeric(length(gci_samp))
  qoq <- numeric(length(gci_samp))
  h1h2 <- numeric(length(gci_samp))
  hrf <- numeric(length(gci_samp))
  psp <- numeric(length(gci_samp))

  glot_shift <- round(0.5 / 1000 * fs)
  qoq_level <- 0.5
  t0_num <- 3
  min_harm_num <- 5
  hrf_freq_max <- 5000
  psp_fft_size <- 2048

  # Helper function
  find_amid_t <- function(glot_adj, amid, tz) {
    t1 <- 0
    t2 <- 0
    if (tz != 0 && tz >= 1 && tz <= length(glot_adj)) {
      n <- tz
      while (n > 3 && glot_adj[n] > amid) {
        n <- n - 1L
      }
      t1 <- n

      n <- tz
      while (n < length(glot_adj) - 2 && glot_adj[n] > amid) {
        n <- n + 1L
      }
      t2 <- n
    }
    c(t1, t2)
  }

  # PSP helper
  psp_get_a <- function(x) {
    ne_min <- 0.01
    flag <- FALSE
    n <- 3L
    x <- as.numeric(x)
    x2 <- x^2
    k2 <- (0:100)^2

    k2_sum <- sum(k2[1:(n - 1L)])
    k4_sum <- sum(k2[1:(n - 1L)]^2)
    x_k1_sum <- sum(x[1:(n - 1L)])
    x_k1_k2_sum <- sum(x[1:(n - 1L)] * k2[1:(n - 1L)])
    x2_k1_sum <- sum(x2[1:(n - 1L)])

    while (!flag) {
      if (n > length(k2)) {
        k2 <- (0:(length(k2) * 2))^2
      }

      k2_sum <- k2_sum + k2[n]
      k4_sum <- k4_sum + k2[n] * k2[n]
      x_k1_sum <- x_k1_sum + x[n]
      x_k1_k2_sum <- x_k1_k2_sum + x[n] * k2[n]
      x2_k1_sum <- x2_k1_sum + x2[n]

      a <- (n * x_k1_k2_sum - x_k1_sum * k2_sum) / (n * k4_sum - k2_sum^2)
      x_k1_a_k2 <- x[1:n] - a * k2[1:n]
      b <- mean(x_k1_a_k2)
      ne <- sum((x_k1_a_k2 - b)^2) / x2_k1_sum

      if (ne < ne_min) {
        n <- n + 1L
      } else {
        flag <- TRUE
      }
    }
    a
  }

  # Magnitude to dB conversion
  mag2db <- function(x) {
    20 * log10(pmax(as.numeric(x), .Machine$double.eps))
  }

  # Peak finding helper
  find_peaks <- function(x, min_distance = 1) {
    if (length(x) < 3) {
      return(list(index = integer(), value = numeric()))
    }
    candidates <- which(diff(sign(diff(x))) < 0) + 1L
    if (!length(candidates)) {
      return(list(index = integer(), value = numeric()))
    }
    if (min_distance <= 1) {
      return(list(index = candidates, value = x[candidates]))
    }
    selected <- integer()
    for (idx in candidates[order(x[candidates], decreasing = TRUE)]) {
      if (!length(selected) || all(abs(idx - selected) >= min_distance)) {
        selected <- c(selected, idx)
      }
    }
    selected <- sort(selected)
    list(index = selected, value = x[selected])
  }

  # Process each GCI
  for (n in seq_along(gci_samp)) {
    if (n == 1L) {
      start <- 1L
      stop <- gci_samp[n]
      if (length(gci_samp) < 2L) {
        next
      }
      t0 <- gci_samp[n + 1L] - gci_samp[n]
    } else {
      start <- gci_samp[n - 1L]
      stop <- gci_samp[n]
      t0 <- gci_samp[n] - gci_samp[n - 1L]
    }

    f0 <- fs / t0

    if (is.finite(f0) && t0 != 0 && f0 > f0min && f0 < f0max && start >= 1 && stop <= length(gf)) {
      # NAQ and QOQ
      gf_vals <- c(gf[start], gf[stop])
      if (start != stop && length(gf_vals) > 1L) {
        line <- stats::approx(
          x = c(1, stop - start + 1L),
          y = gf_vals,
          xout = 1:(stop - start + 1L),
          method = "linear",
          rule = 2
        )$y
      } else {
        line <- gf_vals[1]
      }

      gf_seg <- gf[start:stop]
      gf_seg_comp <- gf_seg - line

      d_peak <- max(abs(gfd[start:min(stop + glot_shift, length(gfd))]))
      if (d_peak == 0) d_peak <- 1
      max_idx <- which.max(gf_seg_comp)
      f_ac <- gf_seg_comp[max_idx]
      amid <- f_ac * qoq_level
      t_amid <- find_amid_t(gf_seg_comp, amid, max_idx)

      naq[n] <- (f_ac / d_peak) / t0
      qoq[n] <- (t_amid[2] - t_amid[1]) / (fs / f0)

      # H1H2 and HRF
      f_start <- max(gci_samp[n] - round((t0 * t0_num) / 2), 1L)
      f_stop <- min(gci_samp[n] + round((t0 * t0_num) / 2), length(gfd))
      if (f_stop > f_start) {
        f_frame <- gfd[f_start:f_stop] * (2^15)
        ham_win <- as.numeric(av::hamming(length(f_frame)))
        f_win <- f_frame * ham_win
        length(f_win) <- fs
        f_win[is.na(f_win)] <- 0
        f_spec <- mag2db(Mod(stats::fft(f_win)))

        peaks <- find_peaks(f_spec, min_distance = max(1, round(f0 / 2)))
        hrf_harm_num <- floor(hrf_freq_max / f0)
        if (length(peaks$index) >= min_harm_num && hrf_harm_num >= 2) {
          f0_idx <- vapply(
            (1:hrf_harm_num) * f0,
            function(target) which.min(abs(peaks$index - target)),
            integer(1)
          )
          if (length(f0_idx) >= 2) {
            h1h2[n] <- peaks$value[f0_idx[1]] - peaks$value[f0_idx[2]]
            hrf[n] <- sum(peaks$value[f0_idx[-1]]) / peaks$value[f0_idx[1]]
          }
        }
      }

      # PSP
      if (n > 1L) {
        start_psp <- gci_samp[n - 1L]
        stop_psp <- gci_samp[n]
        if (start_psp >= 1 && stop_psp <= length(gf)) {
          gf_seg <- gf[start_psp:stop_psp]
          gf_seg <- gf_seg - min(gf_seg)
          max_gf <- max(gf_seg)
          if (max_gf > 0) {
            gf_seg <- gf_seg / max_gf
          }
          length(gf_seg) <- psp_fft_size
          gf_seg[is.na(gf_seg)] <- 0
          x_spec <- mag2db(Mod(stats::fft(gf_seg)))
          x_spec <- x_spec[seq_along(x_spec) <= floor(fs / 2)]
          a <- psp_get_a(x_spec)

          a_max_sig <- rep(1 / round(fs / f0), round(fs / f0))
          length(a_max_sig) <- psp_fft_size
          a_max_sig[is.na(a_max_sig)] <- 0
          x_a_max <- mag2db(Mod(stats::fft(a_max_sig)))
          x_a_max <- x_a_max[seq_along(x_a_max) <= floor(fs / 2)]
          a_max <- psp_get_a(x_a_max)

          if (a_max != 0) {
            psp[n] <- a / a_max
          }
        }
      }
    }
  }

  list(
    naq = naq,
    qoq = qoq,
    h1h2 = h1h2,
    hrf = hrf,
    psp = psp
  )
}

# Helper: Interpolate GCI-based values to regular frame grid
.vqg_interpolate_gci_values <- function(gci_times, gci_values, frame_times, fs) {
  # Remove NaN/Inf values
  valid_idx <- is.finite(gci_values) & gci_values != 0
  if (!any(valid_idx)) {
    # No valid values, return zeros
    return(rep(0, length(frame_times)))
  }

  gci_valid <- gci_times[valid_idx]
  val_valid <- gci_values[valid_idx]

  # Linear interpolation
  interp_result <- stats::approx(
    x = gci_valid,
    y = val_valid,
    xout = frame_times,
    method = "linear",
    rule = 2
  )

  interp_result$y[is.na(interp_result$y)] <- 0
  as.numeric(interp_result$y)
}

# Set function attributes
attr(trk_covarep_vq_gci, "ext") <- "vqg"
attr(trk_covarep_vq_gci, "tracks") <- c("naq", "qoq", "h1h2", "hrf", "psp")
attr(trk_covarep_vq_gci, "outputType") <- "SSFF"
attr(trk_covarep_vq_gci, "nativeFiletypes") <- character()

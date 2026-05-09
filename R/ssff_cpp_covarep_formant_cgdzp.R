#' Chirp Group Delay Zero-Phase formant tracker
#'
#' Computes formant frequencies using the Chirp Group Delay Zero-Phase (CGDZP)
#' method from COVAREP. This algorithm applies zero-phase processing and exponential
#' envelope warping to estimate formants from the chirp group delay.
#'
#' @param listOfFiles Vector of file paths (WAV, MP3, MP4, etc.) to analyze
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param frameSize Frame length in milliseconds (default: 30)
#' @param frameShift Frame shift in milliseconds (default: 10)
#' @param toFile Write output to file (TRUE) or return object (FALSE). Default: FALSE
#' @param explicitExt Output file extension (default: "cgf")
#' @param outputDirectory Output directory (NULL for same as input file)
#' @param verbose Show progress messages (default: TRUE)
#'
#' @return
#' If `toFile=FALSE` (default): AsspDataObj with columns F1–F5 (formant frequencies in Hz)
#'
#' If `toFile=TRUE`: invisibly returns vector of output file paths
#'
#' @details
#' The CGDZP formant tracker uses:
#' \itemize{
#'   \item Blackman window (frame analysis)
#'   \item Zero-phase processing (magnitude FFT inverse)
#'   \item Exponential envelope warping (smooth pitch-dependent emphasis)
#'   \item Chirp group delay detection (formant frequency localization)
#'   \item Temporal continuity constraint (250 Hz max frame-to-frame jump)
#' }
#'
#' **Algorithm**: Zero-phase processing eliminates phase distortion, making the
#' chirp group delay a reliable formant frequency estimator. The exponential
#' envelope warping (controlled by warping factor r) adjusts for pitch-dependent
#' spectral variation. Blackman window (from av package) provides spectral
#' sidelobe control. Works well on voiced speech, particularly in nasals and
#' fricative transitions where LPC-based methods (Burg) struggle.
#'
#' **Outputs**: 5 formant frequency tracks (F1–F5) in Hz, or zero if unvoiced/unreliable.
#'
#' @examples
#' \dontrun{
#' # Single file, return object
#' formants <- trk_formant_cgdzp("speech.wav", toFile = FALSE)
#'
#' # Batch process multiple files
#' files <- c("file1.wav", "file2.wav")
#' trk_formant_cgdzp(files, toFile = TRUE, outputDirectory = "output/")
#'
#' # Custom frame parameters
#' formants <- trk_formant_cgdzp("speech.wav",
#'                              frameSize = 40,
#'                              frameShift = 5,
#'                              toFile = FALSE)
#' }
#'
#' @export
trk_formant_cgdzp <- function(listOfFiles,
                             beginTime = 0.0,
                             endTime = 0.0,
                             frameSize = 30,
                             frameShift = 10,
                             toFile = FALSE,
                             explicitExt = "cgf",
                             outputDirectory = NULL,
                             verbose = TRUE) {

  # Validate inputs
  if (!is.numeric(frameSize) || frameSize <= 0) {
    cli::cli_abort("frameSize must be positive")
  }
  if (!is.numeric(frameShift) || frameShift <= 0) {
    cli::cli_abort("frameShift must be positive")
  }

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

    # Run CGDZP algorithm
    cgdzp_res <- .cgdzp_formant_analysis(wave, fs, frameSize, frameShift)

    # Convert to AsspDataObj
    assp_obj <- list(
      F1 = cgdzp_res$formant_peaks[, 1],
      F2 = cgdzp_res$formant_peaks[, 2],
      F3 = cgdzp_res$formant_peaks[, 3],
      F4 = cgdzp_res$formant_peaks[, 4],
      F5 = cgdzp_res$formant_peaks[, 5]
    )

    class(assp_obj) <- "AsspDataObj"
    attr(assp_obj, "sampleRate") <- fs
    attr(assp_obj, "startTime") <- as.numeric(beginTime[idx])
    attr(assp_obj, "startRecord") <- 1L
    attr(assp_obj, "endRecord") <- length(cgdzp_res$formant_peaks[, 1])
    attr(assp_obj, "trackFormats") <- rep("SHORT", 5)

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

# Internal: CGDZP formant analysis
.cgdzp_formant_analysis <- function(wave, fs, frame_size_ms = 30, frame_shift_ms = 10) {
  num_formants <- 5L
  frame_size <- round(fs / 1000 * frame_size_ms)
  frame_shift <- round(fs / 1000 * frame_shift_ms)
  fs_lr <- 2048L
  view_range <- round(fs_lr / 3.2)
  r_fix <- 1.12
  max_formant_delta <- 250

  n <- 0:(frame_size - 2L)
  size_wave <- length(wave)
  num_frames <- floor((size_wave - frame_size) / frame_shift)
  t_analysis <- numeric(num_frames)

  # Blackman window from av package
  black_win <- av::blackman(frame_size)
  formant_peaks <- matrix(0, nrow = num_frames, ncol = num_formants)

  for (kk in 0:(num_frames - 1L)) {
    speech_data <- wave[(kk * frame_shift + 1L):(kk * frame_shift + frame_size)]
    windowed_data <- speech_data * black_win
    zero_phase_data <- Re(.cgdzp_fft_vector(Mod(.cgdzp_fft_vector(diff(windowed_data))), inverse = TRUE))

    num_peaks <- 0L
    r_cur <- r_fix
    peak_index <- integer()
    while (num_peaks != num_formants && r_cur > 1.01 && r_cur < 1.25) {
      exponential_envelope <- exp(log(1 / r_cur) * n)
      fourier_trans <- .cgdzp_fft_vector(zero_phase_data * exponential_envelope, n = fs_lr)
      ang_fft <- Arg(fourier_trans[1:view_range])
      chirp_group_delay <- -diff(ang_fft)

      peak_index <- .cgdzp_formant_peak_pick(chirp_group_delay, 1)
      num_peaks <- length(peak_index)
      if (num_peaks > num_formants && r_cur >= r_fix) {
        r_cur <- r_cur + 0.01
        peak_index <- peak_index[1:num_formants]
      } else if (num_peaks < num_formants && r_cur <= r_fix) {
        r_cur <- r_cur - 0.01
        peak_index <- c(peak_index, rep(0, num_formants - num_peaks))
      } else {
        break
      }
    }

    if (num_peaks > num_formants) {
      peak_index <- peak_index[1:num_formants]
    } else if (num_peaks < num_formants) {
      peak_index <- c(peak_index, rep(0, num_formants - length(peak_index)))
    }

    formant_peaks[kk + 1L, ] <- peak_index
    t_analysis[kk + 1L] <- round((kk * frame_shift + 1L + kk * frame_shift + frame_size) / 2) / fs
  }

  formant_peaks <- round(formant_peaks * fs / fs_lr)
  formant_peaks_cost <- formant_peaks * 0

  if (num_frames >= 5L) {
    for (kk in 3:(num_frames - 2L)) {
      pre_pre <- formant_peaks[kk - 2L, ]
      post_post <- formant_peaks[kk + 2L, ]
      pre <- formant_peaks[kk - 1L, ]
      post <- formant_peaks[kk + 1L, ]
      current <- formant_peaks[kk, ]
      current_cost <- current * 0

      for (mm in seq_len(num_formants)) {
        if (current[mm] == 0) {
          current_cost[mm] <- fs / 2
        } else {
          distance_pre <- sort(abs(pre - current[mm]))
          distance_post <- sort(abs(post - current[mm]))
          distance_pre_pre <- sort(abs(pre_pre - current[mm]))
          distance_post_post <- sort(abs(post_post - current[mm]))
          all_distances <- c((distance_pre[1] + distance_post[1]) / 2, distance_pre[1], distance_post[1])
          all_distances2 <- c((distance_pre_pre[1] + distance_post_post[1]) / 2, distance_pre_pre[1], distance_post_post[1])
          current_cost[mm] <- min(c(all_distances, all_distances2))
        }
      }
      formant_peaks_cost[kk, ] <- current_cost
    }

    for (kk in seq_len(num_frames)) {
      for (mm in seq_len(num_formants)) {
        if (formant_peaks_cost[kk, mm] > max_formant_delta) {
          formant_peaks[kk, mm] <- 0
        }
      }
    }

    for (kk in 2:(num_frames - 1L)) {
      current <- formant_peaks[kk, ]
      index <- which(current != 0)
      non_zero <- current[index]
      num_non_zero <- length(index)
      num_zero <- num_formants - num_non_zero

      if (num_non_zero < num_formants && length(index)) {
        possible_values <- sort(c(formant_peaks[kk - 1L, ], formant_peaks[kk + 1L, ]))
        while (length(possible_values) > 0 && possible_values[1] == 0) {
          possible_values <- possible_values[-1]
        }

        possible_candidates <- numeric()
        for (mm in seq_along(possible_values)) {
          distance_array <- sort(abs(non_zero - possible_values[mm]))
          if (!length(distance_array) || distance_array[1] > max_formant_delta) {
            possible_candidates <- c(possible_candidates, possible_values[mm])
          }
        }

        len_candidates <- length(possible_candidates)
        if (len_candidates <= num_zero) {
          current <- sort(c(non_zero, possible_candidates, rep(0, num_zero - len_candidates)))
        } else if (num_zero == 1L) {
          index2 <- order(diff(possible_candidates))
          current <- sort(c(non_zero, possible_candidates[index2[1]]))
        } else if (num_zero < len_candidates) {
          current <- sort(c(non_zero, possible_candidates[1:num_zero]))
        }
        formant_peaks[kk, ] <- current
      }
    }
  }

  list(formant_peaks = formant_peaks, t_analysis = t_analysis)
}

# Helpers
.cgdzp_fft_vector <- function(x, n = length(x), inverse = FALSE) {
  x <- as.complex(x)
  length(x) <- n
  x[is.na(x)] <- 0
  y <- stats::fft(x, inverse = inverse)
  if (inverse) {
    y <- y / n
  }
  y
}

.cgdzp_formant_peak_pick <- function(diff_phase, min_peak_dist = 1) {
  peak_index <- integer()
  lendiff_phase <- length(diff_phase)

  if (lendiff_phase < 13) {
    return(peak_index)
  }

  for (kk in 6:(lendiff_phase - 6)) {
    if ((diff_phase[kk] >= diff_phase[kk - 1] && diff_phase[kk] >= diff_phase[kk + 1]) &&
        (diff_phase[kk] >= diff_phase[kk - 2] && diff_phase[kk] >= diff_phase[kk + 2]) &&
        (diff_phase[kk] > diff_phase[kk - 3] && diff_phase[kk] > diff_phase[kk + 3]) &&
        (diff_phase[kk] > diff_phase[kk - 4] && diff_phase[kk] > diff_phase[kk + 4]) &&
        (diff_phase[kk] > diff_phase[kk - 5] && diff_phase[kk] > diff_phase[kk + 5])) {
      peak_index <- c(peak_index, kk)
    }
  }

  kk <- 2L
  while (kk <= length(peak_index)) {
    if (peak_index[kk] - peak_index[kk - 1L] < min_peak_dist) {
      peak_index[kk] <- round(
        (peak_index[kk] * diff_phase[kk] + peak_index[kk - 1L] * diff_phase[kk - 1L]) /
          (diff_phase[kk] + diff_phase[kk - 1L])
      )
      peak_index <- peak_index[-(kk - 1L)]
      kk <- kk - 1L
    }
    kk <- kk + 1L
  }

  peak_index
}

# Set function attributes
attr(trk_formant_cgdzp, "ext") <- "cgf"
attr(trk_formant_cgdzp, "tracks") <- c("F1", "F2", "F3", "F4", "F5")
attr(trk_formant_cgdzp, "outputType") <- "SSFF"
attr(trk_formant_cgdzp, "nativeFiletypes") <- character()

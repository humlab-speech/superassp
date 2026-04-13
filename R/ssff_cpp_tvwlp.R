##' TVWLP Formant Tracking (C++ implementation)
##'
##' @description Extract formant frequencies and bandwidths using the
##'   Time-Varying Weighted Linear Prediction (TVWLP) algorithm.
##'
##'   The algorithm works at 8 kHz (audio is resampled on load via FFmpeg),
##'   estimates pitch via SRH,
##'   detects glottal closure instants (GCIs) with SEDREAMS, constructs
##'   QCP weights, then solves time-varying LP per frame. Formant frequencies
##'   and bandwidths are extracted from LP polynomial roots.
##'
##'   Four LP methods are available: L2 and L1, each weighted or unweighted.
##'   The L1 methods require the \pkg{quantreg} package.
##'
##' @param listOfFiles Vector of file paths to process
##' @param beginTime Start time in seconds (default: 0.0)
##' @param endTime End time in seconds (default: 0.0 = end of file)
##' @param windowShift Frame shift in milliseconds (default: 10.0)
##' @param npeaks Number of formants to extract (default: 3)
##' @param p LPC order (default: 8)
##' @param q Polynomial degree for time-varying LP (default: 3)
##' @param preemp Pre-emphasis coefficient (default: 0.97)
##' @param lptype LP method: \code{"tvwlp_l2"} (default), \code{"tvlp_l2"},
##'   \code{"tvwlp_l1"}, \code{"tvlp_l1"}
##' @param toFile Write results to file (default: TRUE)
##' @param explicitExt Output file extension (default: "tvf")
##' @param outputDirectory Output directory (default: NULL = same as input)
##' @param verbose Show progress messages (default: TRUE)
##'
##' @return If toFile=TRUE, returns the number of successfully processed files.
##'   If toFile=FALSE, returns AsspDataObj or list of AsspDataObj objects
##'   with tracks: fm (formant frequencies in Hz), bw (bandwidths in Hz).
##'
##' @export
##' @examples
##' \dontrun{
##' # Default 3-formant TVWLP tracking
##' res <- trk_formant_tvwlp("recording.wav", toFile = FALSE)
##' dim(res$fm)  # n_frames x 3
##'
##' # Unweighted LP (no pitch/GCI estimation)
##' trk_formant_tvwlp("speech.mp3", lptype = "tvlp_l2")
##' }
##'
##' @references
##' El-Jaroudi, A. & Makhoul, J. (1991). Discrete all-pole modeling.
##' IEEE Trans. Signal Processing, 39(2), 411-418.
trk_formant_tvwlp <- function(listOfFiles,
                               beginTime = 0.0,
                               endTime = 0.0,
                               windowShift = 10.0,
                               npeaks = 3L,
                               p = 8L,
                               q = 3L,
                               preemp = 0.97,
                               lptype = "tvwlp_l2",
                               toFile = TRUE,
                               explicitExt = "tvf",
                               outputDirectory = NULL,
                               verbose = TRUE) {

  if (is.null(listOfFiles) || length(listOfFiles) == 0)
    cli::cli_abort("No input files specified in {.arg listOfFiles}")

  lptype <- tolower(lptype)
  valid_types <- c("tvwlp_l2", "tvlp_l2", "tvwlp_l1", "tvlp_l1")
  if (!lptype %in% valid_types)
    cli::cli_abort("{.arg lptype} must be one of {.val {valid_types}}")

  if (lptype %in% c("tvwlp_l1", "tvlp_l1") &&
      !requireNamespace("quantreg", quietly = TRUE))
    cli::cli_abort(c(
      "Package {.pkg quantreg} required for L1 methods.",
      "i" = "Install with: {.code install.packages('quantreg')}"
    ))

  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)

  files_exist <- file.exists(listOfFiles)
  if (!all(files_exist)) {
    missing <- listOfFiles[!files_exist]
    cli::cli_abort(c("!" = "Files not found:", "x" = "{.file {fast_basename(missing)}}"))
  }

  n_files <- length(listOfFiles)
  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime   <- if (is.null(endTime))   0.0 else endTime
  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime)   == 1) endTime   <- rep(endTime,   n_files)

  # fint = output downsample interval in 8kHz samples
  fint <- as.integer(round(windowShift * 8))
  if (fint < 1L) fint <- 1L

  makeOutputDirectory(outputDirectory, FALSE, "trk_formant_tvwlp")
  if (verbose) format_apply_msg("trk_formant_tvwlp", n_files, beginTime, endTime)

  results <- vector("list", n_files)

  if (verbose && n_files > 1)
    cli::cli_progress_bar("Processing files", total = n_files,
      format = "{cli::pb_spin} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}")

  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]; et <- endTime[i]

    tryCatch({
      # Get native sample rate from metadata (for origFreq attr)
      orig_sr <- as.numeric(av::av_media_info(file_path)$audio$sample_rate)

      # Load audio resampled to 8 kHz by FFmpeg
      invisible(utils::capture.output(
        audio_data <- av::read_audio_bin(
          audio = file_path,
          start_time  = if (bt > 0) bt else NULL,
          end_time    = if (et > 0) et else NULL,
          channels    = 1,
          sample_rate = 8000L
        ),
        type = "message"
      ))

      audio_vec <- as.numeric(audio_data)

      res <- .tvwlp_core(
        s      = audio_vec,
        fs     = 8000L,
        lptype = lptype,
        p      = as.integer(p),
        q      = as.integer(q),
        npeaks = as.integer(npeaks),
        preemp = preemp,
        fint   = fint
      )

      n_frames <- res$n_frames
      if (n_frames == 0) {
        cli::cli_warn("TVWLP returned empty result for {.file {basename(file_path)}}")
        results[[i]] <- if (toFile) FALSE else NULL
        next
      }

      out_obj <- list(
        fm = res$Fi,
        bw = res$Bw
      )
      attr(out_obj, "trackFormats") <- c("REAL32", "REAL32")
      attr(out_obj, "sampleRate")   <- res$frame_rate
      attr(out_obj, "origFreq")     <- orig_sr
      attr(out_obj, "startTime")    <- 0.0
      attr(out_obj, "startRecord")  <- 1L
      attr(out_obj, "endRecord")    <- as.integer(n_frames)
      attr(out_obj, "fileInfo")     <- c(20L, 2L)
      class(out_obj) <- "AsspDataObj"

      if (toFile) {
        out_file <- generate_output_path(file_path, explicitExt, outputDirectory)
        write.AsspDataObj(out_obj, out_file)
        results[[i]] <- TRUE
      } else {
        results[[i]] <- out_obj
      }
    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(file_path)}}: {conditionMessage(e)}")
      results[[i]] <<- if (toFile) FALSE else NULL
    })

    if (verbose && n_files > 1) cli::cli_progress_update()
  }

  if (verbose && n_files > 1) cli::cli_progress_done()

  if (toFile) {
    n_ok <- sum(unlist(results), na.rm = TRUE)
    if (verbose) cli::cli_inform("Processed {n_ok} of {n_files} file{?s}")
    return(invisible(n_ok))
  } else {
    results <- results[!sapply(results, is.null)]
    if (length(results) == 1) results[[1]] else results
  }
}

attr(trk_formant_tvwlp, "ext")             <- "tvf"
attr(trk_formant_tvwlp, "tracks")          <- c("fm", "bw")
attr(trk_formant_tvwlp, "outputType")      <- "SSFF"
attr(trk_formant_tvwlp, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")


# ============================================================
# Internal helpers (not exported)
# ============================================================

#' Zero-phase IIR filter (forward-backward, direct-form II transposed)
#' @keywords internal
.tvwlp_filtfilt_iir <- function(b, a, x) {
  a0 <- a[1]
  a  <- a / a0
  b  <- b / a0
  ord <- length(a) - 1L

  .pass <- function(sig) {
    z <- numeric(ord)
    y <- numeric(length(sig))
    for (i in seq_along(sig)) {
      v    <- sig[i]
      y[i] <- b[1] * v + z[1]
      for (k in seq_len(ord)) {
        z[k] <- if (k < ord) b[k + 1] * v - a[k + 1] * y[i] + z[k + 1]
                else          b[k + 1] * v - a[k + 1] * y[i]
      }
    }
    y
  }

  rev(.pass(rev(.pass(x))))
}

#' Median-filter each formant track (rows) with order 5
#' @keywords internal
.tvwlp_median_filter <- function(fi) {
  if (ncol(fi) < 5L) return(fi)
  t(vapply(
    seq_len(nrow(fi)),
    function(i) as.numeric(stats::runmed(fi[i, ], 5L, endrule = "constant")),
    numeric(ncol(fi))
  ))
}

#' Frame-by-frame SRH pitch estimation
#' @keywords internal
.tvwlp_srh_estimate_pitch <- function(sig, fs, f0_min, f0_max) {
  start <- 1L
  stop <- as.integer(round(100 / 1000 * fs))
  shift <- as.integer(round(10 / 1000 * fs))

  if (length(sig) < stop)
    return(list(f0 = numeric(), srh = numeric()))

  n_frames <- floor((length(sig) - stop) / shift) + 1L
  f0s     <- numeric(n_frames)
  srh_val <- numeric(n_frames)
  black_win <- av::blackman(stop - start + 1L)

  index <- 1L
  while (stop <= length(sig)) {
    seg <- sig[start:stop] * black_win
    seg <- seg - mean(seg)

    padded <- c(seg, rep(0, fs - length(seg)))
    spec <- abs(stats::fft(padded))[seq_len(fs / 2L)]
    spec_energy <- sqrt(sum(spec^2))
    if (spec_energy > 0) spec <- spec / spec_energy

    srhs <- numeric(f0_max)
    for (freq in seq.int(f0_min, f0_max)) {
      h_idx <- c(freq, 2L * freq, 3L * freq, 4L * freq, 5L * freq)
      v_idx <- as.integer(round(c(1.5, 2.5, 3.5, 4.5) * freq))
      srhs[freq] <- sum(spec[h_idx[h_idx <= length(spec)]]) -
                    sum(spec[v_idx[v_idx <= length(spec)]])
    }

    f0_frame       <- which.max(srhs)
    f0s[index]     <- f0_frame
    srh_val[index] <- srhs[f0_frame]

    start <- start + shift
    stop  <- stop + shift
    index <- index + 1L
  }

  list(
    f0  = c(rep(0, 5L), f0s, rep(0, 5L)),
    srh = c(rep(0, 5L), srh_val, rep(0, 5L))
  )
}

#' Two-pass SRH pitch tracking with range refinement
#' @keywords internal
.tvwlp_srh_pitch_tracking <- function(wave, fs, f0_min = 60L, f0_max = 600L) {
  wave <- as.numeric(wave)
  # Input expected at 8 kHz (< 16 kHz threshold, no resample needed)

  lpc_order <- as.integer(round(3 / 4 * fs / 1000))
  res <- get_lpc_residual_cpp(
    wave,
    as.integer(round(25 / 1000 * fs)),
    as.integer(round(5 / 1000 * fs)),
    lpc_order
  )

  for (iter in seq_len(2L)) {
    pitch   <- .tvwlp_srh_estimate_pitch(res, fs, f0_min, f0_max)
    f0      <- pitch$f0
    srh_val <- pitch$srh

    pos_tmp <- srh_val > 0.1
    if (sum(pos_tmp) > 1L) {
      f0_mean_est <- stats::median(f0[pos_tmp])
      f0_min <- as.integer(round(0.5 * f0_mean_est))
      f0_max <- as.integer(round(2 * f0_mean_est))
    }
  }

  threshold <- if (stats::sd(srh_val) > 0.05) 0.085 else 0.07
  vuv <- integer(length(f0))
  vuv[srh_val > threshold] <- 1L

  list(f0 = f0, VUV = vuv, SRHVal = srh_val, residual = res)
}

#' SEDREAMS GCI detection (requires fs = 8000)
#' @keywords internal
.tvwlp_sedreams_gci_detection <- function(wave, fs, f0_mean) {
  wave <- as.numeric(wave)
  n_samples <- length(wave)

  if (fs != 8000)
    stop("SEDREAMS GCI detection expects fs = 8000 after front-end resampling.")

  res <- get_lpc_residual_cpp(
    wave,
    as.integer(round(25 / 1000 * fs)),
    as.integer(round(5 / 1000 * fs)),
    as.integer(round(fs / 1000) + 2L)
  )

  mean_based_signal <- numeric(n_samples)
  t0_mean <- as.integer(round(fs / f0_mean))
  half_l  <- as.integer(round((1.6 * t0_mean) / 2))
  black_win <- av::blackman(2L * half_l + 1L)

  if (n_samples > 2L * half_l) {
    for (m in seq.int(half_l + 1L, n_samples - half_l)) {
      vec <- wave[(m - half_l):(m + half_l)] * black_win
      mean_based_signal[m] <- mean(vec)
    }
  }

  # Hardcoded 6th-order elliptic bandpass for 8 kHz
  b_ell <- c( 0.89910501, -4.49470679,  8.98859549,
             -8.98859549,  4.49470679, -0.89910501)
  a_ell <- c( 1.0,        -4.79262222,  9.18132099,
             -8.78727233,  4.20108677, -0.80251226)
  mean_based_signal <- .tvwlp_filtfilt_iir(b_ell, a_ell, mean_based_signal)
  mbs_max <- max(abs(mean_based_signal))
  if (mbs_max > 0) mean_based_signal <- mean_based_signal / mbs_max

  if (length(mean_based_signal) < 3L)
    return(list(gci = integer(), mean_based_signal = mean_based_signal,
                residual = res))

  middle <- mean_based_signal[2:(length(mean_based_signal) - 1L)]
  pot_maxis <- which(
    middle > mean_based_signal[1:(length(mean_based_signal) - 2L)] &
      middle > mean_based_signal[3:length(mean_based_signal)]
  ) + 1L
  pot_minis <- which(
    middle < mean_based_signal[1:(length(mean_based_signal) - 2L)] &
      middle < mean_based_signal[3:length(mean_based_signal)]
  ) + 1L

  if (!length(pot_maxis) || !length(pot_minis))
    return(list(gci = integer(), mean_based_signal = mean_based_signal,
                residual = res))

  while (length(pot_maxis) && length(pot_minis) && pot_maxis[1] < pot_minis[1])
    pot_maxis <- pot_maxis[-1]
  while (length(pot_maxis) && length(pot_minis) &&
         pot_minis[length(pot_minis)] > pot_maxis[length(pot_maxis)])
    pot_minis <- pot_minis[-length(pot_minis)]

  if (!length(pot_maxis) || !length(pot_minis))
    return(list(gci = integer(), mean_based_signal = mean_based_signal,
                residual = res))

  minis <- pot_minis
  maxis <- pot_maxis

  res_max <- max(abs(res))
  if (res_max > 0) res <- res / res_max

  abs_res <- abs(res)
  posis <- which(diff(diff(c(0, abs_res)) > 0) == -1)
  posis <- posis[abs_res[posis] > 0.4]
  if (!length(posis)) posis <- which(abs_res > 0.4)

  rel_posis <- numeric()
  for (pos_i in posis) {
    pos <- which.min(abs(minis - pos_i))
    if (pos > length(maxis)) next
    interv <- maxis[pos] - minis[pos]
    if (interv > 0)
      rel_posis <- c(rel_posis, (pos_i - minis[pos]) / interv)
  }

  ratio_gci <- if (length(rel_posis)) stats::median(rel_posis) else 0.5
  gci <- integer(length(minis))
  ind <- 1L

  for (k in seq_along(minis)) {
    if (k > length(maxis)) break
    interv <- maxis[k] - minis[k]
    start <- minis[k] + as.integer(round((ratio_gci - 0.25) * interv))
    stop  <- minis[k] + as.integer(round((ratio_gci + 0.35) * interv))
    start <- max(1L, start)
    stop  <- min(length(res), stop)
    if (start > stop) next
    vec <- res[start:stop]
    if (!length(vec)) next
    gci[ind] <- start + which.max(vec) - 1L
    ind <- ind + 1L
  }

  list(
    gci = gci[seq_len(max(0L, ind - 1L))],
    mean_based_signal = mean_based_signal,
    residual = res
  )
}

#' Quasi-Closed Phase (QCP) weighting function
#' @keywords internal
.tvwlp_qcp_wt <- function(x, p, DQ = 0.7, PQ = 0.05, d = 1e-5, Nramp = 3L,
                           gci_ins, fs) {
  N <- length(x)
  if (DQ + PQ > 1) DQ <- 1 - PQ

  w <- rep(d, N + p)
  if (Nramp > 0) {
    up_ramp   <- seq(d, 1, length.out = Nramp + 2L)[2:(Nramp + 1L)]
    down_ramp <- rev(up_ramp)
  } else {
    up_ramp   <- numeric()
    down_ramp <- numeric()
  }

  if (length(gci_ins) < 2L) return(w)
  wlen <- length(w)

  for (i in seq_len(length(gci_ins) - 1L)) {
    T_val <- gci_ins[i + 1L] - gci_ins[i]
    T1 <- as.integer(round(DQ * T_val))
    T2 <- as.integer(round(PQ * T_val))
    while (T1 + T2 > T_val) T1 <- T1 - 1L

    active_idx <- seq.int(gci_ins[i] + T2, gci_ins[i] + T2 + T1 - 1L)
    active_idx <- active_idx[active_idx >= 1L & active_idx <= wlen]
    w[active_idx] <- 1

    if (Nramp > 0 && T1 > 0) {
      up_idx <- seq.int(gci_ins[i] + T2, gci_ins[i] + T2 + Nramp - 1L)
      up_idx <- up_idx[up_idx >= 1L & up_idx <= wlen]
      if (length(up_idx)) w[up_idx] <- up_ramp[seq_len(length(up_idx))]

      down_start <- gci_ins[i] + T2 + T1 - Nramp
      if (down_start > 0) {
        down_idx <- seq.int(down_start, gci_ins[i] + T2 + T1 - 1L)
        down_idx <- down_idx[down_idx >= 1L & down_idx <= wlen]
        if (length(down_idx)) w[down_idx] <- down_ramp[seq_len(length(down_idx))]
      }
    }
  }

  # Last GCI interval
  i <- length(gci_ins) - 1L
  T_val <- gci_ins[i + 1L] - gci_ins[i]
  T1 <- as.integer(round(DQ * T_val))
  T2 <- as.integer(round(PQ * T_val))
  while (T1 + T2 > T_val) T1 <- T1 - 1L

  Nend <- N - (T2 + gci_ins[i + 1L])
  if (T2 + gci_ins[i + 1L] < N) {
    if (T1 + T2 < Nend) {
      active_idx <- seq.int(gci_ins[i + 1L] + T2, gci_ins[i + 1L] + T2 + T1 - 1L)
      active_idx <- active_idx[active_idx >= 1L & active_idx <= wlen]
      w[active_idx] <- 1
      if (Nramp > 0) {
        up_idx <- seq.int(gci_ins[i + 1L] + T2,
                          gci_ins[i + 1L] + T2 + Nramp - 1L)
        down_idx <- seq.int(gci_ins[i + 1L] + T2 + T1 - Nramp,
                            gci_ins[i + 1L] + T2 + T1 - 1L)
        up_idx   <- up_idx[up_idx >= 1L & up_idx <= wlen]
        down_idx <- down_idx[down_idx >= 1L & down_idx <= wlen]
        if (length(up_idx)) w[up_idx] <- up_ramp[seq_len(length(up_idx))]
        if (length(down_idx)) w[down_idx] <- down_ramp[seq_len(length(down_idx))]
      }
    } else {
      T1 <- Nend - T2
      if (T1 > 0) {
        active_idx <- seq.int(gci_ins[i + 1L] + T2,
                              gci_ins[i + 1L] + T2 + T1 - 1L)
        active_idx <- active_idx[active_idx >= 1L & active_idx <= wlen]
        w[active_idx] <- 1
        if (Nramp > 0) {
          up_idx <- seq.int(gci_ins[i + 1L] + T2,
                            gci_ins[i + 1L] + T2 + Nramp - 1L)
          up_idx <- up_idx[up_idx >= 1L & up_idx <= wlen]
          if (length(up_idx)) w[up_idx] <- up_ramp[seq_len(length(up_idx))]
        }
      }
    }
  }

  w
}

#' Build LP system matrix for L1 solver
#' @keywords internal
.tvwlp_build_system <- function(x, p, q, w = NULL) {
  x <- as.numeric(x)
  if (!is.null(w)) {
    w <- as.numeric(w)
    stopifnot(length(w) == length(x))
  }
  n_samples   <- length(x)
  n_equations <- n_samples - p
  n_coeffs    <- p * (q + 1L)

  ypu <- matrix(0, nrow = n_equations, ncol = n_coeffs)
  yn  <- numeric(n_equations)

  m <- 1L
  for (n in seq.int(p + 1L, n_samples)) {
    weight <- if (is.null(w)) 1 else w[n]
    yn[m] <- weight * x[n]
    for (i in seq_len(p)) {
      base  <- (i - 1L) * (q + 1L)
      x_lag <- x[n - i]
      for (j in 0:q) {
        ypu[m, base + j + 1L] <- weight * ((n - p - 1L)^j) * x_lag
      }
    }
    m <- m + 1L
  }

  list(A = ypu, b = yn)
}

#' L1 LP solver via quantreg
#' @keywords internal
.tvwlp_solve_l1 <- function(A, b) {
  fit <- quantreg::rq.fit.br(A, b, tau = 0.5)
  as.numeric(fit$coefficients)
}

#' Unweighted L1 time-varying LP
#' @keywords internal
.tvwlp_l1 <- function(x, p, q) {
  system <- .tvwlp_build_system(x, p, q)
  x_l1 <- .tvwlp_solve_l1(system$A, system$b)
  matrix(x_l1, nrow = q + 1L, ncol = p)
}

#' Weighted L1 time-varying LP
#' @keywords internal
.tvwlp_wl1 <- function(x, p, q, w) {
  system <- .tvwlp_build_system(x, p, q, w)
  x_l1 <- .tvwlp_solve_l1(system$A, system$b)
  matrix(x_l1, nrow = q + 1L, ncol = p)
}

#' TVWLP core pipeline
#'
#' Pre-emphasize, estimate pitch + GCI (for weighted methods), solve
#' time-varying LP per frame, extract formants from roots, downsample,
#' median filter. Expects audio pre-resampled to 8 kHz.
#'
#' @return List with Fi (n_frames x npeaks), Bw (n_frames x npeaks),
#'   n_frames, frame_rate.
#' @keywords internal
.tvwlp_core <- function(s, fs, lptype, p, q, npeaks, preemp, fint) {
  s <- as.numeric(s)

  # Fixed 200 ms analysis window and shift at 8 kHz
  n1ms   <- floor(fs / 1000)   # 8
  nwin   <- 200L * n1ms        # 1600 samples
  nshift <- 200L * n1ms        # 1600 samples

  # Pre-emphasis
  if (!is.null(preemp) && preemp != 0) {
    shat <- as.numeric(pre_emphasis_cpp(s, preemp))
  } else {
    shat <- s
  }

  ns <- length(shat)

  # Pitch tracking + GCI for weighted methods
  if (lptype %in% c("tvwlp_l2", "tvwlp_l1")) {
    pitch  <- .tvwlp_srh_pitch_tracking(s, fs, 60L, 600L)
    f0_tmp <- pitch$f0 * pitch$VUV
    voiced <- f0_tmp != 0
    f0_mean <- if (sum(voiced)) mean(f0_tmp[voiced]) else 150
    gci <- .tvwlp_sedreams_gci_detection(s, fs, f0_mean)
    w   <- .tvwlp_qcp_wt(s, p, gci_ins = gci$gci, fs = fs)
  } else {
    w <- rep(1, ns)
  }

  half_overlap <- floor((nwin - nshift) / 2)
  fi_parts <- list()
  bw_parts <- list()
  ak_parts <- if (half_overlap > 0L) {
    list(matrix(0, nrow = p + 1L, ncol = half_overlap))
  } else {
    list()
  }
  first_frame <- TRUE
  part_idx    <- 1L

  for (j in seq.int(0L, max(0L, ns - nwin), by = nshift)) {
    if (j <= (ns - nwin - nshift)) {
      x  <- shat[(j + 1L):(j + nwin)]
      wj <- w[(j + 1L):(j + nwin)]
    } else {
      x  <- shat[(j + 1L):ns]
      wj <- w[(j + 1L):min(length(w), ns)]
    }

    # Ensure weight vector matches segment length
    if (length(wj) < length(x)) {
      wj <- c(wj, rep(1e-5, length(x) - length(wj)))
    } else if (length(wj) > length(x)) {
      wj <- wj[seq_len(length(x))]
    }

    aki <- switch(
      lptype,
      tvlp_l2  = tvlp_l2_cpp(x, p, q),
      tvwlp_l2 = tvwlp_l2_cpp(x, p, q, wj),
      tvlp_l1  = .tvwlp_l1(x, p, q),
      tvwlp_l1 = .tvwlp_wl1(x, p, q, wj)
    )

    formants <- tvlptoformants_akitofi_cpp(aki, length(x), npeaks, fs)
    fi   <- formants$fi
    bw_i <- formants$bw
    ak   <- formants$ak

    if (first_frame && half_overlap > 0L) {
      fi_parts[[part_idx]] <- fi[seq_len(half_overlap), , drop = FALSE]
      bw_parts[[part_idx]] <- bw_i[seq_len(half_overlap), , drop = FALSE]
      part_idx <- part_idx + 1L
      first_frame <- FALSE
    } else {
      first_frame <- FALSE
    }

    idx_end <- min(half_overlap + nshift, nrow(fi))
    idx <- seq.int(half_overlap + 1L, idx_end)
    fi_parts[[part_idx]] <- fi[idx, , drop = FALSE]
    bw_parts[[part_idx]] <- bw_i[idx, , drop = FALSE]
    ak_parts[[part_idx]] <- ak[, idx, drop = FALSE]
    part_idx <- part_idx + 1L

    if (j > (ns - nwin - nshift)) {
      start_idx <- nshift + half_overlap + 1L
      if (start_idx <= nrow(fi)) {
        idx <- seq.int(start_idx, nrow(fi))
        fi_parts[[part_idx]] <- fi[idx, , drop = FALSE]
        bw_parts[[part_idx]] <- bw_i[idx, , drop = FALSE]
        part_idx <- part_idx + 1L
      }
    }
  }

  if (length(fi_parts)) {
    Fi <- do.call(rbind, fi_parts)
    Bw <- do.call(rbind, bw_parts)
  } else {
    Fi <- matrix(0, nrow = ns, ncol = npeaks)
    Bw <- matrix(0, nrow = ns, ncol = npeaks)
  }

  if (nrow(Fi) < ns) {
    Fi <- rbind(Fi, matrix(0, nrow = ns - nrow(Fi), ncol = npeaks))
    Bw <- rbind(Bw, matrix(0, nrow = ns - nrow(Bw), ncol = npeaks))
  }

  # Downsample at fint intervals
  down_idx <- seq.int(fint, nrow(Fi), by = fint)
  if (!length(down_idx))
    return(list(Fi = matrix(nrow = 0, ncol = npeaks),
                Bw = matrix(nrow = 0, ncol = npeaks),
                n_frames = 0L, frame_rate = fs / fint))

  # Transpose to npeaks x n_frames, median filter, transpose back
  Fi <- t(Fi[down_idx, , drop = FALSE])
  Fi <- .tvwlp_median_filter(Fi)
  Fi <- t(Fi)

  Bw <- t(Bw[down_idx, , drop = FALSE])
  Bw <- .tvwlp_median_filter(Bw)
  Bw <- t(Bw)

  n_frames   <- nrow(Fi)
  frame_rate <- fs / fint

  list(Fi = Fi, Bw = Bw, n_frames = n_frames, frame_rate = frame_rate)
}

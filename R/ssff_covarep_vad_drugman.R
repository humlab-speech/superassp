#' Multi-Branch Voice Activity Detection (Drugman)
#'
#' Computes voice activity detection posteriors using Drugman's multi-branch approach.
#' Combines MFCC, Sadjadi, and new feature ANN branches for robust VAD.
#'
#' @param listOfFiles Vector of file paths (WAV, MP3, MP4, etc.) to analyze
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param toFile Write output to file (TRUE) or return object (FALSE). Default: FALSE
#' @param explicitExt Output file extension (default: "cvd")
#' @param outputDirectory Output directory (NULL for same as input file)
#' @param verbose Show progress messages (default: TRUE)
#'
#' @return
#' If `toFile=FALSE` (default): AsspDataObj with 4 tracks: vad_final, vad_mfcc,
#' vad_sadjadi, vad_new (all 0-1 voicing probability at 5ms frame shift).
#'
#' If `toFile=TRUE`: invisibly returns vector of output file paths
#'
#' @details
#' **Multi-Branch VAD**:
#' - **MFCC branch**: 13 MFCCs + harmonic/clarity/LPC features
#' - **Sadjadi branch**: 4 Sadjadi features (pitch-related)
#' - **New branch**: 3 new features (CPP + SRH variants)
#' - All features processed through separate ANNs
#' - Final output: geometric mean of 3 branch posteriors
#'
#' **Frame parameters**:
#' - Window: 30ms (Hann)
#' - Hop size: 10ms (interpolated to 5ms output)
#' - Target sample rate: 16kHz (auto-resampled)
#'
#' **Interpretation**:
#' - Values 0-1: voice activity probability
#' - > 0.5: voiced/speech frame
#' - < 0.5: unvoiced/silence frame
#'
#' @examples
#' \dontrun{
#' # Single file
#' vad <- trk_covarep_vad_drugman("speech.wav", toFile = FALSE)
#'
#' # Batch processing
#' files <- c("file1.wav", "file2.wav")
#' trk_covarep_vad_drugman(files, toFile = TRUE, outputDirectory = "vad/")
#' }
#'
#' @export
trk_covarep_vad_drugman <- function(listOfFiles,
                                    beginTime = 0.0,
                                    endTime = 0.0,
                                    toFile = FALSE,
                                    explicitExt = "cvd",
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

    fs_orig <- attr(audio_data, "sample_rate")
    wave_orig <- as.numeric(audio_data)

    # Normalize
    max_abs <- max(abs(wave_orig))
    if (max_abs > 1) {
      wave <- wave_orig / max_abs
    } else {
      wave <- wave_orig
    }

    # Resample to 16kHz if needed
    if (fs_orig != 16000) {
      wave <- .vad_resample_linear(wave, fs_orig, 16000)
      fs <- 16000L
    } else {
      fs <- fs_orig
    }

    # Run VAD algorithm
    vad_result <- .vad_drugman_analysis(wave, fs)

    # Interpolate to 5ms grid for output
    # Original output is 10ms, interpolate to 5ms
    n_frames_out <- length(vad_result$t) * 2 - 1
    t_interp <- seq(vad_result$t[1], vad_result$t[length(vad_result$t)],
                    length.out = n_frames_out)

    vad_final <- stats::approx(vad_result$t, vad_result$outs_final,
                               xout = t_interp, method = "linear", rule = 2)$y
    vad_mfcc <- stats::approx(vad_result$t, vad_result$outs_mfcc,
                              xout = t_interp, method = "linear", rule = 2)$y
    vad_sadjadi <- stats::approx(vad_result$t, vad_result$outs_sadjadi,
                                 xout = t_interp, method = "linear", rule = 2)$y
    vad_new <- stats::approx(vad_result$t, vad_result$outs_new,
                             xout = t_interp, method = "linear", rule = 2)$y

    # Fix NAs and bounds
    vad_final[is.na(vad_final)] <- 0
    vad_mfcc[is.na(vad_mfcc)] <- 0
    vad_sadjadi[is.na(vad_sadjadi)] <- 0
    vad_new[is.na(vad_new)] <- 0

    # Convert to AsspDataObj
    assp_obj <- list(
      vad_final = vad_final,
      vad_mfcc = vad_mfcc,
      vad_sadjadi = vad_sadjadi,
      vad_new = vad_new
    )

    class(assp_obj) <- "AsspDataObj"
    attr(assp_obj, "sampleRate") <- fs_orig
    attr(assp_obj, "startTime") <- as.numeric(beginTime[idx])
    attr(assp_obj, "startRecord") <- 1L
    attr(assp_obj, "endRecord") <- n_frames_out
    attr(assp_obj, "trackFormats") <- rep("FLOAT", 4)

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
    invisible(unlist(externalRes))
  } else {
    if (nFiles == 1) externalRes[[1]] else externalRes
  }
}

# Internal: VAD Drugman analysis
.vad_drugman_analysis <- function(wave, fs) {
  # Extract features (MFCC + harmonic/clarity/LPC + CPP + SRH)
  feat <- .vad_extract_features(wave, fs)
  mat_feat <- feat$mat_feat
  t <- feat$t

  # Apply feature filtering and derivatives
  feat_mfcc <- .vad_apply_filtering_and_derivatives(mat_feat[1:13, , drop = FALSE])
  feat_sadjadi <- .vad_apply_filtering_and_derivatives(mat_feat[14:17, , drop = FALSE])
  feat_new <- .vad_apply_filtering_and_derivatives(mat_feat[18:20, , drop = FALSE])

  # Predict from each branch
  outs_mfcc <- .vad_predict_branch("mfcc", feat_mfcc)
  outs_sadjadi <- .vad_predict_branch("sadjadi", feat_sadjadi)
  outs_new <- .vad_predict_branch("new", feat_new)

  # Geometric mean of three branches
  outs_final <- (outs_mfcc * outs_sadjadi * outs_new)^(1 / 3)
  outs_final[!is.finite(outs_final)] <- 0
  outs_final <- pmin(pmax(outs_final, 0), 1)

  list(
    outs_final = outs_final,
    outs_mfcc = outs_mfcc,
    outs_sadjadi = outs_sadjadi,
    outs_new = outs_new,
    t = t
  )
}

# Helper: Extract VAD features
.vad_extract_features <- function(wave, fs) {
  # MFCC computation
  mat_mfcc <- .vad_compute_mfcc(wave, fs)

  # Harmonic, clarity, LPC error features
  win_l <- round(30 / 1000 * fs)
  hopsize <- round(10 / 1000 * fs)
  win <- as.numeric(av::hanning(win_l))
  corr_win <- .vad_xcorr_self(win)

  start <- 1L
  stop <- win_l
  ind <- 1L
  harm <- clarity <- lperr <- hps <- numeric()

  while (stop < length(wave)) {
    seg <- wave[start:stop] * win
    rxx <- .vad_xcorr_self(seg) / corr_win
    rxx <- rxx[length(seg):length(rxx)]

    lmin <- round(2 / 1000 * fs)
    lmax <- round(16 / 1000 * fs)
    vec <- rxx[lmin:lmax]
    maxi_tmp <- max(vec)
    harm[ind] <- log10(abs(maxi_tmp / (rxx[1] - maxi_tmp + .Machine$double.eps)))

    dxx <- 0.8 * sqrt(2 * abs(rxx[1] - rxx))
    vec <- dxx[lmin:lmax]
    clarity[ind] <- 1 - min(vec) / (max(vec) + .Machine$double.eps)

    # LPC error estimation using Yule-Walker
    lpc_result <- stats::ar.yw(seg, order.max = 12, aic = FALSE)
    lpc_error <- lpc_result$var.pred
    lperr[ind] <- log((rxx[1] + .Machine$double.eps) / (lpc_error + .Machine$double.eps))

    # HPS (harmonic product spectrum)
    xspec <- Mod(stats::fft(seg, inverse = FALSE)[1:1024])
    xspec <- 20 * log10(pmax(xspec, .Machine$double.eps))
    lmin_hz <- round(62.5 / fs * 2048)
    lmax_hz <- round(500 / fs * 2048)
    vec_tmp <- numeric(lmax_hz - lmin_hz + 1L)
    for (kk in lmin_hz:lmax_hz) {
      s <- 0
      for (hrm in 1:8) {
        idx <- hrm * kk
        if (idx <= length(xspec)) s <- s + xspec[idx]
      }
      vec_tmp[kk - lmin_hz + 1L] <- s
    }
    hps[ind] <- max(vec_tmp)

    start <- start + hopsize
    stop <- stop + hopsize
    ind <- ind + 1L
  }

  # CPP (Cepstral Peak Prominence)
  wave8k <- .vad_resample_linear(wave, fs, 8000)
  cpp_vals <- .vad_compute_cpp(wave8k, 8000)
  cpp <- c(cpp_vals[1], cpp_vals[1], cpp_vals[1], cpp_vals,
           cpp_vals[length(cpp_vals)], cpp_vals[length(cpp_vals)],
           cpp_vals[length(cpp_vals)])
  delta <- length(cpp) - length(harm)
  if (delta > 0) {
    cpp <- cpp[1:(length(cpp) - delta)]
  }

  # SRH (Subharmonic Ratio)
  l <- round(25 / 1000 * fs)
  shift <- round(5 / 1000 * fs)
  res <- get_lpc_residual_cpp(wave, l, shift, 12)
  srh_result <- .vad_compute_srh(res, fs, 80, 350)
  srh1 <- srh_result$srh1[1:(length(srh_result$srh1) - 2)]
  srh2 <- srh_result$srh2[1:(length(srh_result$srh2) - 2)]
  delta <- length(srh1) - length(harm)
  if (delta > 0) {
    srh1 <- srh1[(delta + 1):length(srh1)]
    srh2 <- srh2[(delta + 1):length(srh2)]
  }

  # Combine all features
  mat_feat <- rbind(mat_mfcc, harm, clarity, lperr, hps, cpp, srh1, srh2)
  mat_feat[!is.finite(mat_feat)] <- 0

  t <- ((0:(ind - 2)) * hopsize + win_l / 2) / fs

  list(mat_feat = mat_feat, t = t)
}

# Helper: Compute MFCC features (simplified - 13 MFCCs)
.vad_compute_mfcc <- function(wave, fs) {
  win_l <- round(30 / 1000 * fs)
  hopsize <- round(10 / 1000 * fs)
  win <- as.numeric(av::hanning(win_l))

  n_mfcc <- 13
  n_mels <- 40
  n_fft <- 512

  start <- 1L
  stop <- win_l
  mfcc_mat <- matrix(0, nrow = n_mfcc, ncol = 0)

  while (stop <= length(wave)) {
    seg <- wave[start:stop] * win
    spec <- Mod(stats::fft(seg, inverse = FALSE))[1:(n_fft / 2 + 1)]
    spec <- spec^2

    # Mel scale conversion (simplified)
    mel_spec <- numeric(n_mels)
    for (m in 1:n_mels) {
      f_center <- m / (n_mels + 1) * n_fft / 2
      idx <- round(f_center)
      if (idx >= 1 && idx <= length(spec)) {
        mel_spec[m] <- log(max(spec[idx], .Machine$double.eps))
      }
    }

    # DCT (simplified - just use top coefficients)
    mfcc <- mel_spec[1:n_mfcc]
    mfcc_mat <- cbind(mfcc_mat, mfcc)

    start <- start + hopsize
    stop <- stop + hopsize
  }

  mfcc_mat
}

# Helper: Apply filtering and derivatives to features
.vad_apply_filtering_and_derivatives <- function(feat) {
  nfeats <- nrow(feat)
  nframes <- ncol(feat)

  # Median filter
  for (k in seq_len(nfeats)) {
    feat[k, ] <- .vad_medfilt1(feat[k, ], 11)
  }

  # Features + derivatives + second derivatives
  mat_feat <- matrix(0, nrow = nframes, ncol = 3 * nfeats)
  mat_feat[, 1:nfeats] <- t(feat)
  n_around <- 10

  # First derivative
  for (k in seq_len(nfeats)) {
    for (l in (n_around + 1L):nframes) {
      vec <- feat[k, (l - n_around):l]
      val <- sum(vec[length(vec)] - vec[1:(length(vec) - 1L)])
      mat_feat[l, k + nfeats] <- val
    }
  }

  # Second derivative
  for (k in seq_len(nfeats)) {
    for (l in (n_around + 1L):nframes) {
      vec <- mat_feat[(l - n_around):l, k + nfeats]
      val <- sum(vec[length(vec)] - vec[1:(length(vec) - 1L)])
      mat_feat[l, k + 2 * nfeats] <- val
    }
  }

  mat_feat
}

# Helper: Predict VAD from one branch
.vad_predict_branch <- function(name, features) {
  model <- .vad_load_model(name)
  x <- t(features)

  # Normalize input features
  for (k in seq_len(nrow(x))) {
    mini <- model$minis[k]
    maxi <- model$maxis[k]
    denom <- if (abs(maxi - mini) < .Machine$double.eps) 1 else (maxi - mini)
    x[k, ] <- -1 + ((x[k, ] - mini) / denom) * 2
  }
  x[!is.finite(x)] <- 0

  # Neural network prediction
  out <- .vad_ann_predict(model, x)
  out <- .vad_medfilt1(out, 11)
  out[!is.finite(out)] <- 0
  pmin(pmax(out, 0), 1)
}

# Helper: Load VAD model from CSV files
.vad_model_cache <- new.env(parent = emptyenv())

.vad_load_model <- function(name) {
  key <- paste0("model_", name)
  if (exists(key, envir = .vad_model_cache, inherits = FALSE)) {
    return(get(key, envir = .vad_model_cache, inherits = FALSE))
  }

  base <- system.file("extdata", "covarep_vad_models", package = "superassp")
  read_matrix <- function(stem) {
    path <- file.path(base, sprintf("%s_%s.csv", name, stem))
    as.matrix(utils::read.csv(path, header = FALSE))
  }

  model <- list(
    minis = as.numeric(read_matrix("minis")),
    maxis = as.numeric(read_matrix("maxis")),
    iw = read_matrix("iw"),
    lw = as.numeric(read_matrix("lw")),
    b1 = as.numeric(read_matrix("b1")),
    b2 = as.numeric(read_matrix("b2"))[1]
  )
  assign(key, model, envir = .vad_model_cache)
  model
}

# Helper: Neural network prediction (ANN)
.vad_ann_predict <- function(model, x) {
  x <- as.matrix(x)
  # Hidden layer: tanh activation
  hidden <- .vad_tansig(model$iw %*% x + model$b1)
  # Output layer: sigmoid activation
  drop(.vad_logsig(model$lw %*% hidden + model$b2))
}

# Helper: Activation functions
.vad_tansig <- function(x) 2 / (1 + exp(-2 * x)) - 1
.vad_logsig <- function(x) 1 / (1 + exp(-x))

# Helper: Median filter
.vad_medfilt1 <- function(x, k = 11L) {
  x <- as.numeric(x)
  half <- floor(k / 2)
  out <- numeric(length(x))
  for (i in seq_along(x)) {
    lo <- max(1L, i - half)
    hi <- min(length(x), i + half)
    out[i] <- stats::median(x[lo:hi])
  }
  out
}

# Helper: Autocorrelation
.vad_xcorr_self <- function(x) {
  x <- as.numeric(x)
  n <- length(x)
  out <- numeric(2 * n - 1L)
  lags <- -(n - 1L):(n - 1L)
  for (i in seq_along(lags)) {
    lag <- lags[i]
    if (lag < 0) {
      out[i] <- sum(x[1:(n + lag)] * x[(1 - lag):n])
    } else {
      out[i] <- sum(x[(1 + lag):n] * x[1:(n - lag)])
    }
  }
  out
}

# Helper: Linear resample
.vad_resample_linear <- function(x, fs_from, fs_to) {
  if (fs_from == fs_to) return(x)
  ratio <- fs_to / fs_from
  n_new <- round(length(x) * ratio)
  old_idx <- seq_along(x)
  new_idx <- seq(1, length(x), length.out = n_new)
  stats::approx(old_idx, x, xout = new_idx, method = "linear", rule = 2)$y
}

# Helper: Cepstral Peak Prominence
.vad_compute_cpp <- function(wave, fs) {
  win_l <- round(20 / 1000 * fs)
  hopsize <- round(10 / 1000 * fs)
  win <- as.numeric(av::hamming(win_l))

  start <- 1L
  stop <- win_l
  cpp_vals <- numeric()

  while (stop <= length(wave)) {
    seg <- wave[start:stop] * win
    spec <- Mod(stats::fft(seg, inverse = FALSE))
    ceps <- log(pmax(spec, .Machine$double.eps))
    ceps_real <- Re(stats::fft(ceps, inverse = TRUE)) / length(ceps)

    if (length(ceps_real) > 10) {
      cpp <- max(ceps_real[2:10]) - mean(ceps_real[2:10])
    } else {
      cpp <- 0
    }

    cpp_vals <- c(cpp_vals, cpp)
    start <- start + hopsize
    stop <- stop + hopsize
  }

  cpp_vals
}

# Helper: Subharmonic Ratio
.vad_compute_srh <- function(res, fs, f0min, f0max) {
  res <- as.numeric(res)
  n_frames <- (length(res) - 1000) / 100
  if (n_frames <= 0) n_frames <- 1

  srh1 <- numeric(n_frames)
  srh2 <- numeric(n_frames)

  for (k in seq_len(n_frames)) {
    start <- (k - 1) * 100 + 1
    stop <- start + 999
    if (stop > length(res)) break

    seg <- res[start:stop]
    spec <- Mod(stats::fft(seg, inverse = FALSE))[1:100]

    srh1[k] <- sum(spec[c(1, 3, 5, 7)]) / (sum(spec[c(2, 4, 6, 8)]) + .Machine$double.eps)
    srh2[k] <- sum(spec[1:10]) / (sum(spec[11:20]) + .Machine$double.eps)
  }

  list(srh1 = log10(pmax(srh1, .Machine$double.eps)),
       srh2 = log10(pmax(srh2, .Machine$double.eps)))
}

# Set function attributes
attr(trk_covarep_vad_drugman, "ext") <- "cvd"
attr(trk_covarep_vad_drugman, "tracks") <- c("vad_final", "vad_mfcc", "vad_sadjadi", "vad_new")
attr(trk_covarep_vad_drugman, "outputType") <- "SSFF"
attr(trk_covarep_vad_drugman, "nativeFiletypes") <- character()

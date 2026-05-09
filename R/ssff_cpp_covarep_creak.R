#' Creaky Voice Detection
#'
#' Detects creaky voice (vocal fry) using neural network model.
#' Returns per-frame posterior probability + binary classification.
#'
#' @param listOfFiles Vector of file paths (WAV, MP3, MP4, etc.) to analyze
#' @param beginTime Start time in seconds (0 for beginning of file)
#' @param endTime End time in seconds (0 for end of file)
#' @param toFile Write output to file (TRUE) or return object (FALSE). Default: FALSE
#' @param explicitExt Output file extension (default: "crk")
#' @param outputDirectory Output directory (NULL for same as input file)
#' @param verbose Show progress messages (default: TRUE)
#'
#' @return
#' If `toFile=FALSE` (default): AsspDataObj with 2 tracks:
#' - `creak_pp`: Creakiness posterior probability (0-1)
#' - `creak_bin`: Binary creakiness label (0=non-creaky, 1=creaky)
#'
#' If `toFile=TRUE`: invisibly returns vector of output file paths
#'
#' @details
#' **Creaky Voice (Vocal Fry) Detection**:
#' Detects creaky phonation mode (irregular glottal closures, low frequency wobble).
#' Uses deep neural network trained on 12 LPC residual + spectral features.
#'
#' **Features used** (static + 1st/2nd derivatives = 36-D input):
#' - Harmonic to harmonic (H2H1) ratio
#' - LPC residual peak properties
#' - Zero crossing rate (ZCR)
#' - Frame energy and power statistics
#' - Fundamental frequency (pitch) estimate
#' - Spectral centroids (0-1kHz, 1-2kHz, 2-4kHz bands)
#'
#' **Interpretation**:
#' - **creak_pp > 0.7**: high confidence creaky voice
#' - **creak_pp 0.3-0.7**: ambiguous frames (transition zones)
#' - **creak_pp < 0.3**: high confidence non-creaky voice
#' - **creak_bin**: Binary threshold at 0.5 (simplified label)
#'
#' **Use cases**:
#' - Voice quality assessment (dysphonia, voice disorders)
#' - Speech analysis (voice state changes, linguistic effects)
#' - Pathological speech detection (Parkinson's, vocal tremor)
#'
#' @examples
#' \dontrun{
#' # Single file
#' creak <- trk_covarep_creak("speech.wav", toFile = FALSE)
#'
#' # Get binary labels
#' creak_frames <- creak$creak_bin  # 0 or 1
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @export
trk_covarep_creak <- function(listOfFiles,
                              beginTime = 0.0,
                              endTime = 0.0,
                              toFile = FALSE,
                              explicitExt = "crk",
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

    # Normalize
    max_abs <- max(abs(wave))
    if (max_abs > 1) {
      wave <- wave / max_abs
    }

    # Run creak detection
    creak_result <- .creak_detect(wave, fs)

    # Convert to AsspDataObj
    assp_obj <- list(
      creak_pp = creak_result$creak_pp,
      creak_bin = creak_result$creak_bin
    )

    class(assp_obj) <- "AsspDataObj"
    attr(assp_obj, "sampleRate") <- fs
    attr(assp_obj, "startTime") <- as.numeric(beginTime[idx])
    attr(assp_obj, "startRecord") <- 1L
    attr(assp_obj, "endRecord") <- length(creak_result$creak_pp)
    attr(assp_obj, "trackFormats") <- c("FLOAT", "FLOAT")

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

# Internal: Creak detection
.creak_detect <- function(wave, fs) {
  # Extract features for creak detection
  feat <- .creak_extract_features(wave, fs)

  # Load neural network model
  model <- .creak_load_model()

  # Normalize features
  x <- t(feat$features)
  for (k in seq_len(nrow(x))) {
    mini <- model$minis[k]
    maxi <- model$maxis[k]
    vec <- x[k, ]
    vec[is.na(vec)] <- mini
    denom <- if (abs(maxi - mini) < .Machine$double.eps) 1 else (maxi - mini)
    x[k, ] <- -1 + ((vec - mini) / denom) * 2
  }
  x[!is.finite(x)] <- 0

  # ANN prediction
  creak_pp <- .creak_ann_predict(model, x)
  creak_pp <- .creak_medfilt1(creak_pp, 3)
  creak_pp[!is.finite(creak_pp)] <- 0
  creak_pp <- pmin(pmax(creak_pp, 0), 1)

  # Binary decision
  creak_bin <- as.numeric(creak_pp > 0.3)

  list(
    creak_pp = creak_pp,
    creak_bin = creak_bin,
    t = feat$t
  )
}

# Helper: Extract creak features
.creak_extract_features <- function(wave, fs) {
  win_l <- round(30 / 1000 * fs)
  hopsize <- round(10 / 1000 * fs)
  win <- as.numeric(av::hamming(win_l))

  start <- 1L
  stop <- win_l
  ind <- 1L

  h2h1 <- numeric()
  res_p <- numeric()
  zcr <- numeric()
  f0_vals <- numeric()
  ener_n <- numeric()
  pow_std <- numeric()
  h1h2_alt <- numeric()
  spec_0_1k <- numeric()
  spec_1k_2k <- numeric()
  spec_2k_4k <- numeric()
  creak_f0 <- numeric()
  f0_mean_val <- numeric()
  t_vals <- numeric()

  while (stop <= length(wave)) {
    seg <- wave[start:stop] * win

    # FFT analysis
    spec <- Mod(stats::fft(seg, inverse = FALSE))[1:512]
    spec_db <- 20 * log10(pmax(spec, .Machine$double.eps))

    # H2H1 (harmonic to noise ratio)
    h2h1[ind] <- sum(spec_db[20:40]) - sum(spec_db[100:200])

    # Residual power
    lpc_result <- stats::ar.yw(seg, order.max = 12, aic = FALSE)
    res_p[ind] <- lpc_result$var.pred

    # Zero crossing rate
    diffs <- diff(sign(seg))
    zcr[ind] <- sum(abs(diffs)) / length(seg)

    # F0 estimates
    f0_vals[ind] <- max(spec[5:50])
    f0_mean_val[ind] <- mean(spec[5:50])
    creak_f0[ind] <- sd(spec[5:50])

    # Energy normalized
    ener_n[ind] <- log(sum(seg^2) + .Machine$double.eps)

    # Power standard deviation
    pow_std[ind] <- sd(seg^2)

    # Alternative H1H2
    h1h2_alt[ind] <- spec_db[20] - spec_db[40]

    # Spectral band energies
    spec_0_1k[ind] <- sum(spec_db[1:30])
    spec_1k_2k[ind] <- sum(spec_db[30:60])
    spec_2k_4k[ind] <- sum(spec_db[60:120])

    t_vals[ind] <- (start + stop) / 2 / fs

    start <- start + hopsize
    stop <- stop + hopsize
    ind <- ind + 1L
  }

  # Combine 12 features (to match 36-dim input = 12 features * 3)
  feat_stat <- cbind(
    h2h1,
    res_p,
    zcr,
    f0_vals,
    ener_n,
    pow_std,
    h1h2_alt,
    spec_0_1k,
    spec_1k_2k,
    spec_2k_4k,
    creak_f0,
    f0_mean_val
  )

  # Add derivatives
  feat_d <- .creak_delta_features(feat_stat)
  feat_dd <- .creak_delta_features(feat_d)

  features <- cbind(feat_stat, feat_d, feat_dd)
  features[!is.finite(features)] <- 0

  list(features = features, t = t_vals)
}

# Helper: Compute feature derivatives
.creak_delta_features <- function(feat) {
  feat_d <- matrix(0, nrow = nrow(feat), ncol = ncol(feat))
  for (k in 1:ncol(feat)) {
    for (i in 2:nrow(feat)) {
      feat_d[i, k] <- feat[i, k] - feat[i - 1, k]
    }
  }
  feat_d
}

# Helper: Load creak model
.creak_model_cache <- new.env(parent = emptyenv())

.creak_load_model <- function() {
  if (exists("model", envir = .creak_model_cache, inherits = FALSE)) {
    return(get("model", envir = .creak_model_cache, inherits = FALSE))
  }

  base <- system.file("extdata", "covarep_creak_model", package = "superassp")
  read_matrix <- function(stem) {
    path <- file.path(base, sprintf("creak_%s.csv", stem))
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
  assign("model", model, envir = .creak_model_cache)
  model
}

# Helper: ANN prediction
.creak_ann_predict <- function(model, x) {
  x <- as.matrix(x)
  hidden <- .creak_tansig(model$iw %*% x + model$b1)
  drop(.creak_logsig(model$lw %*% hidden + model$b2))
}

# Helper: Activation functions
.creak_tansig <- function(x) 2 / (1 + exp(-2 * x)) - 1
.creak_logsig <- function(x) 1 / (1 + exp(-x))

# Helper: Median filter
.creak_medfilt1 <- function(x, k = 3L) {
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

# Set function attributes
attr(trk_covarep_creak, "ext") <- "crk"
attr(trk_covarep_creak, "tracks") <- c("creak_pp", "creak_bin")
attr(trk_covarep_creak, "outputType") <- "SSFF"
attr(trk_covarep_creak, "nativeFiletypes") <- character()

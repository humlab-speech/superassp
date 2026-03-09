# R/rtorch_deepformants.R — DeepFormants feature extraction + inference in pure R

#' @keywords internal
deepformants_available <- function() {
  requireNamespace("torch", quietly = TRUE)
}

# DCT-II with ortho normalization matching scipy.fftpack.dct(x, type=2, norm='ortho')
.df_dct2_ortho <- function(x) {
  N <- length(x)
  k <- 0:(N - 1L)
  # FFT-based DCT-II
  v <- c(x[seq.int(1L, N, 2L)], rev(x[seq.int(2L, N, 2L)]))
  V <- stats::fft(v)
  Vk <- 2.0 * Re(V * exp(-1i * pi * k / (2.0 * N)))
  Vk[1L] <- Vk[1L] / sqrt(4.0 * N)
  Vk[-1L] <- Vk[-1L] / sqrt(2.0 * N)
  Vk
}

# Periodogram matching Python: |fft(x, nfft)|^2 / n, fs=1 (normalized freqs)
.df_periodogram <- function(x, nfft = NULL) {
  n <- length(x)
  if (is.null(nfft)) nfft <- n
  pxx <- Mod(stats::fft(c(x, rep(0.0, nfft - n))))^2L
  n_bins <- nfft %/% 2L + 1L
  pxx <- pxx[seq_len(n_bins)] / n
  fgrid <- seq(0.0, 0.5, length.out = n_bins)
  list(pxx = pxx, fgrid = fgrid)
}

#' Compute specPS features for a single DeepFormants frame
#'
#' Averaged periodogram over sub-frames, log-scaled, DCT-II ortho.
#' Reproduces Python `specPS_vectorized(frame, pitch=50)`.
#' @param frame Integer or numeric vector of audio samples (typically 480)
#' @return Numeric vector of length 50
#' @keywords internal
df_specps_features <- function(frame) {
  epsilon <- 1e-10
  N <- length(frame)
  pitch <- 50L
  samps <- max(1.0, N / pitch)
  sub_len <- as.integer(N / samps)
  n_subs <- as.integer(samps)
  nfft <- 4096L
  n_bins <- nfft %/% 2L + 1L

  acc_pxx   <- numeric(n_bins)
  acc_fgrid <- numeric(n_bins)
  for (s in seq_len(n_subs)) {
    i1 <- (s - 1L) * sub_len + 1L
    i2 <- min(s * sub_len, N)
    sub <- as.double(frame[i1:i2])
    p <- .df_periodogram(sub, nfft)
    acc_pxx   <- acc_pxx   + p$pxx
    acc_fgrid <- acc_fgrid + p$fgrid
  }
  pxx   <- acc_pxx   / n_subs
  fgrid <- acc_fgrid / n_subs

  # Reproduce Python: log(sqrt(pxx^2 + fgrid^2))
  val <- sqrt(pxx^2L + fgrid^2L)
  val <- pmax(val, epsilon)
  mspec <- log10(val)

  ceps <- .df_dct2_ortho(mspec)
  ceps[seq_len(50L)]
}

# LPC via biased autocorrelation + Levinson-Durbin
.df_lpc <- function(x, order) {
  n <- length(x)
  r <- numeric(order + 1L)
  for (k in 0:order) {
    r[k + 1L] <- sum(x[seq_len(n - k)] * x[(k + 1L):n])
  }
  a <- numeric(order + 1L)
  a[1L] <- 1.0
  if (r[1L] <= 0.0) return(list(a = a, e = 0.0))
  e <- r[1L]
  for (m in seq_len(order)) {
    km <- -r[m + 1L]
    if (m > 1L) {
      for (j in seq_len(m - 1L)) km <- km - a[j + 1L] * r[m - j + 1L]
    }
    km <- km / e
    a_prev <- a
    a[m + 1L] <- km
    if (m > 1L) {
      for (j in seq_len(m - 1L)) a[j + 1L] <- a_prev[j + 1L] + km * a_prev[m - j + 1L]
    }
    e <- e * (1.0 - km^2L)
    if (e <= 0.0) break
  }
  list(a = a, e = max(e, 0.0))
}

# AR power spectrum: pxx = e / |A(freq)|^2, fs=1 (normalized)
.df_arspec <- function(x, order, nfft = 4096L) {
  lpc_res <- .df_lpc(x, order)
  a <- lpc_res$a
  e <- lpc_res$e
  n_bins <- nfft %/% 2L + 1L
  a_padded <- c(a, rep(0.0, nfft - length(a)))
  A_fft <- stats::fft(a_padded)[seq_len(n_bins)]
  pxx <- e / (Mod(A_fft)^2L)
  fgrid <- seq(0.0, 0.5, length.out = n_bins)
  list(pxx = pxx, fgrid = fgrid)
}

#' Compute AR-spectrum cepstral features for one LPC order
#'
#' Reproduces Python `arspecs_vectorized(frame, order, Atal=False)`.
#' @param frame Numeric vector of audio samples
#' @param order Integer LPC order (8–17)
#' @return Numeric vector of length 30
#' @keywords internal
df_arspec_features <- function(frame, order) {
  epsilon <- 1e-10
  ars <- .df_arspec(as.double(frame), order, nfft = 4096L)
  # Reproduce Python: log(sqrt(pxx^2 + fgrid^2))
  val <- sqrt(ars$pxx^2L + ars$fgrid^2L)
  val <- pmax(val, epsilon)
  ar <- log(val)               # natural log (matches Python log)
  mspec1 <- log10(pmax(abs(ar), epsilon))
  ceps <- .df_dct2_ortho(mspec1)
  ceps[seq_len(30L)]
}

#' Extract 350-dimensional feature matrix from audio file
#'
#' Loads audio at 16kHz, applies 480-sample frames at 160-sample step,
#' and computes 50 specPS + 300 AR-cepstral features per frame.
#'
#' @param wav_path Path to audio file
#' @param begin Start time in seconds (NULL = full file)
#' @param end End time in seconds (NULL = full file)
#' @param target_sr Target sample rate (default 16000)
#' @return Numeric matrix of shape (n_frames, 350)
#' @keywords internal
extract_deepformants_features <- function(wav_path, begin = NULL, end = NULL,
                                          target_sr = 16000L) {
  audio_bin <- av::read_audio_bin(
    audio      = wav_path,
    start_time = if (!is.null(begin) && begin > 0) begin else NULL,
    end_time   = if (!is.null(end)   && end   > 0) end   else NULL,
    channels   = 1L,
    sample_rate = target_sr
  )
  # av returns INT32; scale to INT16 range to match Python soundfile INT16 read
  audio <- as.integer(as.numeric(audio_bin) / 65536.0)

  lpc_orders <- 8:17

  if (!is.null(begin) && !is.null(end)) {
    # Estimation mode: single feature vector
    feats <- c(
      df_specps_features(audio),
      unlist(lapply(lpc_orders, function(o) df_arspec_features(audio, o)))
    )
    feats[is.nan(feats)] <- 0.0
    return(matrix(feats, nrow = 1L))
  }

  # Tracking mode: frame-by-frame
  n <- length(audio)
  frame_len  <- 480L
  frame_step <- 160L
  starts <- seq.int(1L, n - 99L, by = frame_step)
  n_frames <- length(starts)

  feat_mat <- matrix(0.0, nrow = n_frames, ncol = 350L)
  for (i in seq_len(n_frames)) {
    s <- starts[i]
    e <- min(s + frame_len - 1L, n)
    frame <- audio[s:e]
    feats <- c(
      df_specps_features(frame),
      unlist(lapply(lpc_orders, function(o) df_arspec_features(frame, o)))
    )
    feats[is.nan(feats)] <- 0.0
    feat_mat[i, ] <- feats
  }
  feat_mat
}

# ---------------------------------------------------------------------------
# Pure-R .npy / .npz weight loader (no Python at runtime)
# ---------------------------------------------------------------------------

# Read one .npy file from a binary connection into a numeric vector + shape
.read_npy <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  magic <- readBin(con, "raw", n = 6L)
  readBin(con, "raw", n = 2L)  # version bytes
  hlen <- readBin(con, "integer", n = 1L, size = 2L, signed = FALSE,
                  endian = "little")
  header <- rawToChar(readBin(con, "raw", n = hlen))
  # shape
  sm <- regmatches(header, regexpr("'shape':\\s*\\(([^)]*)\\)", header))
  shape_str <- gsub("'shape':\\s*\\(|\\)", "", sm)
  shape <- if (nzchar(trimws(shape_str))) {
    as.integer(Filter(nzchar, trimws(strsplit(shape_str, ",")[[1]])))
  } else integer(0L)
  n_el <- max(1L, prod(shape))
  data <- readBin(con, "numeric", n = n_el, size = 4L, endian = "little")
  list(data = data, shape = shape)
}

# Unzip .npz and return a named list of torch tensors
.load_npz_state_dict <- function(npz_path) {
  tmp <- tempfile()
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  utils::unzip(npz_path, exdir = tmp)
  files <- list.files(tmp, pattern = "\\.npy$", full.names = TRUE)
  result <- lapply(files, function(f) {
    arr <- .read_npy(f)
    t <- torch::torch_tensor(arr$data, dtype = torch::torch_float32())
    if (length(arr$shape) > 1L) t <- t$reshape(arr$shape)
    t
  })
  names(result) <- tools::file_path_sans_ext(basename(files))
  result
}

# ---------------------------------------------------------------------------
# Model architectures
# ---------------------------------------------------------------------------

.DeepFormantsEstimator <- torch::nn_module(
  "DeepFormantsEstimator",
  initialize = function() {
    self$Dense1 <- torch::nn_linear(350L, 1024L)
    self$Dense2 <- torch::nn_linear(1024L, 512L)
    self$Dense3 <- torch::nn_linear(512L, 256L)
    self$out    <- torch::nn_linear(256L, 4L)
  },
  forward = function(x) {
    x <- torch::torch_sigmoid(self$Dense1(x))
    x <- torch::torch_sigmoid(self$Dense2(x))
    x <- torch::torch_sigmoid(self$Dense3(x))
    self$out(x)
  }
)

.DeepFormantsTracker <- torch::nn_module(
  "DeepFormantsTracker",
  initialize = function() {
    self$lstm1 <- torch::nn_lstm(350L, 512L, batch_first = TRUE)
    self$lstm2 <- torch::nn_lstm(512L, 256L, batch_first = TRUE)
    self$fc    <- torch::nn_linear(256L, 4L)
  },
  forward = function(x) {
    x <- self$lstm1(x)[[1]]
    x <- self$lstm2(x)[[1]]
    self$fc(x)
  }
)

# ---------------------------------------------------------------------------
# Model loaders (cached per session)
# ---------------------------------------------------------------------------

.df_model_cache <- new.env(parent = emptyenv())

#' Load bundled DeepFormants Estimator (MLP) weights
#' @keywords internal
load_deepformants_estimator <- function() {
  if (!is.null(.df_model_cache$estimator)) return(.df_model_cache$estimator)
  npz <- system.file("python", "DeepFormants", "pytorchFormants",
                     "Estimator", "LPC_NN_scaledLoss.npz",
                     package = "superassp")
  if (!nzchar(npz)) stop("Estimator weights not found in package.", call. = FALSE)
  model <- .DeepFormantsEstimator()
  model$load_state_dict(.load_npz_state_dict(npz))
  model$eval()
  .df_model_cache$estimator <- model
  model
}

#' Load bundled DeepFormants Tracker (LSTM) weights
#' @keywords internal
load_deepformants_tracker <- function() {
  if (!is.null(.df_model_cache$tracker)) return(.df_model_cache$tracker)
  npz <- system.file("python", "DeepFormants", "pytorchFormants",
                     "Tracker", "LPC_RNN.npz",
                     package = "superassp")
  if (!nzchar(npz)) stop("Tracker weights not found in package.", call. = FALSE)
  model <- .DeepFormantsTracker()
  # Python saves LSTM weights with 0-indexed suffix (_l0); R torch expects _l1
  state <- .load_npz_state_dict(npz)
  names(state) <- gsub("_l0$", "_l1", names(state))
  model$load_state_dict(state)
  model$eval()
  .df_model_cache$tracker <- model
  model
}

# ---------------------------------------------------------------------------
# Inference helpers
# ---------------------------------------------------------------------------

#' Run DeepFormants LSTM tracker on feature matrix
#' @param feat_mat (n_frames, 350) numeric matrix
#' @return (n_frames, 4) matrix of formant frequencies in Hz
#' @keywords internal
run_deepformants_tracker <- function(feat_mat) {
  model <- load_deepformants_tracker()
  x <- torch::torch_tensor(feat_mat, dtype = torch::torch_float32())$unsqueeze(1L)
  torch::with_no_grad(out <- model(x))
  as.matrix(out$squeeze(1L)$cpu()) * 1000.0
}

#' Run DeepFormants MLP estimator on single feature vector
#' @param feat_mat (1, 350) numeric matrix
#' @return length-4 numeric vector of formant frequencies in Hz
#' @keywords internal
run_deepformants_estimator <- function(feat_mat) {
  model <- load_deepformants_estimator()
  x <- torch::torch_tensor(feat_mat, dtype = torch::torch_float32())
  torch::with_no_grad(out <- model(x))
  as.numeric(out$cpu()) * 1000.0
}

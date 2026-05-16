
#' Predict formants using the FormantNet neural formant tracker (ONNX)
#'
#' FormantNet \insertCite{Sakamoto.2021}{superassp} is an LSTM-based formant
#' tracker that operates on log-scale smoothed spectral envelopes extracted at
#' 5 ms intervals. Inference uses ONNX Runtime via the same C++ backend as
#' \code{\link{trk_pitch_crepe}} — no Python or TensorFlow required.
#'
#' Pre-processing pipeline (fixed by model training):
#' \enumerate{
#'   \item Resample to 16 kHz.
#'   \item Pre-emphasis (factor 0.98).
#'   \item 512-sample (32 ms) Hann-windowed STFT, 80-sample (5 ms) hop.
#'   \item Spectral envelope smoothing — 6 passes of two-stage binomial
#'         smoothing on the linear-scale magnitude spectrum.
#'   \item Convert to dB: \eqn{20 \log_{10}(c \cdot |X| + 0.001)}, where
#'         \eqn{c = 2 / \sum w} and \eqn{w} is the Hann window.
#'   \item Retain the first 257 frequency bins (0–8 kHz).
#'   \item Normalise by training-set global mean and standard deviation.
#' }
#'
#' Post-processing pipeline (matching the original FormantNet scripts):
#' \enumerate{
#'   \item Rescale sigmoid output to Hz / dB.
#'   \item Sort formants by mean frequency (ascending).
#'   \item Apply binomial smoothing over time (10 passes).
#' }
#'
#' @param listOfFiles Character vector of audio file paths (any format
#'   supported by av).
#' @param beginTime Start time (seconds). Default 0 (beginning of file).
#' @param endTime End time (seconds). Default 0 (end of file).
#' @param numFormants Integer (1–6). Number of formants to return.
#'   The model always predicts 6 formants internally; this selects the lowest
#'   \code{numFormants} after sorting by mean frequency. Default 3.
#' @param toFile Logical. If \code{TRUE}, write SSFF files; if \code{FALSE},
#'   return an \code{AsspDataObj}. Default \code{TRUE}.
#' @param explicitExt Output file extension. Default \code{"fnf"}.
#' @param outputDirectory Output directory. Default \code{NULL} (same as
#'   input file).
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks
#'   \code{fm} (REAL32, Hz, \emph{n\_frames} × \code{numFormants}) and
#'   \code{bw} (REAL32, Hz, \emph{n\_frames} × \code{numFormants}).
#'   Frame rate is 200 Hz (5 ms hop).
#'   If \code{toFile = TRUE}: the number of files written (invisibly).
#'
#' @details
#' ONNX Runtime is installed automatically on first use (~30 MB, cached in
#' the R user directory). The model file (\code{inst/onnx/formantnet/formantnet.onnx},
#' 6.4 MB) is bundled with the package.
#'
#' Track access: \code{result$fm[, 1]} = F1, \code{result$fm[, 2]} = F2, etc.
#'
#' @export
#'
#' @references
#'   \insertAllCited{}
#'
trk_formant_formantnet <- function(listOfFiles,
                                   beginTime       = 0.0,
                                   endTime         = 0.0,
                                   numFormants     = 3L,
                                   toFile          = TRUE,
                                   explicitExt     = "fnf",
                                   outputDirectory = NULL,
                                   verbose         = TRUE) {

  # ── Guards ──────────────────────────────────────────────────────────────────
  ensure_onnx()

  numFormants <- as.integer(numFormants)
  if (numFormants < 1L || numFormants > 6L) {
    cli::cli_abort("{.arg numFormants} must be between 1 and 6, got {numFormants}.")
  }

  model_path <- system.file("onnx", "formantnet", "formantnet.onnx",
                            package = "superassp")
  if (!nzchar(model_path) || !file.exists(model_path)) {
    cli::cli_abort(c(
      "FormantNet ONNX model not found.",
      "i" = "Expected at inst/onnx/formantnet/formantnet.onnx"
    ))
  }

  norm_path <- system.file("onnx", "formantnet", "normstats.txt",
                           package = "superassp")
  if (!nzchar(norm_path) || !file.exists(norm_path)) {
    cli::cli_abort(c(
      "FormantNet normstats file not found.",
      "i" = "Expected at inst/onnx/formantnet/normstats.txt"
    ))
  }

  norm_lines <- readLines(norm_path, warn = FALSE)
  norm_mean  <- as.numeric(norm_lines[1L])
  norm_std   <- as.numeric(norm_lines[2L])

  if (length(listOfFiles) > 1L && !toFile) {
    cli::cli_abort("{.arg toFile = FALSE} only permitted for a single file.")
  }

  n_files   <- length(listOfFiles)
  beginTime <- fast_recycle_times(beginTime, n_files)
  endTime   <- fast_recycle_times(endTime,   n_files)

  missing_files <- !file.exists(listOfFiles)
  if (any(missing_files)) {
    cli::cli_abort("File(s) not found: {.file {listOfFiles[missing_files]}}")
  }

  # ── ORT session (reused across files) ───────────────────────────────────────
  session <- ort_session(model_path)

  # ── Per-file loop ────────────────────────────────────────────────────────────
  outListOfFiles <- character(0L)
  outDataObj     <- NULL

  for (i in seq_along(listOfFiles)) {
    origSoundFile <- normalizePath(listOfFiles[[i]], mustWork = TRUE)
    bt <- beginTime[i]
    et <- endTime[i]

    if (verbose) {
      cli::cli_inform("FormantNet: {basename(origSoundFile)}")
    }

    tryCatch({
      # ── Load audio at 16 kHz ─────────────────────────────────────────────
      invisible(utils::capture.output(
        audio_data <- av::read_audio_bin(
          audio       = origSoundFile,
          start_time  = if (bt > 0) bt else NULL,
          end_time    = if (et > 0) et else NULL,
          channels    = 1L,
          sample_rate = 16000L
        ),
        type = "message"
      ))
      # av returns INT32 PCM; normalise to float matching FormantNet training
      # (training: tf.audio.decode_wav → [-1,1], then ×32768 for int16 scale)
      audio_int16 <- as.numeric(audio_data) / 65536.0 * 32768.0

      # ── Preprocessing ────────────────────────────────────────────────────
      spectra <- .formantnet_preprocess(audio_int16)
      spectra <- (spectra - norm_mean) / norm_std          # normalise

      # ── ONNX inference ───────────────────────────────────────────────────
      raw_params <- .formantnet_infer(spectra, session)    # n_frames × 20

      # ── Postprocessing ───────────────────────────────────────────────────
      params   <- .formantnet_postprocess(raw_params, numFormants)
      n_frames <- nrow(params$fm)

      # ── Build AsspDataObj ─────────────────────────────────────────────────
      sample_rate    <- 200.0   # 1000 / 5 ms hop
      start_time_ssff <- 0.0

      outDataObj <- list()
      attr(outDataObj, "trackFormats") <- c("REAL32", "REAL32")
      attr(outDataObj, "sampleRate")   <- sample_rate
      attr(outDataObj, "origFreq")     <- 16000.0
      attr(outDataObj, "startTime")    <- start_time_ssff
      attr(outDataObj, "startRecord")  <- 1L
      attr(outDataObj, "endRecord")    <- as.integer(n_frames)
      class(outDataObj) <- "AsspDataObj"
      AsspFileFormat(outDataObj) <- "SSFF"
      AsspDataFormat(outDataObj) <- 2L

      outDataObj <- addTrack(outDataObj, "fm", params$fm, "REAL32")
      outDataObj <- addTrack(outDataObj, "bw", params$bw, "REAL32")

      # ── Output ───────────────────────────────────────────────────────────
      base_name <- tools::file_path_sans_ext(basename(origSoundFile))
      out_dir   <- if (is.null(outputDirectory)) dirname(origSoundFile) else outputDirectory
      ssff_file <- file.path(out_dir, paste0(base_name, ".", explicitExt))
      attr(outDataObj, "filePath") <- as.character(ssff_file)

      if (toFile) {
        write.AsspDataObj(dobj = outDataObj, file = ssff_file)
        outListOfFiles <- c(outListOfFiles, ssff_file)
      }

    }, error = function(e) {
      cli::cli_warn(
        "FormantNet failed for {.file {basename(origSoundFile)}}: {conditionMessage(e)}"
      )
    })
  }

  if (toFile) invisible(length(outListOfFiles)) else outDataObj
}

# ── Function attributes ──────────────────────────────────────────────────────
attr(trk_formant_formantnet, "ext")             <- "fnf"
attr(trk_formant_formantnet, "tracks")          <- c("fm", "bw")
attr(trk_formant_formantnet, "outputType")      <- "SSFF"
attr(trk_formant_formantnet, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")


# ── Internal: STFT + spectral envelope preprocessing ────────────────────────
#
# Converts a float audio vector (int16 scale, 16 kHz) to a normalised
# log-scale smoothed spectral envelope matrix (n_frames × 257).
# All constants match the FormantNet training configuration:
#   PREEMPH = 0.98, WIN = 512, HOP = 80, SMOOTH_LINEAR = TRUE,
#   ENV_SMOOTH_PASSES = 6, FLOOR = 0.001, SPECTRUM_NPOINTS = 257.
#
.formantnet_preprocess <- function(audio) {
  win_len  <- 512L
  hop      <- 80L
  n_bins   <- 257L
  preemph  <- 0.98

  n <- length(audio)
  if (n < win_len) {
    cli::cli_abort("Audio too short for FormantNet (minimum {win_len} samples at 16 kHz).")
  }

  # 1. Pre-emphasis: y[k] = x[k] - 0.98 * x[k-1]
  audio_pe <- c(audio[1L], audio[-1L] - preemph * audio[-n])

  # 2. Symmetric padding: 256 zeros left, 79 zeros right
  audio_padded <- c(rep(0.0, win_len %/% 2L), audio_pe, rep(0.0, hop - 1L))

  # 3. Frame count matches TF tf.signal.stft with pad_end=FALSE
  n_frames <- 1L + (n - 177L) %/% hop
  if (n_frames < 1L) n_frames <- 1L

  # 4. Build frame matrix (vectorised: no explicit loop over frames)
  hann       <- 0.5 * (1.0 - cos(2.0 * pi * seq(0L, win_len - 1L) / (win_len - 1L)))
  offsets    <- (seq_len(n_frames) - 1L) * hop          # frame start indices (0-based)
  row_idx    <- seq(0L, win_len - 1L)                   # within-frame offsets
  # idx[j, i] = start of frame i + row offset j + 1 (1-based for R)
  idx        <- outer(row_idx, offsets, "+") + 1L        # win_len × n_frames
  frames     <- matrix(audio_padded[idx], nrow = win_len) * hann  # win_len × n_frames

  # 5. Batch FFT (mvfft acts on columns) → magnitude, first n_bins bins
  fft_mag    <- Mod(mvfft(frames)[seq_len(n_bins), ])    # n_bins × n_frames
  magnitude  <- t(fft_mag)                               # n_frames × n_bins

  # 6. Spectral envelope smoothing (linear scale, 6 passes, two-stage per pass)
  magnitude  <- .formantnet_smooth_spectral(magnitude, passes = 6L)

  # 7. dB conversion  c = 2 / sum(hann)
  scaling    <- 2.0 / sum(hann)
  20.0 * log10(scaling * magnitude + 0.001)              # n_frames × n_bins
}


# ── Internal: two-stage spectral envelope smoother ──────────────────────────
#
# Vectorised translation of FormantNet smooth_spenvl():
#   Pass 1 per pass: strict local minimum → replace with 0.5*(left + right)
#   Pass 2 per pass: any point ≤ either neighbour → weighted 0.25/0.5/0.25
# Both passes run sequentially within each of `passes` outer iterations.
# Smoothing is along the FREQUENCY axis (within each frame row).
#
.formantnet_smooth_spectral <- function(spec, passes = 6L) {
  n_bins <- ncol(spec)
  inner  <- seq(2L, n_bins - 1L)   # interior column indices

  for (p in seq_len(passes)) {
    # ── Stage 1: strict local minimum → average of neighbours ────────────
    left_s1  <- spec[, inner - 1L]
    right_s1 <- spec[, inner + 1L]
    avg_s1   <- 0.5 * (left_s1 + right_s1)
    is_min   <- spec[, inner] < left_s1 & spec[, inner] < right_s1
    spec[, inner][is_min] <- avg_s1[is_min]

    # ── Stage 2: at-or-below either neighbour → binomial blend ───────────
    left_s2  <- spec[, inner - 1L]
    right_s2 <- spec[, inner + 1L]
    blend_s2 <- 0.25 * left_s2 + 0.5 * spec[, inner] + 0.25 * right_s2
    at_or_below <- spec[, inner] <= left_s2 | spec[, inner] <= right_s2
    spec[, inner][at_or_below] <- blend_s2[at_or_below]
  }
  spec
}


# ── Internal: ONNX inference ─────────────────────────────────────────────────
#
# Sends a normalised spectrogram matrix (n_frames × 257) to the FormantNet
# ONNX model and returns the raw sigmoid output (n_frames × 20).
# ORT requires a row-major flat vector; R matrices are column-major, so
# we transpose before flattening.
#
.formantnet_infer <- function(spectra, session) {
  n_frames   <- nrow(spectra)
  input_flat <- c(t(spectra))          # row-major: features vary fastest
  result     <- ort_run(
    session,
    inputs       = list(input = input_flat),
    shapes       = list(c(1L, as.integer(n_frames), 257L)),
    output_names = NULL
  )
  # result[[1]] is a flat vector with shape attribute (1, n_frames, 20)
  matrix(result[[1L]], nrow = n_frames, ncol = 20L, byrow = TRUE)
}


# ── Internal: postprocessing ─────────────────────────────────────────────────
#
# Converts raw sigmoid output (n_frames × 20) to formant frequencies and
# bandwidths in Hz, following FN_model.get_rescale_fn() and track_files():
#   Cols  1– 7: frequencies F1…F6, Fz  (raw × 8000 Hz)
#   Cols  8–14: bandwidths  B1…B6, Bz  (20 + raw × 4980 Hz)
#   Cols 15–20: amplitudes  A1…A6      (-100 + raw × 200 dB, unused here)
# After rescaling: sort F1–F6 by mean frequency, then 10-pass time smoothing.
#
.formantnet_postprocess <- function(raw_params, numFormants = 3L) {
  # Rescale
  freq <- raw_params[, 1L:7L,  drop = FALSE] * 8000.0
  bw   <- 20.0 + raw_params[, 8L:14L, drop = FALSE] * 4980.0

  # Negate antiformant frequency (col 7 = Fz → negative by convention)
  freq[, 7L] <- -abs(freq[, 7L])

  # Sort formants 1–6 by mean frequency (ascending)
  ord        <- order(colMeans(freq[, 1L:6L, drop = FALSE]))
  freq[, 1L:6L] <- freq[, ord, drop = FALSE]
  bw[,   1L:6L] <- bw[,   ord, drop = FALSE]

  # Binomial smoothing over the TIME axis (10 passes)
  freq[, 1L:6L] <- .formantnet_smooth_time(freq[, 1L:6L, drop = FALSE], passes = 10L)
  bw[,   1L:6L] <- .formantnet_smooth_time(bw[,   1L:6L, drop = FALSE], passes = 10L)

  list(
    fm = freq[, seq_len(numFormants), drop = FALSE],
    bw = bw[,   seq_len(numFormants), drop = FALSE]
  )
}


# ── Internal: binomial time smoother ─────────────────────────────────────────
#
# Applies `passes` iterations of weighted-average smoothing along the TIME
# axis (rows) of a matrix. Edge frames use asymmetric weights (0.75/0.25).
#
.formantnet_smooth_time <- function(mat, passes = 10L) {
  n <- nrow(mat)
  if (n < 2L) return(mat)

  for (p in seq_len(passes)) {
    if (n >= 2L) {
      mat[1L, ]      <- 0.75 * mat[1L, ]      + 0.25 * mat[2L, ]
      mat[n,  ]      <- 0.75 * mat[n,  ]      + 0.25 * mat[n - 1L, ]
    }
    if (n >= 3L) {
      mat[2L:(n-1L), ] <- 0.5  * mat[2L:(n-1L), ] +
                          0.25 * mat[1L:(n-2L), ] +
                          0.25 * mat[3L:n,      ]
    }
  }
  mat
}

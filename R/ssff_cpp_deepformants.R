
#' Predict formants using the DeepFormants neural formant tracker (ONNX)
#'
#' DeepFormants \insertCite{Dissen.2016}{superassp} is an LSTM-based formant
#' tracker that operates on composite LPC/cepstral features extracted at 10 ms
#' intervals. Inference uses ONNX Runtime via the same C++ backend as
#' \code{\link{trk_pitch_crepe}} — no Python required.
#'
#' Pre-processing pipeline (fixed by model training):
#' \enumerate{
#'   \item Resample to 16 kHz, treat as raw int16 values (no normalisation).
#'   \item Extract 480-sample (30 ms) frames with 160-sample (10 ms) hop.
#'   \item For each frame, compute a 350-dimensional feature vector:
#'     \itemize{
#'       \item \strong{specPS} (50 coefficients): averaged periodogram across
#'             9 sub-frames → \eqn{\log\sqrt{P^2+f^2}} → \eqn{\log_{10}} →
#'             DCT-II ortho → first 50 coefficients.
#'       \item \strong{arspecs × 10} (30 × 10 = 300 coefficients): biased
#'             autocorrelation → Levinson-Durbin LPC (orders 8–17) → AR
#'             spectrum → same transform → first 30 coefficients each.
#'     }
#' }
#'
#' Post-processing: multiply raw model output by 1000 to obtain Hz.
#' No bandwidth output; model predicts F1–F4 only.
#'
#' @param listOfFiles Character vector of audio file paths (any format
#'   supported by av).
#' @param beginTime Start time (seconds). Default 0 (beginning of file).
#' @param endTime End time (seconds). Default 0 (end of file).
#' @param numFormants Integer (1–4). Number of formants to return.
#'   The model always predicts 4 formants; this selects the first
#'   \code{numFormants}. Default 3.
#' @param toFile Logical. If \code{TRUE}, write SSFF files; if \code{FALSE},
#'   return an \code{AsspDataObj}. Default \code{TRUE}.
#' @param explicitExt Output file extension. Default \code{"dff"}.
#' @param outputDirectory Output directory. Default \code{NULL} (same as
#'   input file).
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track
#'   \code{fm} (REAL32, Hz, \emph{n\_frames} × \code{numFormants}).
#'   Frame rate is 100 Hz (10 ms hop). No bandwidth track is produced.
#'   If \code{toFile = TRUE}: the number of files written (invisibly).
#'
#' @details
#' ONNX Runtime is installed automatically on first use (~30 MB, cached in
#' the R user directory). The model file
#' (\code{inst/onnx/deepformants/lpc_tracker.onnx}, ~10 MB) is bundled with
#' the package.
#'
#' Track access: \code{result$fm[, 1]} = F1, \code{result$fm[, 2]} = F2, etc.
#'
#' @export
#'
#' @references
#'   \insertAllCited{}
#'
trk_formant_deepformants <- function(listOfFiles,
                                     beginTime       = 0.0,
                                     endTime         = 0.0,
                                     numFormants     = 3L,
                                     toFile          = TRUE,
                                     explicitExt     = "dff",
                                     outputDirectory = NULL,
                                     verbose         = TRUE) {

  # ── Guards ──────────────────────────────────────────────────────────────────
  ensure_onnx()

  numFormants <- as.integer(numFormants)
  if (numFormants < 1L || numFormants > 4L) {
    cli::cli_abort("{.arg numFormants} must be between 1 and 4, got {numFormants}.")
  }

  model_path <- system.file("onnx", "deepformants", "lpc_tracker.onnx",
                            package = "superassp")
  if (!nzchar(model_path) || !file.exists(model_path)) {
    cli::cli_abort(c(
      "DeepFormants ONNX model not found.",
      "i" = "Expected at inst/onnx/deepformants/lpc_tracker.onnx"
    ))
  }

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
      cli::cli_inform("DeepFormants: {basename(origSoundFile)}")
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
      # av returns INT32 PCM; divide by 65536 to recover int16 scale
      # matching Python np.frombuffer(dstr, np.int16)
      audio_int16 <- as.numeric(audio_data) / 65536.0

      # ── Frame extraction ─────────────────────────────────────────────────
      N <- length(audio_int16)
      if (N < 480L) {
        cli::cli_warn(
          "DeepFormants: {.file {basename(origSoundFile)}} too short (< 30 ms), skipping."
        )
        next
      }
      frame_starts <- seq(0L, N - 480L, by = 160L)
      n_frames     <- length(frame_starts)
      idx          <- outer(seq(0L, 479L), frame_starts, "+") + 1L
      frames       <- matrix(audio_int16[idx], nrow = 480L)   # 480 × n_frames

      # ── Feature extraction ───────────────────────────────────────────────
      features <- .deepformants_build_features(frames)   # 350 × n_frames

      # ── ONNX inference ───────────────────────────────────────────────────
      raw_params <- .deepformants_infer(features, session)   # n_frames × 4

      # ── Postprocessing ───────────────────────────────────────────────────
      formant_hz <- raw_params * 1000.0
      fm         <- formant_hz[, seq_len(numFormants), drop = FALSE]

      # ── Build AsspDataObj ─────────────────────────────────────────────────
      outDataObj <- list()
      attr(outDataObj, "trackFormats") <- "REAL32"
      attr(outDataObj, "sampleRate")   <- 100.0    # 1000 / 10 ms
      attr(outDataObj, "origFreq")     <- 16000.0
      attr(outDataObj, "startTime")    <- 0.0
      attr(outDataObj, "startRecord")  <- 1L
      attr(outDataObj, "endRecord")    <- as.integer(n_frames)
      class(outDataObj) <- "AsspDataObj"
      AsspFileFormat(outDataObj) <- "SSFF"
      AsspDataFormat(outDataObj) <- 2L

      outDataObj <- addTrack(outDataObj, "fm", fm, "REAL32")

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
        "DeepFormants failed for {.file {basename(origSoundFile)}}: {conditionMessage(e)}"
      )
    })
  }

  if (toFile) invisible(length(outListOfFiles)) else outDataObj
}

# ── Function attributes ──────────────────────────────────────────────────────
attr(trk_formant_deepformants, "ext")             <- "dff"
attr(trk_formant_deepformants, "tracks")          <- "fm"
attr(trk_formant_deepformants, "outputType")      <- "SSFF"
attr(trk_formant_deepformants, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")


# ── Internal: build 350-dim feature matrix ───────────────────────────────────
#
# Applies build_single_feature_row() from DeepFormants extract_features.py
# to every column (frame) of the 480 × n_frames matrix.
# Returns a 350 × n_frames matrix (50 specPS + 10 × 30 arspecs).
#
.deepformants_build_features <- function(frames) {
  lpc_orders <- 8L:17L
  apply(frames, 2L, function(x) {
    spec <- .deepformants_specPS(x)
    ars  <- unlist(lapply(lpc_orders, function(ord) .deepformants_arspecs(x, ord)))
    feat <- c(spec, ars)
    feat[!is.finite(feat)] <- 0.0
    feat
  })
}


# ── Internal: specPS ─────────────────────────────────────────────────────────
#
# Replicates specPS(data, pitch=50) from extract_features.py.
# Splits a 480-sample frame into 9 sub-frames of 50 samples,
# computes 4096-pt periodograms, averages them (dividing by 9.6 not 9),
# applies log(sqrt(pxx^2 + fgrid^2)), log10, DCT-II ortho, returns first 50.
#
.deepformants_specPS <- function(x) {
  samps    <- 480.0 / 50.0           # = 9.6 (float — used as divisor)
  n_subfr  <- as.integer(samps)      # = 9  (integer loop count)
  samp_len <- as.integer(480.0 / samps)  # = 50 samples per sub-frame
  pxx_avg  <- numeric(2049L)

  for (i in seq(0L, n_subfr - 1L)) {
    seg     <- x[(samp_len * i + 1L):(samp_len * (i + 1L))]
    sp      <- Mod(fft(c(seg, rep(0.0, 4096L - samp_len))))[1L:2049L]^2 /
                 (samp_len * 1.0)
    pxx_avg <- pxx_avg + sp
  }
  pxx_avg <- pxx_avg / samps          # divide by 9.6

  fgrid <- seq(0.0, 0.5, length.out = 2049L)
  peri  <- log(sqrt(pxx_avg^2 + fgrid^2))
  peri[!is.finite(peri) | peri == 0.0] <- log(1e-10)
  # log10 of negative peri values (log of magnitudes < 1) produces NaN;
  # this matches Python: np.log10(peri) → NaN → replaced with 0.0 upstream
  mspec <- suppressWarnings(log10(peri))
  .deepformants_dct(mspec)[1L:50L]
}


# ── Internal: arspecs ────────────────────────────────────────────────────────
#
# Replicates arspecs(data, order, Atal=False) from extract_features.py.
# Biased autocorrelation → Levinson-Durbin LPC → 4096-pt AR spectrum →
# log(sqrt(pxx^2 + fgrid^2)) → log10 → DCT-II ortho → first 30 coefficients.
#
.deepformants_arspecs <- function(x, order) {
  n    <- length(x)
  # Biased autocorrelation via FFT (matches acorr_lpc from levinson_lpc.py)
  nfft <- 2L^as.integer(ceiling(log2(2.0 * n - 1.0)))
  raw  <- Re(fft(Mod(fft(c(x, rep(0.0, nfft - n))))^2, inverse = TRUE))
  r    <- raw[1L:(n + 1L)] / n      # biased: divide by n

  ld   <- .deepformants_levinson(r, order)
  a    <- ld$a
  e    <- ld$e

  # AR spectrum: px = 1 / FFT(a, 4096), pxx = |px|^2 * e (since fs=1)
  pn  <- 2049L
  px  <- 1.0 / fft(c(a, rep(0.0, 4096L - length(a))))[1L:pn]
  pxx <- Re(Conj(px) * px) * e

  fgrid <- seq(0.0, 0.5, length.out = pn)
  ar    <- log(sqrt(pxx^2 + fgrid^2))
  ar[!is.finite(ar) | ar == 0.0] <- log(1e-10)
  mspec <- suppressWarnings(log10(ar))
  .deepformants_dct(mspec)[1L:30L]
}


# ── Internal: Levinson-Durbin recursion ──────────────────────────────────────
#
# Exact port of levinson_1d() from levinson_lpc.py (DeepFormants repo).
# Returns list(a, e): LPC coefficients (length order+1, a[1]=1) and
# residual prediction error.
#
.deepformants_levinson <- function(r, order) {
  a    <- numeric(order + 1L); a[1L] <- 1.0
  t    <- numeric(order + 1L)
  k    <- numeric(order)
  e    <- r[1L]
  for (i in seq_len(order)) {
    acc <- r[i + 1L]
    if (i > 1L) {
      for (j in seq_len(i - 1L)) acc <- acc + a[j + 1L] * r[i - j + 1L]
    }
    k[i]     <- -acc / e
    a[i + 1L] <- k[i]
    t[]       <- a
    if (i > 1L) {
      for (j in seq_len(i - 1L)) a[j + 1L] <- a[j + 1L] + k[i] * t[i - j + 1L]
    }
    e <- e * (1.0 - k[i]^2)
  }
  list(a = a, e = e)
}


# ── Internal: DCT-II ortho ───────────────────────────────────────────────────
#
# Replicates scipy.fftpack.dct(x, type=2, norm='ortho').
# Uses even-symmetry FFT method: reorder → FFT → phase-shift → normalise.
# y[0] /= sqrt(4n),  y[k>0] /= sqrt(2n).
#
.deepformants_dct <- function(x) {
  n <- length(x)
  w <- c(x[seq(1L, n, 2L)], x[seq(n, 2L, -2L)])
  y <- Re(fft(w) * exp(-1i * pi * seq(0L, n - 1L) / (2.0 * n)))
  y[1L]  <- y[1L]  / sqrt(4.0 * n)
  y[-1L] <- y[-1L] / sqrt(2.0 * n)
  y
}


# ── Internal: ONNX inference ─────────────────────────────────────────────────
#
# Sends a 350 × n_frames feature matrix to the DeepFormants LSTM tracker
# and returns raw output (n_frames × 4). ORT requires row-major input;
# R matrices are column-major, so transpose before flattening.
#
.deepformants_infer <- function(features, session) {
  n_frames   <- ncol(features)
  input_flat <- c(t(features))     # row-major: time varies slowest
  result     <- ort_run(
    session,
    inputs       = list(input = input_flat),
    shapes       = list(c(1L, as.integer(n_frames), 350L)),
    output_names = NULL
  )
  # result[[1]] is flat (1 × n_frames × 4) → n_frames × 4
  matrix(result[[1L]], nrow = n_frames, ncol = 4L, byrow = TRUE)
}

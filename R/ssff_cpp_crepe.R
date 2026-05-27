
#' Track fundamental frequency and periodicity using CREPE (ONNX)
#'
#' Estimates F0 by applying a deep convolutional neural network directly to the
#' time-domain waveform
#' (CREPE; \insertCite{Kim.2018.10.1109/icassp.2018.8461329}{superassp}). Two
#' model sizes are available: \code{"tiny"} (~1.9 MB, fast) and \code{"full"}
#' (~85 MB, more accurate). No Python or PyTorch required — inference uses the
#' bundled ONNX Runtime.
#'
#' @param listOfFiles Character vector of audio file paths. Any format
#'   supported by \pkg{av} is accepted.
#' @param beginTime Numeric. Start of analysis window in seconds. Default 0
#'   (file start).
#' @param endTime Numeric. End of analysis window in seconds. Default 0
#'   (file end).
#' @param windowShift Numeric. Frame shift in milliseconds; sets output frame
#'   rate (\code{1000 / windowShift} Hz). Default 10 ms (100 Hz).
#' @param windowSize Numeric. Smoothing filter window size in milliseconds,
#'   applied to both median (periodicity) and mean (F0) post-processing filters.
#'   Default 15 ms.
#' @param minF Numeric. Minimum F0 in Hz. Frames with estimated F0 below this
#'   value are set to 0 (unvoiced). Default 50 Hz.
#' @param maxF Numeric. Maximum F0 in Hz. Frames above this value are set to 0.
#'   Default 550 Hz.
#' @param voicing.threshold Numeric. Periodicity threshold for voicing
#'   decisions. Frames with periodicity below this value are set to unvoiced
#'   (F0 = 0). Default 0.21.
#' @param silence.threshold Numeric. A-weighted dB threshold relative to global
#'   maximum. Frames below this are treated as silent (periodicity = 0).
#'   Default -60 dB.
#' @param model Character. Model variant: \code{"tiny"} (fast) or
#'   \code{"full"} (more accurate). Default \code{"tiny"}.
#' @param decoder Character. Decoding method: \code{"viterbi"} (recommended,
#'   smoother) or \code{"argmax"} (faster but noisier). Default
#'   \code{"viterbi"}.
#' @param batch_size Integer. Number of frames per ONNX inference batch.
#'   Default 512.
#' @param explicitExt Character. Output file extension. Default \code{"crp"}.
#' @param outputDirectory Character. Directory for output files. \code{NULL}
#'   (default) writes alongside the input file.
#' @param toFile Logical. If \code{TRUE}, write SSFF output files and return
#'   the count written (invisibly). If \code{FALSE}, return an
#'   \code{AsspDataObj}. Default \code{TRUE}.
#' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
#'
#' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
#'   \describe{
#'     \item{\code{f0}}{REAL32, Hz, \emph{n\_frames} × 1. Fundamental
#'       frequency; 0 in unvoiced/silent frames.}
#'     \item{\code{periodicity}}{REAL32, 0–1, \emph{n\_frames} × 1. Model
#'       confidence; values below \code{voicing.threshold} are treated as
#'       unvoiced.}
#'   }
#'   Frame rate: \code{1000 / windowShift} Hz (default 100 Hz, 10 ms hop).
#'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
#'
#' @details
#' ONNX Runtime is installed automatically on first use (~30 MB, cached in
#' the R user directory) and persists across R sessions and package reinstalls.
#'
#' Post-processing (matching torchcrepe): median filter on periodicity →
#' A-weighted silence detection → voicing threshold → NaN-aware mean filter
#' on F0 → F0 range clamping.
#'
#' @examples
#' \dontrun{
#' trk_pitch_crepe(
#'   system.file("samples", "sustained", "a1.wav", package = "superassp"),
#'   toFile = FALSE
#' )
#' }
#' @export
#'
#' @references
#'   \insertAllCited{}
#'
trk_pitch_crepe <- function(listOfFiles,
                      beginTime = 0,
                      endTime = 0,
                      windowShift = 10,
                      windowSize = 15,
                      minF = 50,
                      maxF = 550,
                      voicing.threshold = 0.21,
                      silence.threshold = -60.0,
                      model = c("tiny", "full"),
                      decoder = c("viterbi", "argmax"),
                      batch_size = 512L,
                      explicitExt = "crp",
                      outputDirectory = NULL,
                      toFile = TRUE,
                      verbose = TRUE) {

  model <- match.arg(model)
  decoder <- match.arg(decoder)

  # Ensure ONNX Runtime is available (auto-installs if needed)
  ensure_onnx()

  # Locate ONNX model
  model_path <- system.file("onnx", "crepe", paste0(model, ".onnx"),
                            package = "superassp")
  if (!nzchar(model_path) || !file.exists(model_path)) {
    stop("CREPE ", model, " ONNX model not found. Expected at inst/onnx/crepe/",
         model, ".onnx", call. = FALSE)
  }

  # Input validation
  if (length(listOfFiles) > 1 && !toFile) {
    stop("toFile=FALSE only permitted for single files.", call. = FALSE)
  }

  n_files <- length(listOfFiles)
  beginTime <- fast_recycle_times(beginTime, n_files)
  endTime <- fast_recycle_times(endTime, n_files)

  filesEx <- file.exists(listOfFiles)
  if (!all(filesEx)) {
    stop("File(s) not found: ",
         paste(listOfFiles[!filesEx], collapse = ", "), call. = FALSE)
  }

  outListOfFiles <- c()
  outDataObj <- NULL
  target_sr <- 16000L  # CREPE expects 16 kHz

  for (i in seq_along(listOfFiles)) {
    origSoundFile <- normalizePath(listOfFiles[i], mustWork = TRUE)
    bt <- beginTime[i]
    et <- endTime[i]

    if (verbose) {
      cli::cli_inform("Processing {basename(origSoundFile)} ({model} model, {decoder})")
    }

    # ── Load audio ──────────────────────────────────────────────────────────
    invisible(utils::capture.output(
      audio_data <- av::read_audio_bin(
        audio = origSoundFile,
        start_time = if (bt > 0) bt else NULL,
        end_time = if (et > 0) et else NULL,
        channels = 1,
        sample_rate = target_sr
      ),
      type = "message"
    ))
    sr <- attr(audio_data, "sample_rate")
    audio_float <- as.numeric(audio_data) / 2147483647.0  # INT32_MAX

    # ── C++ ONNX inference ──────────────────────────────────────────────────
    hop_length <- as.integer(sr * windowShift / 1000.0)

    raw <- crepe_inference_cpp(
      audio = audio_float,
      sample_rate = as.double(sr),
      model_path = model_path,
      hop_length = hop_length,
      batch_size = as.integer(batch_size),
      use_viterbi = (decoder == "viterbi")
    )

    n_frames <- raw$n_frames
    f0 <- raw$f0
    periodicity <- raw$periodicity

    # ── Post-processing (matching torchcrepe pipeline) ──────────────────────

    # 1. Median filter on periodicity
    filter_win <- max(3L, as.integer(ceiling(windowSize / windowShift)))
    if (filter_win %% 2 == 0) filter_win <- filter_win + 1L  # must be odd
    periodicity <- stats::runmed(periodicity, k = filter_win, endrule = "constant")

    # 2. A-weighted silence detection
    periodicity <- .crepe_silence_mask(
      periodicity, audio_float, sr, hop_length,
      silence.threshold, n_frames
    )

    # 3. Voicing threshold: set F0 to 0 where periodicity < threshold
    unvoiced <- periodicity < voicing.threshold
    f0[unvoiced] <- 0.0

    # 4. NaN-aware mean filter on F0
    # Replace 0 with NA, apply mean filter, restore 0
    f0_na <- ifelse(f0 == 0, NA_real_, f0)
    f0_na <- .nanmean_filter(f0_na, filter_win)
    f0 <- ifelse(is.na(f0_na), 0.0, f0_na)

    # 5. F0 range clamping
    out_of_range <- f0 > 0 & (f0 < minF | f0 > maxF)
    f0[out_of_range] <- 0.0
    periodicity[out_of_range] <- 0.0

    # ── Build AsspDataObj ───────────────────────────────────────────────────
    sampleRate <- 1000.0 / windowShift
    startTime_ssff <- 1.0 / sampleRate

    outDataObj <- list()
    attr(outDataObj, "trackFormats") <- c("REAL32", "REAL32")
    attr(outDataObj, "sampleRate") <- sampleRate
    attr(outDataObj, "origFreq") <- as.numeric(sr)
    attr(outDataObj, "startTime") <- startTime_ssff
    attr(outDataObj, "startRecord") <- 1L
    attr(outDataObj, "endRecord") <- as.integer(n_frames)
    class(outDataObj) <- "AsspDataObj"

    AsspFileFormat(outDataObj) <- "SSFF"
    AsspDataFormat(outDataObj) <- as.integer(2)  # binary

    outDataObj <- addTrack(outDataObj, "f0",
                           matrix(f0, ncol = 1), "REAL32")
    outDataObj <- addTrack(outDataObj, "periodicity",
                           matrix(periodicity, ncol = 1), "REAL32")

    # ── File output ─────────────────────────────────────────────────────────
    base_name <- tools::file_path_sans_ext(basename(origSoundFile))
    out_dir <- if (is.null(outputDirectory)) dirname(origSoundFile) else outputDirectory
    ssff_file <- file.path(out_dir, paste0(base_name, ".", explicitExt))

    attr(outDataObj, "filePath") <- as.character(ssff_file)

    if (toFile) {
      write.AsspDataObj(dobj = outDataObj, file = ssff_file)
      outListOfFiles <- c(outListOfFiles, ssff_file)
    }
  }

  if (toFile) {
    return(invisible(length(outListOfFiles)))
  } else {
    return(outDataObj)
  }
}

# ── Function attributes ──────────────────────────────────────────────────────
attr(trk_pitch_crepe, "ext") <- "crp"
attr(trk_pitch_crepe, "tracks") <- c("f0", "periodicity")
attr(trk_pitch_crepe, "outputType") <- "SSFF"
attr(trk_pitch_crepe, "nativeFiletypes") <- c("wav")


# ── Internal: A-weighted silence detection ───────────────────────────────────
#
# Computes A-weighted power per frame using pladdrr spectrogram, then zeros
# periodicity for frames whose dB level is below threshold relative to max.
#
# Falls back to simple peak-amplitude method if pladdrr is not available.
#
.crepe_silence_mask <- function(periodicity, audio_float, sr, hop_length,
                                silence_db, n_frames) {
  # Try pladdrr spectrogram approach
  if (requireNamespace("pladdrr", quietly = TRUE)) {
    tryCatch({
      return(.crepe_silence_pladdrr(periodicity, audio_float, sr,
                                     hop_length, silence_db, n_frames))
    }, error = function(e) {
      # Fall back to simple method
    })
  }

  # Fallback: simple per-frame RMS in dB
  .crepe_silence_simple(periodicity, audio_float, sr, hop_length,
                        silence_db, n_frames)
}


# A-weighted silence via pladdrr spectrogram
.crepe_silence_pladdrr <- function(periodicity, audio_float, sr, hop_length,
                                    silence_db, n_frames) {
  # Write audio to temp WAV for pladdrr Sound
  tmp <- tempfile(fileext = ".wav")
  .write_minimal_wav(audio_float, sr, tmp)

  snd <- pladdrr::Sound$new(path = tmp)
  on.exit({
    rm(snd)      # Release pladdrr R6 object first (drops C++ xptr ref)
    gc(FALSE)    # Run finalizers before deleting backing file
    unlink(tmp)
  }, add = TRUE)
  time_step <- hop_length / sr

  spec <- snd$to_spectrogram(
    window_length = 0.064,  # 64ms Gaussian window
    time_step = time_step
  )
  mat <- spec$as_matrix()  # [n_freq x n_time] power spectral density

  freqs <- as.numeric(rownames(mat))

  # IEC 61672 A-weighting transfer function
  # A(f) = 12194^2 * f^4 / ((f^2+20.6^2) * sqrt((f^2+107.7^2)*(f^2+737.9^2)) * (f^2+12194^2))
  # Weight in dB = 20*log10(A(f)) + 2.0 (normalization constant)
  f2 <- freqs^2
  a_num <- 12194^2 * f2^2
  a_den <- (f2 + 20.6^2) * sqrt((f2 + 107.7^2) * (f2 + 737.9^2)) * (f2 + 12194^2)
  a_weight_linear <- a_num / a_den
  # Normalize so A(1000 Hz) = 1
  f1k <- 1000^2
  a_1k <- 12194^2 * f1k^2 / ((f1k + 20.6^2) * sqrt((f1k + 107.7^2) * (f1k + 737.9^2)) * (f1k + 12194^2))
  a_weight_linear <- a_weight_linear / a_1k

  # Apply A-weighting to power spectrum (power = linear^2)
  a_weight_power <- a_weight_linear^2
  weighted_mat <- mat * a_weight_power  # broadcasts over columns

  # Sum weighted power per frame → total A-weighted power
  a_power <- colSums(weighted_mat)

  # Convert to dB
  max_power <- max(a_power[a_power > 0], na.rm = TRUE)
  if (max_power <= 0) return(periodicity)  # all silent, leave as-is

  a_db <- 10 * log10(pmax(a_power, .Machine$double.eps) / max_power)

  # Interpolate spectrogram frames to CREPE frames
  spec_times <- as.numeric(colnames(mat))
  crepe_times <- (seq_len(n_frames) - 1) * (hop_length / sr)

  a_db_interp <- stats::approx(spec_times, a_db, xout = crepe_times,
                                rule = 2)$y

  # Zero periodicity for silent frames
  silent <- a_db_interp < silence_db
  periodicity[silent] <- 0.0

  periodicity
}


# Fallback: simple per-frame RMS silence detection
.crepe_silence_simple <- function(periodicity, audio_float, sr, hop_length,
                                   silence_db, n_frames) {
  n_samples <- length(audio_float)

  # Compute per-frame RMS
  rms <- numeric(n_frames)
  pad <- 512L  # match CREPE padding
  for (t in seq_len(n_frames)) {
    start <- (t - 1L) * hop_length + 1L
    end_idx <- min(start + 1023L, n_samples)
    if (start > n_samples) {
      rms[t] <- 0
    } else {
      frame <- audio_float[start:end_idx]
      rms[t] <- sqrt(mean(frame^2))
    }
  }

  # Convert to dB relative to max
  max_rms <- max(rms, na.rm = TRUE)
  if (max_rms <= 0) return(periodicity)

  rms_db <- 20 * log10(pmax(rms, .Machine$double.eps) / max_rms)

  silent <- rms_db < silence_db
  periodicity[silent] <- 0.0

  periodicity
}


# Write minimal 16-bit PCM WAV file
.write_minimal_wav <- function(audio_float, sr, path) {
  # Convert to 16-bit signed integers
  samples <- as.integer(round(pmin(pmax(audio_float, -1), 1) * 32767))

  n <- length(samples)
  data_bytes <- n * 2L  # 16-bit = 2 bytes per sample
  byte_rate <- as.integer(sr) * 2L
  file_size <- 36L + data_bytes

  con <- file(path, "wb")
  on.exit(close(con), add = TRUE)

  # RIFF header
  writeBin(charToRaw("RIFF"), con)
  writeBin(as.integer(file_size), con, size = 4, endian = "little")
  writeBin(charToRaw("WAVE"), con)

  # fmt chunk
  writeBin(charToRaw("fmt "), con)
  writeBin(16L, con, size = 4, endian = "little")       # chunk size
  writeBin(1L, con, size = 2, endian = "little")         # PCM format
  writeBin(1L, con, size = 2, endian = "little")         # mono
  writeBin(as.integer(sr), con, size = 4, endian = "little")  # sample rate
  writeBin(byte_rate, con, size = 4, endian = "little")  # byte rate
  writeBin(2L, con, size = 2, endian = "little")         # block align
  writeBin(16L, con, size = 2, endian = "little")        # bits per sample

  # data chunk
  writeBin(charToRaw("data"), con)
  writeBin(data_bytes, con, size = 4, endian = "little")
  writeBin(samples, con, size = 2, endian = "little")
}


# ── Internal: NaN-aware running mean filter ──────────────────────────────────
#
# Applies a centered running mean of width k, ignoring NA values.
# Returns NA for frames where all values in the window are NA.
#
.nanmean_filter <- function(x, k) {
  n <- length(x)
  if (n == 0 || k <= 1) return(x)

  half <- (k - 1L) %/% 2L
  result <- numeric(n)

  for (i in seq_len(n)) {
    lo <- max(1L, i - half)
    hi <- min(n, i + half)
    window_vals <- x[lo:hi]
    non_na <- window_vals[!is.na(window_vals)]
    if (length(non_na) == 0) {
      result[i] <- NA_real_
    } else {
      result[i] <- mean(non_na)
    }
  }
  result
}

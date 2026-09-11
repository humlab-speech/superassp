
#' Track fundamental frequency using SwiftF0 (ONNX)
#'
#' Estimates F0 by applying a convolutional neural network directly to a
#' short-time spectrogram computed inside the ONNX graph
#' (SwiftF0; \insertCite{nieradzik2025swiftf0}{superassp}). SwiftF0 targets
#' real-time speed (~130 ms for 5 s of audio on CPU) at accuracy competitive
#' with CREPE. No Python or PyTorch required — inference uses ONNX Runtime.
#'
#' @param listOfFiles Character vector of audio file paths. Any format
#'   supported by \pkg{av} is accepted.
#' @param beginTime Numeric. Start of analysis window in seconds. Default 0
#'   (file start).
#' @param endTime Numeric. End of analysis window in seconds. Default 0
#'   (file end).
#' @param minF Numeric. Minimum F0 in Hz to treat as voiced. Default 75 Hz
#'   (speech). Must be >= 46.875 Hz (model minimum; G1). For music, use
#'   46.875.
#' @param maxF Numeric. Maximum F0 in Hz to treat as voiced. Default 400 Hz
#'   (speech). Must be <= 2093.75 Hz (model maximum; C7). For music, use
#'   2093.75.
#' @param confidence_threshold Numeric (0–1). Frames with model confidence
#'   below this value are marked unvoiced. Default 0.9.
#' @param toFile Logical. If \code{TRUE}, write SSFF output files and return
#'   the count written (invisibly). If \code{FALSE}, return an
#'   \code{AsspDataObj}. Default \code{TRUE}.
#' @param explicitExt Character. Output file extension. Default \code{"sf0"}.
#' @param outputDirectory Character. Directory for output files. \code{NULL}
#'   (default) writes alongside the input file.
#' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
#'
#' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
#'   \describe{
#'     \item{\code{f0}}{REAL32, Hz, \emph{n_frames} × 1. Fundamental
#'       frequency; 0 in unvoiced frames.}
#'     \item{\code{confidence}}{REAL32, 0–1, \emph{n_frames} × 1. Model
#'       voicing confidence.}
#'   }
#'   Frame rate: fixed 62.5 Hz (16 ms hop; not configurable — fixed by the
#'   ONNX graph's internal STFT). If \code{toFile = TRUE}: integer count of
#'   files written, returned invisibly.
#'
#' @details
#' ONNX Runtime is installed automatically on first use (~30 MB, cached in
#' the R user directory). The model file (~400 KB) is downloaded from the
#' \href{https://huggingface.co/FredrikKarlssonSpeech/swift-f0-onnx}{
#' swift-f0-onnx Hugging Face Hub repo} on first use (requires the
#' \pkg{huggingfaceR} package and a network connection) and cached in the R
#' user directory; subsequent calls read the cached copy with no network
#' access.
#'
#' Pre-processing (fixed by model training): resample to 16 kHz mono,
#' normalise to full-scale float32 in \verb{[-1, 1]}. Framing, STFT (1024-sample
#' window, 256-sample hop), and pitch/confidence estimation all happen
#' inside the ONNX graph — no manual windowing is needed. Post-processing
#' (applied outside the graph, matching the reference implementation): a
#' frame is voiced when \code{confidence > confidence_threshold} AND its
#' pitch estimate falls within \code{[minF, maxF]}; unvoiced frames get
#' \code{f0 = 0}.
#'
#' @examples
#' \dontrun{
#' trk_pitch_swiftf0(
#'   system.file("samples", "sustained", "a1.wav", package = "superassp"),
#'   toFile = FALSE
#' )
#' }
#' @usage trk_pitch_swiftf0(listOfFiles, beginTime = 0, endTime = 0, minF = 75, maxF = 400, confidence_threshold = 0.9, toFile = TRUE, explicitExt = "sf0", outputDirectory = NULL, verbose = TRUE)
#' @export
#'
#' @references
#'   \insertAllCited{}
#'
trk_pitch_swiftf0 <- function(listOfFiles,
                              beginTime             = 0.0,
                              endTime               = 0.0,
                              minF                  = 75.0,
                              maxF                  = 400.0,
                              confidence_threshold  = 0.9,
                              toFile                = TRUE,
                              explicitExt           = "sf0",
                              outputDirectory       = NULL,
                              verbose               = TRUE) {

  # -- Guards ------------------------------------------------------------------
  ensure_onnx()

  minF <- as.numeric(minF)
  maxF <- as.numeric(maxF)
  if (minF < 46.875 || maxF > 2093.75 || minF >= maxF) {
    cli::cli_abort(c(
      "{.arg minF}/{.arg maxF} out of range.",
      "i" = "Model supports 46.875-2093.75 Hz; got minF={minF}, maxF={maxF}."
    ))
  }

  confidence_threshold <- as.numeric(confidence_threshold)
  if (confidence_threshold < 0 || confidence_threshold > 1) {
    cli::cli_abort("{.arg confidence_threshold} must be between 0 and 1, got {confidence_threshold}.")
  }

  model_path <- .hf_get_cached_model(
    repo_id   = "FredrikKarlssonSpeech/swift-f0-onnx",
    filename  = "onnx/model.onnx",
    subdir    = "swift-f0",
    revision  = "ab37191652a37a99cef118bc37715e20048e0fb4"
  )

  if (length(listOfFiles) > 1L && !toFile) {
    cli::cli_abort("{.arg toFile = FALSE} only permitted for a single file.")
  }

  n_files   <- length(listOfFiles)
  .warn_if_lossy_input(listOfFiles)
  beginTime <- fast_recycle_times(beginTime, n_files)
  endTime   <- fast_recycle_times(endTime,   n_files)

  missing_files <- !file.exists(listOfFiles)
  if (any(missing_files)) {
    cli::cli_abort("File(s) not found: {.file {listOfFiles[missing_files]}}")
  }

  # -- ORT session (reused across files) ---------------------------------------
  session <- ort_session(model_path)
  hop     <- 256L
  sample_rate_ssff <- 16000.0 / hop

  # -- Per-file loop ------------------------------------------------------------
  outListOfFiles <- character(0L)
  outDataObj     <- NULL

  for (i in seq_along(listOfFiles)) {
    origSoundFile <- normalizePath(listOfFiles[[i]], mustWork = TRUE)
    bt <- beginTime[i]
    et <- endTime[i]

    if (verbose) {
      cli::cli_inform("SwiftF0: {basename(origSoundFile)}")
    }

    tryCatch({
      # -- Load audio at 16 kHz ---------------------------------------------
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
      audio_float <- as.numeric(audio_data) / 2147483647.0  # INT32_MAX -> [-1, 1]
      n_samples   <- length(audio_float)

      if (n_samples < hop) {
        cli::cli_warn(
          "SwiftF0: {.file {basename(origSoundFile)}} too short (< {hop} samples), skipping."
        )
        next
      }

      # -- ONNX inference (STFT + CNN run inside the graph) ------------------
      result <- ort_run(
        session,
        inputs       = list(input_audio = audio_float),
        shapes       = list(c(1L, as.integer(n_samples))),
        output_names = c("pitch_hz", "confidence")
      )
      pitch_hz   <- as.numeric(result$pitch_hz)
      confidence <- as.numeric(result$confidence)
      n_frames   <- length(pitch_hz)

      # -- Post-processing (voicing decision outside the graph) -------------
      voiced <- confidence > confidence_threshold & pitch_hz >= minF & pitch_hz <= maxF
      f0     <- ifelse(voiced, pitch_hz, 0.0)

      # -- Build AsspDataObj -------------------------------------------------
      start_time_ssff <- 1.0 / sample_rate_ssff

      outDataObj <- list()
      attr(outDataObj, "trackFormats") <- c("REAL32", "REAL32")
      attr(outDataObj, "sampleRate")   <- sample_rate_ssff
      attr(outDataObj, "origFreq")     <- 16000.0
      attr(outDataObj, "startTime")    <- start_time_ssff
      attr(outDataObj, "startRecord")  <- 1L
      attr(outDataObj, "endRecord")    <- as.integer(n_frames)
      class(outDataObj) <- "AsspDataObj"
      AsspFileFormat(outDataObj) <- "SSFF"
      AsspDataFormat(outDataObj) <- 2L

      outDataObj <- addTrack(outDataObj, "f0", matrix(f0, ncol = 1), "REAL32")
      outDataObj <- addTrack(outDataObj, "confidence", matrix(confidence, ncol = 1), "REAL32")

      # -- Output -----------------------------------------------------------
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
        "SwiftF0 failed for {.file {basename(origSoundFile)}}: {conditionMessage(e)}"
      )
    })
  }

  if (toFile) invisible(length(outListOfFiles)) else outDataObj
}

# -- Function attributes ------------------------------------------------------
attr(trk_pitch_swiftf0, "ext")             <- "sf0"
attr(trk_pitch_swiftf0, "tracks")          <- c("f0", "confidence")
attr(trk_pitch_swiftf0, "outputType")      <- "SSFF"
attr(trk_pitch_swiftf0, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")

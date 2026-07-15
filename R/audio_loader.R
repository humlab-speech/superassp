#' Uniform audio-loading entry point for DSP wrappers
#'
#' Single helper used by `trk_*` and `lst_*` wrappers to obtain audio data
#' for a downstream DSP framework. Encapsulates the package-wide fallback
#' contract:
#'
#' \enumerate{
#'   \item Try the framework's native loader (libassp via
#'     [processMediaFiles_LoadAndProcess()], pladdrr via
#'     [av_load_for_pladdrr()], etc.).
#'   \item On failure, on unsupported native format, or on any
#'     file-format-not-supported condition, fall back to [read_audio()] to
#'     obtain an [AsspDataObj] and adapt that to the framework's
#'     in-memory representation.
#'   \item If the framework cannot consume an [AsspDataObj] at all (some
#'     external binaries, file-path-only APIs), dump the AsspDataObj to a
#'     temporary WAV file via the package's WAV writer and pass the path.
#' }
#'
#' This helper is the **only** loader new DSP wrappers should call. Existing
#' wrappers may still call the lower-level loaders directly during the
#' transition period.
#'
#' @param file Path to the audio or media file.
#' @param begin Start of region to read (seconds, or samples if
#'   `samples = TRUE`). Default 0 = file start.
#' @param end End of region to read (seconds, or samples if
#'   `samples = TRUE`). Default 0 = file end.
#' @param samples Logical. If TRUE, `begin`/`end` are sample indices.
#' @param framework One of `"assp"`, `"sptk"`, `"snack"`, `"pladdrr"`,
#'   `"raw"`. Selects the native loader to try first; `"raw"` skips the
#'   native attempt and goes straight to [read_audio()].
#' @param ... Framework-specific arguments forwarded to the native loader.
#'
#' @return The framework's native audio container:
#'   \itemize{
#'     \item `assp`/`sptk`/`snack`/`raw` — an [AsspDataObj].
#'     \item `pladdrr` — a pladdrr `Sound` (R6) object.
#'   }
#'
#' @keywords internal
#' @noRd
assp_load_audio_for_dsp <- function(file,
                                    begin = 0,
                                    end = 0,
                                    samples = FALSE,
                                    framework = c("assp", "sptk", "snack",
                                                  "pladdrr", "raw"),
                                    ...) {
  framework <- match.arg(framework)

  # Lossy-input warning (design goal: applying a DSP routine to a lossy-encoded
  # signal presents a warning). Reuses knownLossless() so behavior matches the
  # batch path in processMediaFiles_LoadAndProcess(). Deduped once per file per
  # session via .frequency so parallel workers / repeated calls do not spam.
  # Lives here (the DSP-loading helper) rather than in read_audio(), which is
  # raw I/O a user may intentionally call on lossy audio.
  .warn_if_lossy_input(file)

  if (framework == "pladdrr") {
    # pladdrr Sound: use the dedicated loader (already implements the
    # native-then-av-transcode fallback contract).
    start_sec <- if (samples && begin > 0) .samples_to_seconds(file, begin) else begin
    end_sec   <- if (samples && end   > 0) .samples_to_seconds(file, end)   else end
    return(av_load_for_pladdrr(file_path = file,
                               start_time = start_sec,
                               end_time   = end_sec,
                               ...))
  }

  # All AsspDataObj-consuming frameworks share the same fallback chain:
  # read_audio() already implements (libassp native) -> (av fallback) with
  # sample-accurate windowing for variable-rate containers.
  obj <- read_audio(fname = file, begin = begin, end = end, samples = samples)

  # DSP kernels (e.g. SPTK REAPER, which casts the waveform to int16_t) assume
  # int16-range samples. read_audio() returns raw PCM at the source bit depth,
  # so a 24/32-bit file overflows the cast. Normalize to int16 here — the DSP
  # boundary — leaving the public read_audio()/read_ssff() faithful.
  .normalize_dsp_audio_int16(obj)
}

#' Downscale >16-bit integer audio to the int16 range.
#'
#' AsspDataObj audio tracks feeding int16-based DSP kernels must be in int16
#' range. INT24/INT32 sources are downscaled so every framework sees consistent
#' magnitude (mirrors the /65536 normalization the av fallback applies).
#'
#' @param obj An AsspDataObj (or NULL).
#' @return The object with its `audio` track normalized to INT16, unchanged for
#'   INT16 / non-audio inputs.
#' @keywords internal
#' @noRd
.normalize_dsp_audio_int16 <- function(obj) {
  if (is.null(obj) || !("audio" %in% names(obj))) return(obj)
  tf    <- attr(obj, "trackFormats")
  fmt   <- if (length(tf)) as.character(tf[[1]]) else NA_character_
  shift <- switch(fmt, INT32 = 65536, INT24 = 256, 1)
  if (shift != 1) {
    a <- round(obj[["audio"]] / shift)
    storage.mode(a) <- "integer"
    obj[["audio"]] <- a
    attr(obj, "trackFormats")[1] <- "INT16"
  }
  obj
}

#' Warn once per file when a DSP wrapper is handed a lossy-encoded input.
#'
#' Extension-based check against [knownLossless()] — the same source of truth
#' the batch path (`processMediaFiles_LoadAndProcess()`) uses, keeping behavior
#' consistent. Emitted unconditionally (design goal) but deduplicated to once
#' per file per session via cli/rlang `.frequency` so parallel workers and
#' repeated calls on the same recording do not spam the console.
#'
#' @param file One or more file paths.
#' @return Invisibly `NULL`; called for its side effect (a warning).
#' @keywords internal
#' @noRd
.warn_if_lossy_input <- function(file) {
  paths <- as.character(file)
  paths <- paths[nzchar(paths)]
  if (length(paths) == 0L) return(invisible(NULL))

  known_lossless <- knownLossless()
  for (p in unique(paths)) {
    ext <- tolower(tools::file_ext(p))
    if (nzchar(ext) && !(ext %in% tolower(known_lossless))) {
      cli::cli_warn(
        c(
          "!" = "{.file {basename(p)}} is in a lossy-compressed format ({.val {ext}}).",
          "i" = "Lossy compression may reduce DSP accuracy.",
          "x" = "For faithful analysis, use a lossless format (e.g. {.val {c('wav', 'flac')}})."
        ),
        .frequency = "once",
        .frequency_id = paste0("superassp_lossy_", normalizePath(p, mustWork = FALSE))
      )
    }
  }
  invisible(NULL)
}

#' Convert a sample index to seconds via av media metadata.
#' @keywords internal
#' @noRd
.samples_to_seconds <- function(file, n_samples) {
  info <- tryCatch(media_info(file), error = function(e) NULL)
  if (is.null(info) || is.null(info$audio)) {
    cli::cli_abort("Could not determine sample rate for sample-based indexing of {.file {basename(file)}}.")
  }
  n_samples / info$audio$sample_rate
}

#' Pre-fetch a media manifest for batch processing.
#'
#' Builds a one-row-per-file data frame summarising every input to a batch
#' DSP wrapper. Uses [media_info()] (cached) so subsequent worker invocations
#' for the same file do not re-probe FFmpeg. The native-vs-av decision is
#' encoded once here so workers can skip the fast-path attempt for known
#' non-native formats.
#'
#' @param listOfFiles Character vector of file paths.
#' @return A data frame with one row per file and columns
#'   `path` (input, untouched), `normalized_path`, `exists`, `native_ext`
#'   (logical), `duration` (seconds, NA if unknown), `sample_rate` (Hz),
#'   `channels`, `codec` (audio codec name).
#' @keywords internal
#' @noRd
build_media_manifest <- function(listOfFiles) {
  native_exts <- c("wav", "au", "kay", "nist", "nsp", "aiff", "aif")
  paths <- as.character(listOfFiles)
  n <- length(paths)
  if (n == 0L) {
    return(data.frame(path = character(0), normalized_path = character(0),
                      exists = logical(0), native_ext = logical(0),
                      duration = numeric(0), sample_rate = numeric(0),
                      channels = integer(0), codec = character(0),
                      stringsAsFactors = FALSE))
  }

  normalized <- normalizePath(path.expand(paths), winslash = "/", mustWork = FALSE)
  exists_vec <- file.exists(normalized)
  exts       <- tolower(tools::file_ext(paths))
  native_vec <- exts %in% native_exts

  duration <- rep(NA_real_, n)
  sr       <- rep(NA_real_, n)
  channels <- rep(NA_integer_, n)
  codec    <- rep(NA_character_, n)

  for (i in seq_len(n)) {
    if (!exists_vec[i]) next
    info <- tryCatch(media_info(normalized[i]), error = function(e) NULL)
    if (is.null(info)) next
    duration[i] <- if (is.null(info$duration)) NA_real_ else info$duration
    if (length(info$audio) > 0L) {
      sr[i]       <- info$audio$sample_rate
      channels[i] <- as.integer(info$audio$channels)
      codec[i]    <- info$audio$codec
    }
  }

  data.frame(
    path            = paths,
    normalized_path = normalized,
    exists          = exists_vec,
    native_ext      = native_vec,
    duration        = duration,
    sample_rate     = sr,
    channels        = channels,
    codec           = codec,
    stringsAsFactors = FALSE
  )
}

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
  read_audio(fname = file, begin = begin, end = end, samples = samples)
}

#' Convert a sample index to seconds via av media metadata.
#' @keywords internal
#' @noRd
.samples_to_seconds <- function(file, n_samples) {
  info <- tryCatch(media_info(file), error = function(e) NULL)
  if (is.null(info) || is.null(info$audio)) {
    stop("Could not determine sample rate for sample-based indexing of '",
         basename(file), "'", call. = FALSE)
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

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
  info <- tryCatch(av::av_media_info(file), error = function(e) NULL)
  if (is.null(info) || is.null(info$audio)) {
    stop("Could not determine sample rate for sample-based indexing of '",
         basename(file), "'", call. = FALSE)
  }
  n_samples / info$audio$sample_rate
}

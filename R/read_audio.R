#' Read an audio file into an AsspDataObj
#'
#' Unified audio reader. Tries the ASSP C-level reader first (fastest, supports
#' WAV, AU, NIST, SSFF and other native formats). Falls back to the \code{av}
#' package for any format av/FFmpeg can decode (MP3, MP4, FLAC, OGG, etc.).
#'
#' For variable-rate encoded formats (e.g. MP3) with \code{samples = TRUE},
#' sample positions are approximated via the nominal sample rate reported by the
#' container.
#'
#' @param fname Path to the audio or SSFF file.
#' @param begin Start of region to read. Default 0 = file start.
#'   In seconds (default) or samples (when \code{samples = TRUE}).
#' @param end End of region to read. Default 0 = file end.
#'   In seconds or samples.
#' @param samples Logical. If \code{TRUE}, \code{begin} and \code{end} are
#'   interpreted as sample indices.
#' @return An \code{AsspDataObj}.
#' @seealso \code{\link{read_ssff}} for SSFF-only files;
#'   \code{\link{read_json_track}} for JSTF files.
#' @export
read_audio <- function(fname, begin = 0, end = 0, samples = FALSE) {

  # --- Path 1: ASSP C-level (fastest, native formats) ---
  result <- tryCatch(
    read_ssff(fname, begin = begin, end = end, samples = samples),
    error = function(e) NULL
  )
  if (!is.null(result)) return(result)

  # --- Path 2: av package fallback (MP3, MP4, FLAC, OGG, ...) ---
  if (!requireNamespace("av", quietly = TRUE)) {
    stop("Could not read '", basename(fname), "' via ASSP C reader, and the ",
         "'av' package is not available as a fallback. Install av: ",
         "install.packages('av')")
  }

  # Convert sample indices to times if needed
  start_sec <- begin
  end_sec   <- end

  if (samples && (begin > 0 || end > 0)) {
    info <- tryCatch(av::av_media_info(fname), error = function(e) NULL)
    if (is.null(info) || is.null(info$audio)) {
      stop("Could not determine sample rate for sample-based indexing of '",
           basename(fname), "'")
    }
    sr        <- info$audio$sample_rate
    start_sec <- begin / sr
    end_sec   <- if (end == 0) 0 else end / sr
  }

  av_to_asspDataObj(
    fname,
    start_time = start_sec,
    end_time   = if (end_sec == 0) NULL else end_sec
  )
}

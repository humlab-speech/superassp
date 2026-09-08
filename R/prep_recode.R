#' Re-encode Media File with Custom Parameters
#'
#' Re-encodes any media file (audio/video) supported by the av package into
#' 16-bit PCM WAV, with optional resampling, channel remixing, and time
#' windowing. Returns the audio data in the same format as
#' \code{av::read_audio_bin}.
#'
#' Re-encoding goes through a temporary WAV file (\code{av::av_audio_convert()}
#' followed by \code{av::read_audio_bin()}); the temp file is removed on exit.
#' It's useful for:
#' \itemize{
#'   \item Converting sample rates for analysis
#'   \item Extracting audio from video files
#'   \item Time-windowing large files
#'   \item Remixing channel count across a corpus
#' }
#'
#' @param listOfFiles Character vector of file paths to media files
#' @param codec Either \code{"none"} (read the file as-is, no re-encoding) or
#'   \code{"pcm_s16le"} (re-encode to 16-bit PCM WAV). Required. These are the
#'   only two re-encoding needs superassp has internally; for anything else
#'   call \code{av::av_audio_convert()} directly.
#' @param sample_rate Target sample rate in Hz (default: NULL keeps original)
#' @param bit_rate Ignored for \code{"pcm_s16le"} (lossless); kept for
#'   interface symmetry with \code{av::av_audio_convert()}.
#' @param start_time Start time in seconds (default: NULL = start of file)
#' @param end_time End time in seconds (default: NULL = end of file)
#' @param channels Number of output channels: 1 (mono), 2 (stereo), or NULL (keep original)
#' @param verbose Logical; show progress messages (default: TRUE)
#' @param ... Additional arguments passed to \code{av::av_audio_convert}
#'
#' @return For single file: Integer vector with audio samples in s32le format
#'   (32-bit signed integers), with attributes:
#'   \itemize{
#'     \item \code{channels}: Number of audio channels (integer)
#'     \item \code{sample_rate}: Sample rate in Hz (integer)
#'   }
#'
#'   For multiple files: List of integer vectors, one per file
#'
#'   This matches the format returned by \code{av::read_audio_bin()}.
#'
#' @references
#' \insertRef{av2024}{superassp}
#'
#' \insertRef{ffmpeg2024}{superassp}
#'
#' @seealso
#' \code{\link[av]{av_audio_convert}}, \code{\link[av]{read_audio_bin}},
#' \code{\link{av_to_asspDataObj}}
#'
#' @examples
#' \dontrun{
#' # Read as-is
#' audio <- prep_recode("speech.wav", codec = "none")
#'
#' # Extract segment from 1-3 seconds
#' audio_segment <- prep_recode("long.wav",
#'                               codec = "pcm_s16le",
#'                               start_time = 1.0,
#'                               end_time = 3.0)
#'
#' # Downsample to 16 kHz
#' audio_16k <- prep_recode("high_res.wav",
#'                          codec = "pcm_s16le",
#'                          sample_rate = 16000)
#'
#' # Convert to mono
#' audio_mono <- prep_recode("stereo.wav",
#'                           codec = "pcm_s16le",
#'                           channels = 1)
#'
#' # Batch processing
#' files <- c("file1.mp4", "file2.wav", "file3.flac")
#' audio_list <- prep_recode(files,
#'                           codec = "pcm_s16le",
#'                           sample_rate = 44100,
#'                           channels = 1)
#'
#' # Access audio data (same as av::read_audio_bin)
#' audio <- prep_recode("test.wav", codec = "pcm_s16le")
#' cat("Channels:", attr(audio, "channels"), "\n")
#' cat("Sample rate:", attr(audio, "sample_rate"), "\n")
#' cat("Duration:", length(audio) / attr(audio, "channels") / attr(audio, "sample_rate"), "s\n")
#' }
#'
#' @keywords internal
prep_recode <- function(listOfFiles,
                        codec,
                        sample_rate = NULL,
                        bit_rate = NULL,
                        start_time = NULL,
                        end_time = NULL,
                        channels = NULL,
                        verbose = TRUE,
                        ...) {

  # Check av package
  if (!requireNamespace("av", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg av} is required but not installed.")
  }

  # Validate codec
  if (missing(codec) || is.null(codec) || codec == "") {
    cli::cli_abort("codec argument is required ('none' or 'pcm_s16le')")
  }
  if (!codec %in% c("none", "direct", "pcm_s16le")) {
    cli::cli_abort(c(
      "Unsupported codec {.val {codec}}.",
      "i" = "prep_recode() only supports {.val none} and {.val pcm_s16le}.",
      "i" = "For other codecs, call {.fun av::av_audio_convert} directly."
    ))
  }

  # Normalize parameters
  n_files <- length(listOfFiles)
  if (!is.null(start_time)) start_time <- rep_len(start_time, n_files)
  if (!is.null(end_time)) end_time <- rep_len(end_time, n_files)

  # Progress bar for multiple files
  if (verbose && n_files > 1) {
    cli::cli_alert_info("Re-encoding {n_files} file{?s} with codec {codec}")
    pb <- cli::cli_progress_bar("Re-encoding", total = n_files)
  }

  # Process files
  results <- vector("list", n_files)

  # Track if any files succeeded (for single file return)
  any_success <- FALSE

  for (i in seq_along(listOfFiles)) {
    file_path <- listOfFiles[i]

    # Validate file exists
    if (!file.exists(file_path)) {
      cli::cli_warn("File not found: {.file {file_path}}")
      results[[i]] <- NULL
      if (verbose && n_files > 1) cli::cli_progress_update()
      next
    }

    # Get file info with error handling
    info <- tryCatch({
      media_info(file_path)
    }, error = function(e) {
      # FFMPEG error - invalid file
      cli::cli_warn("Invalid media file: {.file {basename(file_path)}} ({e$message})")
      return(NULL)
    })

    if (is.null(info)) {
      results[[i]] <- NULL
      if (verbose && n_files > 1) cli::cli_progress_update()
      next
    }

    if (length(info$audio) == 0) {
      cli::cli_warn("No audio stream found in: {.file {basename(file_path)}}")
      results[[i]] <- NULL
      if (verbose && n_files > 1) cli::cli_progress_update()
      next
    }

    tryCatch({
      audio_info <- info$audio

      # Determine start and end times for this file
      file_start <- if (!is.null(start_time)) start_time[i] else NULL
      file_end <- if (!is.null(end_time)) end_time[i] else NULL

      # Determine target parameters
      target_sr <- if (!is.null(sample_rate)) sample_rate else audio_info$sample_rate
      target_ch <- if (!is.null(channels)) channels else audio_info$channels

      # Re-encoding is required whenever the user specified a codec, any
      # output sample-rate / channel / bit-rate constraint, or a time window.
      needs_recode <- !(codec %in% c("none", "direct")) ||
        (!is.null(sample_rate) && sample_rate != audio_info$sample_rate) ||
        (!is.null(channels)    && channels    != audio_info$channels) ||
        !is.null(bit_rate) ||
        !is.null(file_start)  || !is.null(file_end)

      # If codec is "none" or "direct", read directly without re-encoding
      if (codec %in% c("none", "direct") && is.null(file_start) && is.null(file_end) &&
          is.null(sample_rate) && is.null(channels)) {
        if (verbose && n_files == 1) {
          cli::cli_alert_info("Reading directly (no re-encoding needed)")
        }

        audio_data <- av::read_audio_bin(
          file_path,
          channels = target_ch,
          sample_rate = target_sr
        )

        results[[i]] <- audio_data
        any_success <- TRUE

      } else {
        # Need to re-encode or apply time windowing. av (CRAN) has no
        # in-memory transcode, so go through a temp WAV: av_audio_convert()
        # writes it, read_audio_bin() reads it back in the target format.
        if (verbose && n_files == 1) {
          cli::cli_alert_info("Transcoding to 16-bit PCM WAV")
        }

        tmp_wav <- tempfile(fileext = ".wav")

        convert_args <- list(
          audio = file_path,
          output = tmp_wav,
          channels = target_ch,
          sample_rate = target_sr,
          verbose = FALSE
        )

        # Add optional arguments
        if (!is.null(bit_rate)) {
          convert_args$bit_rate <- bit_rate
        }
        if (!is.null(file_start)) convert_args$start_time <- file_start
        if (!is.null(file_end)) {
          # Calculate duration if start_time specified
          if (!is.null(file_start)) {
            convert_args$total_time <- file_end - file_start
          } else {
            convert_args$total_time <- file_end
          }
        }

        # Add any additional arguments
        extra_args <- list(...)
        if (length(extra_args) > 0) {
          convert_args <- c(convert_args, extra_args)
        }

        do.call(av::av_audio_convert, convert_args)
        audio_data <- av::read_audio_bin(tmp_wav, channels = target_ch, sample_rate = target_sr)
        unlink(tmp_wav)

        results[[i]] <- audio_data
        any_success <- TRUE
      }

    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(file_path)}}: {e$message}")
      results[[i]] <- NULL
    })

    if (verbose && n_files > 1) cli::cli_progress_update()
  }

  if (verbose && n_files > 1) cli::cli_progress_done()

  # Return results
  if (n_files == 1) {
    # For single file, return the result directly (or NULL if failed)
    if (length(results) == 0) {
      return(NULL)
    }
    return(results[[1]])
  } else {
    # For multiple files, return list
    return(results)
  }
}

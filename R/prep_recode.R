#' Re-encode Media File with Custom Parameters
#'
#' Re-encodes any media file (audio/video) supported by the av package into a
#' specified format with custom codec, sample rate, bit rate, and optional time
#' windowing. Returns the audio data in the same format as \code{av::read_audio_bin}.
#'
#' This function performs in-memory transcoding using \code{av::av_audio_transcode()},
#' avoiding intermediate files on disk. It's useful for:
#' \itemize{
#'   \item Converting sample rates for analysis
#'   \item Extracting audio from video files
#'   \item Time-windowing large files
#'   \item Normalizing formats across a corpus
#'   \item Testing codec-specific effects
#' }
#'
#' @param listOfFiles Character vector of file paths to media files
#' @param codec Output codec (e.g., "pcm_s16le", "mp3", "flac", "vorbis"). Required.
#'   See \code{av::av_encoders()} for available codecs.
#' @param sample_rate Target sample rate in Hz (default: NULL keeps original)
#' @param bit_rate Target bit rate for lossy codecs (default: NULL uses codec default).
#'   Specify as integer (bits/second), e.g., 128000, 192000, 320000
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
#' @details
#' **Supported Formats:**
#'
#' The av package supports a wide range of formats through FFmpeg:
#' \itemize{
#'   \item \bold{Lossless:} wav, flac, alac, ape, wv
#'   \item \bold{Lossy:} mp3, ogg, aac, opus, wma
#'   \item \bold{Video:} mp4, mkv, avi, mov, webm (extracts audio)
#' }
#'
#' **Common Codec Examples:**
#' \itemize{
#'   \item \bold{WAV:} "pcm_s16le" (16-bit), "pcm_s24le" (24-bit), "pcm_f32le" (32-bit float)
#'   \item \bold{MP3:} "libmp3lame"
#'   \item \bold{FLAC:} "flac"
#'   \item \bold{OGG:} "libvorbis"
#'   \item \bold{AAC:} "aac"
#'   \item \bold{OPUS:} "libopus"
#' }
#'
#' **Processing Strategy:**
#'
#' 1. If no re-encoding needed (no codec/sample_rate/channels change, no windowing):
#'    - Returns \code{av::read_audio_bin()} result directly
#'
#' 2. If re-encoding or windowing needed:
#'    - Uses \code{av::av_audio_transcode()} for in-memory transcoding
#'    - Returns audio data directly (no temporary files)
#'
#' **Performance:**
#' - Pure in-memory operation (no temporary files)
#' - Fast conversion for compatible codecs
#' - Time windowing reduces memory usage
#'
#' @references
#' \insertRef{av2024}{superassp}
#'
#' \insertRef{ffmpeg2024}{superassp}
#'
#' @seealso
#' \code{\link[av]{av_audio_transcode}}, \code{\link[av]{read_audio_bin}},
#' \code{\link{av_to_asspDataObj}}
#'
#' @examples
#' \dontrun{
#' # Basic usage - convert to WAV PCM
#' audio <- prep_recode("video.mp4", codec = "pcm_s16le")
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
#' # Convert to MP3 with specific bit rate
#' audio_mp3 <- prep_recode("speech.wav",
#'                          codec = "mp3",
#'                          bit_rate = 192000)
#'
#' # Convert to FLAC (lossless compression)
#' audio_flac <- prep_recode("recording.wav",
#'                           codec = "flac")
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
#' @export
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
    stop("Package 'av' is required but not installed.\n",
         "Install with: install.packages('av')",
         call. = FALSE)
  }

  # Validate codec
  if (missing(codec) || is.null(codec) || codec == "") {
    stop("codec argument is required (e.g., 'pcm_s16le', 'mp3', 'flac')", call. = FALSE)
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
      warning("File not found: ", file_path, call. = FALSE)
      results[[i]] <- NULL
      if (verbose && n_files > 1) cli::cli_progress_update()
      next
    }

    # Get file info with error handling
    info <- tryCatch({
      av::av_media_info(file_path)
    }, error = function(e) {
      # FFMPEG error - invalid file
      warning("Invalid media file: ", basename(file_path), " (",  e$message, ")", call. = FALSE)
      return(NULL)
    })

    if (is.null(info)) {
      results[[i]] <- NULL
      if (verbose && n_files > 1) cli::cli_progress_update()
      next
    }

    if (length(info$audio) == 0) {
      warning("No audio stream found in: ", basename(file_path), call. = FALSE)
      results[[i]] <- NULL
      if (verbose && n_files > 1) cli::cli_progress_update()
      next
    }

    tryCatch({
      audio_info <- info$audio

      # Determine start and end times for this file
      file_start <- if (!is.null(start_time)) start_time[i] else NULL
      file_end <- if (!is.null(end_time)) end_time[i] else NULL

      # Determine if we can read directly without re-encoding
      needs_recode <- FALSE

      # Always need to recode if codec is specified
      needs_recode <- TRUE

      # Check if sample rate conversion needed
      if (!is.null(sample_rate) && sample_rate != audio_info$sample_rate) {
        needs_recode <- TRUE
      }

      # Check if channel conversion needed
      if (!is.null(channels) && channels != audio_info$channels) {
        needs_recode <- TRUE
      }

      # Check if bit rate specified
      if (!is.null(bit_rate)) {
        needs_recode <- TRUE
      }

      # Check if time windowing needed
      if (!is.null(file_start) || !is.null(file_end)) {
        needs_recode <- TRUE
      }

      # Determine target parameters
      target_sr <- if (!is.null(sample_rate)) sample_rate else audio_info$sample_rate
      target_ch <- if (!is.null(channels)) channels else audio_info$channels

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
        # Need to re-encode or apply time windowing
        if (verbose && n_files == 1) {
          cli::cli_alert_info("Transcoding in-memory with codec {codec}")
        }

        # Build av_audio_transcode arguments
        transcode_args <- list(
          audio = file_path,
          codec = codec,
          channels = target_ch,
          sample_rate = target_sr,
          verbose = FALSE
        )

        # Add optional arguments
        if (!is.null(bit_rate)) {
          transcode_args$bit_rate <- bit_rate
        }
        if (!is.null(file_start)) transcode_args$start_time <- file_start
        if (!is.null(file_end)) {
          # Calculate duration if start_time specified
          if (!is.null(file_start)) {
            transcode_args$total_time <- file_end - file_start
          } else {
            transcode_args$total_time <- file_end
          }
        }

        # Add any additional arguments
        extra_args <- list(...)
        if (length(extra_args) > 0) {
          transcode_args <- c(transcode_args, extra_args)
        }

        # Perform in-memory transcoding
        audio_data <- do.call(av::av_audio_transcode, transcode_args)

        results[[i]] <- audio_data
        any_success <- TRUE
      }

    }, error = function(e) {
      warning("Error processing ", basename(file_path), ": ",
              e$message, call. = FALSE)
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

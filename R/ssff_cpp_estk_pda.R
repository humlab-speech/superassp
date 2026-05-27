##' Track fundamental frequency using the ESTk PDA algorithm
##'
##' Extracts F0 using the Edinburgh Speech Tools super-resolution Pitch Detection
##' Algorithm (PDA), a normalized cross-correlation tracker with DP-based peak
##' tracking optimized for speech with variable pitch. Offers sub-frame resolution
##' and explicit voiced/unvoiced transition thresholds.
##'
##' @inheritParams trk_acf
##' @param windowShift Numeric. Frame shift in milliseconds; sets output frame rate
##'   (1000 / windowShift Hz). Default 5.0 ms.
##' @param windowSize Numeric. Analysis window length in milliseconds. Default 10.0 ms.
##' @param minF Numeric. Minimum F0 in Hz. Default 40.0 Hz.
##' @param maxF Numeric. Maximum F0 in Hz. Default 400.0 Hz.
##' @param decimation Integer. Correlation decimation factor (higher = faster but coarser).
##'   Default 4.
##' @param noise_floor Numeric. Silence threshold; frames below this RMS are treated
##'   as unvoiced. Default 120.
##' @param min_v2uv_coef_thresh Numeric. Minimum correlation for voiced-to-unvoiced
##'   transition (0–1). Default 0.75.
##' @param v2uv_coef_thresh_ratio Numeric. Correlation ratio threshold for
##'   voiced-to-unvoiced transition (0–1). Default 0.85.
##' @param uv2v_coef_thresh Numeric. Correlation threshold for unvoiced-to-voiced
##'   transition (0–1). Default 0.88.
##' @param anti_doubling_thresh Numeric. Threshold to suppress F0 octave doublings
##'   (0–1). Default 0.77.
##' @param peak_tracking Logical. If \code{TRUE}, enable peak tracking for smoother
##'   F0 contours. Default \code{FALSE}.
##' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
##'   paths written invisibly. If \code{FALSE}, return an \code{AsspDataObj}.
##'   Default \code{FALSE}.
##' @param explicitExt Character. Output file extension. Default \code{"pda"}.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
##'   \describe{
##'     \item{\code{F0}}{REAL32, fundamental frequency in Hz, n_frames × 1.
##'       Zero indicates unvoiced frames.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 200 Hz).
##'   If \code{toFile = TRUE}: output file path(s), returned invisibly.
##'
##' @export
##' @references
##' \insertCite{Medan1991}{superassp}
##'
##' \insertCite{Bagshaw1993}{superassp}
##'
##' @examples
##' \dontrun{
##' # Extract F0 from audio file
##' f0_data <- trk_pitch_pda("recording.wav", toFile = FALSE)
##'
##' # Process with custom parameters
##' trk_pitch_pda("speech.mp3", minF = 75, maxF = 300, windowShift = 10)
##'
##' # Enable peak tracking for smoother contours
##' trk_pitch_pda("recording.wav", peak_tracking = TRUE, toFile = FALSE)
##'
##' # Process video file (extracts audio)
##' trk_pitch_pda("interview.mp4", toFile = FALSE)
##' }
trk_pitch_pda <- function(listOfFiles,
                         beginTime = 0.0,
                         endTime = 0.0,
                         windowShift = 5.0,
                         windowSize = 10.0,
                         minF = 40.0,
                         maxF = 400.0,
                         decimation = 4,
                         noise_floor = 120,
                         min_v2uv_coef_thresh = 0.75,
                         v2uv_coef_thresh_ratio = 0.85,
                         uv2v_coef_thresh = 0.88,
                         anti_doubling_thresh = 0.77,
                         peak_tracking = FALSE,
                         toFile = FALSE,
                         explicitExt = "pda",
                         outputDirectory = NULL,
                         verbose = TRUE) {

  # Validate inputs
  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No input files specified in {.arg listOfFiles}")
  }

  # Normalize paths
  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)

  # Check file existence
  files_exist <- file.exists(listOfFiles)
  if (!all(files_exist)) {
    missing_files <- listOfFiles[!files_exist]
    cli::cli_abort(c(
      "!" = "Some files do not exist:",
      "x" = "{.file {fast_basename(missing_files)}}"
    ))
  }

  n_files <- length(listOfFiles)

  # Normalize time parameters
  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime <- if (is.null(endTime)) 0.0 else endTime

  # Recycle time parameters
  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime) == 1) endTime <- rep(endTime, n_files)

  # Setup output directory
  if (toFile) {
    makeOutputDirectory(outputDirectory, FALSE, "trk_pitch_pda")
  }

  if (verbose) format_apply_msg("trk_pitch_pda", n_files, beginTime, endTime)

  # Process each file
  results <- vector("list", n_files)

  if (verbose && n_files > 1) {
    cli::cli_progress_bar(
      "Processing files",
      total = n_files,
      format = "{cli::pb_spin} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}"
    )
  }

  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[i]

    tryCatch({
      # Load audio with av
      audio_obj <- av_to_asspDataObj(
        file_path,
        start_time = bt,
        end_time = if (et == 0.0) NULL else et
      )

      # Call C++ ESTK PDA
      pda_result <- estk_pda_cpp(
        audio_obj = audio_obj,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        windowSize = windowSize,
        decimation = decimation,
        noise_floor = noise_floor,
        min_v2uv_coef_thresh = min_v2uv_coef_thresh,
        v2uv_coef_thresh_ratio = v2uv_coef_thresh_ratio,
        uv2v_coef_thresh = uv2v_coef_thresh,
        anti_doubling_thresh = anti_doubling_thresh,
        peak_tracking = peak_tracking,
        verbose = FALSE
      )

      # Convert to AsspDataObj with one track (F0)
      out_obj <- create_pda_asspobj(pda_result, windowShift)

      # Handle output
      if (toFile) {
        base_name <- tools::file_path_sans_ext(basename(file_path))
        out_dir <- if (is.null(outputDirectory)) dirname(file_path) else outputDirectory
        out_file <- file.path(out_dir, paste0(base_name, ".", explicitExt))
        write.AsspDataObj(out_obj, out_file)
        results[[i]] <- out_file
      } else {
        results[[i]] <- out_obj
      }

    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(file_path)}}: {conditionMessage(e)}")
      results[[i]] <- if (toFile) NULL else NULL
    })

    if (verbose && n_files > 1) {
      cli::cli_progress_update()
    }
  }

  if (verbose && n_files > 1) {
    cli::cli_progress_done()
  }

  # Return results
  if (toFile) {
    success_count <- sum(!sapply(results, is.null))
    if (verbose) {
      cli::cli_inform("Successfully processed {success_count}/{n_files} file{?s}")
    }
    return(invisible(if (n_files == 1) results[[1]] else results))
  } else {
    return(if (n_files == 1) results[[1]] else results)
  }
}


##' Helper function to create AsspDataObj from ESTK PDA results
##' @keywords internal
##' @noRd
create_pda_asspobj <- function(pda_result, windowShift) {
  n_frames <- pda_result$n_frames
  sample_rate <- pda_result$sample_rate
  frame_rate <- 1000.0 / windowShift  # windowShift is in ms

  # Create AsspDataObj structure with one track
  obj <- list()
  obj$F0 <- pda_result$f0

  # Set attributes matching wrassp format
  attr(obj, "trackFormats") <- c("REAL32")
  attr(obj, "sampleRate") <- frame_rate  # Frames per second
  attr(obj, "origFreq") <- as.numeric(sample_rate)  # Original audio sample rate
  attr(obj, "startTime") <- 0.0
  attr(obj, "startRecord") <- 1L
  attr(obj, "endRecord") <- as.integer(n_frames)
  attr(obj, "fileInfo") <- c(20L, 1L)  # SSFF format, 1 track
  class(obj) <- "AsspDataObj"

  return(obj)
}

# Set function attributes
attr(trk_pitch_pda, "ext") <- "pda"
attr(trk_pitch_pda, "tracks") <- c("F0")
attr(trk_pitch_pda, "outputType") <- "SSFF"
attr(trk_pitch_pda, "nativeFiletypes") <- c("wav")

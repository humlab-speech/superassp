##' YIN Pitch Tracking (C++ implementation)
##'
##' @description Extract F0 (fundamental frequency) using the YIN algorithm.
##'   This is a high-performance C++ implementation that is significantly faster
##'   than the Python/librosa version and requires no Python dependencies.
##'
##'   The YIN algorithm \insertCite{Cheveigné.2002.10.1121/1.1458024}{superassp}
##'   uses a modified autocorrelation approach with cumulative mean normalization
##'   to reliably detect pitch even in noisy conditions.
##'
##'   All input media formats are supported via the av package, including video
##'   files from which audio will be automatically extracted.
##'
##' @param listOfFiles Vector of file paths to process
##' @param beginTime Start time in seconds (default: 0.0)
##' @param endTime End time in seconds (default: 0.0 = end of file)
##' @param windowShift Frame shift in milliseconds (default: 5.0)
##' @param windowSize Window size in milliseconds (default: 30.0)
##' @param minF Minimum F0 in Hz (default: 70.0)
##' @param maxF Maximum F0 in Hz (default: 200.0)
##' @param threshold Voicing threshold (default: 0.1, lower = more permissive)
##' @param toFile Write results to file (default: FALSE)
##' @param explicitExt Output file extension (default: "yip")
##' @param outputDirectory Output directory (default: NULL = same as input)
##' @param verbose Show progress messages (default: TRUE)
##'
##' @return If toFile=TRUE, returns output file path(s) invisibly.
##'   If toFile=FALSE, returns AsspDataObj or list of AsspDataObj objects.
##'   Each object contains two tracks: "F0" (Hz) and "prob" (voicing probability).
##'
##' @export
##' @references
##' \insertAllCited{}
##'
##' @examples
##' \dontrun{
##' # Extract F0 from audio file
##' f0_data <- trk_yin("recording.wav", toFile = FALSE)
##'
##' # Process with custom parameters
##' trk_yin("speech.mp3", minF = 75, maxF = 300, windowShift = 10)
##'
##' # Process video file (extracts audio)
##' trk_yin("interview.mp4", toFile = FALSE)
##' }
trk_yin <- function(listOfFiles,
                    beginTime = 0.0,
                    endTime = 0.0,
                    windowShift = 5.0,
                    windowSize = 30.0,
                    minF = 70.0,
                    maxF = 200.0,
                    threshold = 0.1,
                    toFile = FALSE,
                    explicitExt = "yip",
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
    makeOutputDirectory(outputDirectory, FALSE, "trk_yin")
  }

  if (verbose) format_apply_msg("trk_yin", n_files, beginTime, endTime)

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
      audio_obj <- read_audio(
        file_path,
        begin = bt,
        end   = et
      )

      # Call C++ YIN
      yin_result <- yin_cpp(
        audio_obj = audio_obj,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        windowSize = windowSize,
        threshold = threshold,
        verbose = FALSE
      )

      # Convert to AsspDataObj with two tracks (F0 and probability)
      out_obj <- create_yin_asspobj(yin_result, windowShift)

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


##' Helper function to create AsspDataObj from YIN results
##' @keywords internal
##' @noRd
create_yin_asspobj <- function(yin_result, windowShift) {
  n_frames <- yin_result$n_frames
  sample_rate <- yin_result$sample_rate
  frame_rate <- 1000.0 / windowShift  # windowShift is in ms

  # Create AsspDataObj structure with two tracks
  obj <- list()
  obj$F0 <- yin_result$f0
  obj$prob <- yin_result$probability

  # Set attributes matching wrassp format
  attr(obj, "trackFormats") <- c("REAL32", "REAL32")
  attr(obj, "sampleRate") <- frame_rate  # Frames per second
  attr(obj, "origFreq") <- as.numeric(sample_rate)  # Original audio sample rate
  attr(obj, "startTime") <- 0.0
  attr(obj, "startRecord") <- 1L
  attr(obj, "endRecord") <- as.integer(n_frames)
  attr(obj, "fileInfo") <- c(20L, 2L)  # SSFF format
  class(obj) <- "AsspDataObj"

  return(obj)
}


# Set function attributes
attr(trk_yin, "ext") <- "yip"
attr(trk_yin, "tracks") <- c("F0", "prob")
attr(trk_yin, "outputType") <- "SSFF"
attr(trk_yin, "nativeFiletypes") <- c("wav")

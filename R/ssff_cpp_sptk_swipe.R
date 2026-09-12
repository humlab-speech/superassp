##' Track fundamental frequency using SWIPE (Sawtooth Waveform Inspired Pitch Estimator)
##'
##' Extracts F0 by matching the speech spectrum against sawtooth waveform templates
##' across candidate F0 values via SPTK. SWIPE is particularly effective for noisy
##' or challenging recording conditions where cross-correlation trackers (RAPT, Snack)
##' are less reliable. Its default voicing threshold (0.3) is lower than RAPT's (0.6),
##' yielding more voiced frames.
##'
##' @inheritParams trk_pitch_rapt
##' @param voicing_threshold Numeric. Voicing decision threshold (0–1). Default 0.3
##'   (more permissive than RAPT). Increase toward 0.5 to reduce false voiced frames.
##' @param parallel Logical. Use parallel processing for multiple files. \code{NULL}
##'   (default) enables automatically for 2+ files.
##' @param n_cores Integer. Number of cores for parallel processing. \code{NULL}
##'   (default) uses \code{detectCores() - 1}.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with track:
##'   \describe{
##'     \item{\code{f0}}{REAL32, fundamental frequency in Hz, n_frames × 1.
##'       Zero indicates unvoiced frames.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 100 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @usage trk_pitch_swipe(
##'   listOfFiles,
##'   beginTime = 0,
##'   endTime = 0,
##'   windowShift = 10,
##'   minF = 60,
##'   maxF = 400,
##'   voicing_threshold = 0.3,
##'   toFile = FALSE,
##'   explicitExt = "f0",
##'   outputDirectory = NULL,
##'   verbose = TRUE,
##'   parallel = NULL,
##'   n_cores = NULL
##' )
##' @param beginTime Start time for the extracted portion in seconds. Default: NULL (beginning of signal). Note: uses `beginTime`/`endTime` (seconds) matching DSP function conventions, unlike [read_audio()] which uses `begin`/`end`.
##' @param endTime The end time of the section of the sound files that should be analysed (in seconds). Use 0 for end of file.
##' @param windowShift Numeric. Frame shift in milliseconds; sets output frame rate (\code{1000 / windowShift} Hz). Default 5.0 ms (200 Hz). Must be strictly less than 32 ms (the 512-sample analysis window at 16 kHz). Values other than the training default (5 ms) may slightly reduce accuracy.
##' @param minF Numeric. Minimum F0 in Hz for the internal pitch estimator. Lower values allow lower-pitched voices but may increase false positives. Default 40.0 Hz.
##' @param maxF Numeric. Maximum F0 in Hz to treat as voiced. Default 400 Hz (speech). Must be <= 2093.75 Hz (model maximum; C7). For music, use 2093.75.
##' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the count written. If \code{FALSE}, return an \code{AsspDataObj} (single file only). Default \code{FALSE}.
##' @param explicitExt By default, a character "d" will be prepended to the file name suffix when writing the output to file. The user can also specify an explicit extension which will be used instead.
##' @param outputDirectory The directory where the slice file should be stored. If not defiled (NULL), the sparse slice file will placed in the same folder as the media file.
##' @param verbose Logical. Show a progress bar (sequential path) or a progress-aware parallel apply (`pbapply`/`pbmcapply`, if installed).
##' @export
##' @examples
##' wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
##'
##' \donttest{
##' f0 <- trk_pitch_swipe(wav, toFile = FALSE, verbose = FALSE)
##' head(as.data.frame(f0))
##'
##' # Custom range and voicing threshold
##' trk_pitch_swipe(wav, minF = 100, maxF = 500, voicing_threshold = 0.4,
##'                 toFile = FALSE, verbose = FALSE)
##' }
trk_pitch_swipe <- function(listOfFiles,
                  beginTime = 0.0,
                  endTime = 0.0,
                  windowShift = 10.0,
                  minF = 60.0,
                  maxF = 400.0,
                  voicing_threshold = 0.3,
                  toFile = FALSE,
                  explicitExt = "f0",
                  outputDirectory = NULL,
                  verbose = TRUE,
                  parallel = NULL,
                  n_cores = NULL) {

  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No input files specified in {.arg listOfFiles}")
  }

  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)

  validate_file_paths(listOfFiles, function_name = "trk_pitch_swipe")

  n_files <- length(listOfFiles)

  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime <- if (is.null(endTime)) 0.0 else endTime

  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime) == 1) endTime <- rep(endTime, n_files)

  makeOutputDirectory(outputDirectory, FALSE, "trk_pitch_swipe")

  if (verbose) format_apply_msg("trk_pitch_swipe", n_files, beginTime, endTime)

  process_single_file <- function(i) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[i]

    tryCatch({
      audio_obj <- assp_load_audio_for_dsp(
        file_path,
        begin = bt,
        end   = et,
        framework = "raw"
      )

      swipe_result <- swipe_cpp(
        audio_obj = audio_obj,
        minF = minF,
        maxF = maxF,
        windowShift = windowShift,
        voicing_threshold = voicing_threshold,
        verbose = FALSE
      )

      out_obj <- create_f0_asspobj(swipe_result, windowShift)

      if (toFile) {
        out_file <- generate_output_path(file_path, explicitExt, outputDirectory)
        write.AsspDataObj(out_obj, out_file)
        TRUE
      } else {
        out_obj
      }
    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(file_path)}}: {conditionMessage(e)}")
      if (toFile) FALSE else NULL
    })
  }

  results <- run_parallel_files(
    n_files = n_files,
    process_single_file = process_single_file,
    parallel = parallel,
    n_cores = n_cores,
    verbose = verbose,
    export_vars = c("listOfFiles", "beginTime", "endTime", "minF", "maxF",
                     "windowShift", "voicing_threshold", "toFile",
                     "explicitExt", "outputDirectory"),
    export_env = environment()
  )

  if (toFile) {
    n_success <- sum(unlist(results), na.rm = TRUE)
    if (verbose) {
      cli::cli_inform("Successfully processed {n_success} of {n_files} file{?s}")
    }
    return(invisible(n_success))
  } else {
    results <- results[!sapply(results, is.null)]
    if (length(results) == 1) {
      return(results[[1]])
    } else {
      return(results)
    }
  }
}

attr(trk_pitch_swipe, "ext") <- "f0"
attr(trk_pitch_swipe, "tracks") <- c("f0")
attr(trk_pitch_swipe, "outputType") <- "SSFF"
attr(trk_pitch_swipe, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
attr(trk_pitch_swipe, "suggestCaching") <- FALSE

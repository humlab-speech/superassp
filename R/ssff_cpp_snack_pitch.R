##' Snack Pitch Tracking (C++ implementation)
##'
##' @description Extract F0, voicing probability, RMS energy and
##'   autocorrelation peak using the Snack pitch tracker
##'   (normalized cross-correlation + dynamic programming, Talkin 1995).
##'
##'   This is a direct C++ wrapper around the original Snack/ESPS
##'   \code{dp_f0} algorithm.  Unlike \code{trk_rapt()} which only
##'   returns F0, this function exposes all four output tracks.
##'
##'   All input media formats are supported via the av package.
##'
##' @param listOfFiles Vector of file paths to process
##' @param beginTime Start time in seconds (default: 0.0)
##' @param endTime End time in seconds (default: 0.0 = end of file)
##' @param windowShift Frame shift in milliseconds (default: 10.0)
##' @param minF Minimum F0 in Hz (default: 50.0)
##' @param maxF Maximum F0 in Hz (default: 550.0)
##' @param voiceBias Bias toward voiced hypothesis (default: 0.0,
##'   range approx -0.5 to 0.5)
##' @param toFile Write results to file (default: TRUE)
##' @param explicitExt Output file extension (default: "snackpitch")
##' @param outputDirectory Output directory (default: NULL = same as input)
##' @param verbose Show progress messages (default: TRUE)
##'
##' @return If toFile=TRUE, returns the number of successfully processed files.
##'   If toFile=FALSE, returns AsspDataObj or list of AsspDataObj objects
##'   with tracks: f0, voicing, rms, acpeak.
##'
##' @export
##' @examples
##' \dontrun{
##' # Full 4-track pitch analysis
##' res <- trk_snackp("recording.wav", toFile = FALSE)
##' names(res)  # "f0" "voicing" "rms" "acpeak"
##'
##' # Custom F0 range
##' trk_snackp("speech.mp3", minF = 75, maxF = 300)
##' }
trk_snackp <- function(listOfFiles,
                       beginTime = 0.0,
                       endTime = 0.0,
                       windowShift = 10.0,
                       minF = 50.0,
                       maxF = 550.0,
                       voiceBias = 0.0,
                       toFile = TRUE,
                       explicitExt = "snackpitch",
                       outputDirectory = NULL,
                       verbose = TRUE) {

  if (is.null(listOfFiles) || length(listOfFiles) == 0)
    cli::cli_abort("No input files specified in {.arg listOfFiles}")

  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)

  files_exist <- file.exists(listOfFiles)
  if (!all(files_exist)) {
    missing <- listOfFiles[!files_exist]
    cli::cli_abort(c("!" = "Files not found:", "x" = "{.file {fast_basename(missing)}}"))
  }

  n_files <- length(listOfFiles)
  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime   <- if (is.null(endTime))   0.0 else endTime
  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime)   == 1) endTime   <- rep(endTime,   n_files)

  makeOutputDirectory(outputDirectory, FALSE, "trk_snackp")
  if (verbose) format_apply_msg("trk_snackp", n_files, beginTime, endTime)

  results <- vector("list", n_files)

  if (verbose && n_files > 1)
    cli::cli_progress_bar("Processing files", total = n_files,
      format = "{cli::pb_spin} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}")

  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]; et <- endTime[i]

    tryCatch({
      audio_obj <- read_audio(file_path, begin = bt, end = et)

      res <- snackp_cpp(
        audio_obj   = audio_obj,
        minF        = minF,
        maxF        = maxF,
        windowShift = windowShift,
        voiceBias   = voiceBias,
        verbose     = FALSE
      )

      out_obj <- create_snackp_asspobj(res, windowShift)

      if (toFile) {
        out_file <- generate_output_path(file_path, explicitExt, outputDirectory)
        write.AsspDataObj(out_obj, out_file)
        results[[i]] <- TRUE
      } else {
        results[[i]] <- out_obj
      }
    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(file_path)}}: {conditionMessage(e)}")
      results[[i]] <<- if (toFile) FALSE else NULL
    })

    if (verbose && n_files > 1) cli::cli_progress_update()
  }

  if (verbose && n_files > 1) cli::cli_progress_done()

  if (toFile) {
    n_ok <- sum(unlist(results), na.rm = TRUE)
    if (verbose) cli::cli_inform("Processed {n_ok} of {n_files} file{?s}")
    return(invisible(n_ok))
  } else {
    results <- results[!sapply(results, is.null)]
    if (length(results) == 1) results[[1]] else results
  }
}

attr(trk_snackp, "ext")             <- "snackpitch"
attr(trk_snackp, "tracks")          <- c("f0", "voicing", "rms", "acpeak")
attr(trk_snackp, "outputType")      <- "SSFF"
attr(trk_snackp, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")

##' Snack Formant Tracking (C++ implementation)
##'
##' @description Extract formant frequencies and bandwidths using the
##'   Snack/ESPS LPC + dynamic programming formant tracker
##'   (Talkin / AT&T / KTH).
##'
##'   This is a direct C++ adaptation of the original Snack Sound Toolkit
##'   formant tracking algorithm (\code{jkFormant.c} + \code{sigproc2.c}).
##'
##'   All input media formats are supported via the av package.
##'
##' @param listOfFiles Vector of file paths to process
##' @param beginTime Start time in seconds (default: 0.0)
##' @param endTime End time in seconds (default: 0.0 = end of file)
##' @param numFormants Number of formants to track (default: 4, max: 7)
##' @param lpcOrder LPC order (default: 12)
##' @param windowLength Analysis window duration in seconds (default: 0.049)
##' @param windowShift Frame shift in milliseconds (default: 10.0)
##' @param preEmphasis Pre-emphasis factor (default: 0.7)
##' @param dsFreq Downsample target frequency in Hz (default: 10000)
##' @param nomF1 Nominal F1 for DP cost (default: -10.0 = use built-in defaults)
##' @param lpcType LPC method: 0=autocorrelation, 1=stabilized covariance,
##'   2=covariance (default: 0)
##' @param windowType Window type: 0=rectangular, 1=Hamming, 2=cos^4,
##'   3=Hanning (default: 2)
##' @param toFile Write results to file (default: TRUE)
##' @param explicitExt Output file extension (default: "snackfmt")
##' @param outputDirectory Output directory (default: NULL = same as input)
##' @param verbose Show progress messages (default: TRUE)
##'
##' @return If toFile=TRUE, returns the number of successfully processed files.
##'   If toFile=FALSE, returns AsspDataObj or list of AsspDataObj objects
##'   with tracks: fm (formant frequencies), bw (formant bandwidths).
##'
##' @export
##' @examples
##' \dontrun{
##' # 4-formant tracking
##' res <- trk_snackf("recording.wav", toFile = FALSE)
##' names(res)  # "fm" "bw"
##'
##' # Custom LPC order and formant count
##' trk_snackf("speech.mp3", numFormants = 5, lpcOrder = 14)
##' }
trk_snackf <- function(listOfFiles,
                       beginTime = 0.0,
                       endTime = 0.0,
                       numFormants = 4,
                       lpcOrder = 12,
                       windowLength = 0.049,
                       windowShift = 10.0,
                       preEmphasis = 0.7,
                       dsFreq = 10000.0,
                       nomF1 = -10.0,
                       lpcType = 0,
                       windowType = 2,
                       toFile = TRUE,
                       explicitExt = "snackfmt",
                       outputDirectory = NULL,
                       verbose = TRUE) {

  if (is.null(listOfFiles) || length(listOfFiles) == 0)
    cli::cli_abort("No input files specified in {.arg listOfFiles}")

  if (numFormants < 1 || numFormants > 7)
    cli::cli_abort("{.arg numFormants} must be between 1 and 7")

  if (numFormants > (lpcOrder - 4) / 2)
    cli::cli_abort("{.arg lpcOrder} ({lpcOrder}) too low for {.arg numFormants} ({numFormants}) formants; need at least {numFormants * 2 + 4}")

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

  makeOutputDirectory(outputDirectory, FALSE, "trk_snackf")
  if (verbose) format_apply_msg("trk_snackf", n_files, beginTime, endTime)

  results <- vector("list", n_files)

  if (verbose && n_files > 1)
    cli::cli_progress_bar("Processing files", total = n_files,
      format = "{cli::pb_spin} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}")

  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]; et <- endTime[i]

    tryCatch({
      audio_obj <- read_audio(file_path, begin = bt, end = et)

      res <- snackf_cpp(
        audio_obj    = audio_obj,
        numFormants  = numFormants,
        lpcOrder     = lpcOrder,
        windowLength = windowLength,
        windowShift  = windowShift,
        preEmphasis  = preEmphasis,
        dsFreq       = dsFreq,
        nomF1        = nomF1,
        lpcType      = lpcType,
        windowType   = windowType,
        verbose      = FALSE
      )

      out_obj <- create_formant_asspobj(res, windowShift, numFormants)

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

attr(trk_snackf, "ext")             <- "snackfmt"
attr(trk_snackf, "tracks")          <- c("fm", "bw")
attr(trk_snackf, "outputType")      <- "SSFF"
attr(trk_snackf, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")

##' Track formants and bandwidths using the Snack/ESPS LPC tracker
##'
##' Extracts formant frequencies and bandwidths using the Snack Sound Toolkit
##' dynamic-programming LPC formant tracker (Talkin / AT&T / KTH,
##' \code{jkFormant.c}). This tracker applies normalized cross-correlation
##' for voicing detection and a DP smoother for temporal continuity. It is a
##' strong general-purpose alternative to Burg-method trackers, particularly
##' for lower sample-rate or telephone-bandwidth speech.
##'
##' @param listOfFiles Character vector of audio file paths. Any format supported by
##'   \pkg{av} is accepted; non-native inputs are transcoded automatically.
##' @param beginTime Numeric. Start of analysis window in seconds. Default 0 (file start).
##' @param endTime Numeric. End of analysis window in seconds. Default 0 (file end).
##' @param numFormants Integer. Number of formants to track (1–7). Default 4.
##' @param lpcOrder Integer. LPC order; must satisfy \code{lpcOrder >= numFormants * 2 + 4}.
##'   Default 12.
##' @param windowLength Numeric. Analysis window duration in seconds. Default 0.049 s.
##' @param windowShift Numeric. Frame shift in milliseconds; sets output frame rate
##'   (1000 / windowShift Hz). Default 10.0 ms.
##' @param preEmphasis Numeric. Pre-emphasis filter coefficient (0–1). Default 0.7.
##' @param dsFreq Numeric. Target downsample frequency in Hz; limits formant search
##'   to 0–dsFreq/2 Hz. Default 10000 Hz.
##' @param nomF1 Numeric. Nominal F1 for DP cost function. \code{-10} (default) uses
##'   built-in defaults derived from \code{dsFreq}.
##' @param lpcType Integer. LPC method: 0 = autocorrelation, 1 = stabilized covariance,
##'   2 = covariance. Default 0.
##' @param windowType Integer. Window type: 0 = rectangular, 1 = Hamming, 2 = cos^4,
##'   3 = Hanning. Default 2.
##' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
##'   count written invisibly. If \code{FALSE}, return an \code{AsspDataObj}.
##'   Default \code{TRUE}.
##' @param explicitExt Character. Output file extension. Default \code{"snackfmt"}.
##' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
##'   writes alongside the input file.
##' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
##'   \describe{
##'     \item{\code{fm}}{REAL32, formant frequencies in Hz, n_frames × numFormants.
##'       Columns correspond to F1, F2, …, F\{numFormants\}.}
##'     \item{\code{bw}}{REAL32, formant bandwidths in Hz, n_frames × numFormants.}
##'   }
##'   Frame rate: \code{1000 / windowShift} Hz (default 100 Hz).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @export
##' @examples
##' \dontrun{
##' # 4-formant tracking
##' res <- trk_formant_snack("recording.wav", toFile = FALSE)
##' names(res)  # "fm" "bw"
##'
##' # Custom LPC order and formant count
##' trk_formant_snack("speech.mp3", numFormants = 5, lpcOrder = 14)
##' }
trk_formant_snack <- function(listOfFiles,
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

  makeOutputDirectory(outputDirectory, FALSE, "trk_formant_snack")
  if (verbose) format_apply_msg("trk_formant_snack", n_files, beginTime, endTime)

  results <- vector("list", n_files)

  if (verbose && n_files > 1)
    cli::cli_progress_bar("Processing files", total = n_files,
      format = "{cli::pb_spin} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}")

  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]; et <- endTime[i]

    results[[i]] <- tryCatch({
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
        TRUE
      } else {
        out_obj
      }
    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(file_path)}}: {conditionMessage(e)}")
      if (toFile) FALSE else NULL
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

attr(trk_formant_snack, "ext")             <- "snackfmt"
attr(trk_formant_snack, "tracks")          <- c("fm", "bw")
attr(trk_formant_snack, "outputType")      <- "SSFF"
attr(trk_formant_snack, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")

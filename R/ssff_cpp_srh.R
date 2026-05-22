##' Track fundamental frequency using the Summation of Residual Harmonics (SRH)
##'
##' Extracts F0 and a voiced/unvoiced decision using SRH (Drugman & Alwan 2011),
##' a two-pass harmonic-summation pitch estimator operating on the LPC residual.
##' SRH is robust in noisy conditions and produces an integrated VAD decision. Audio
##' is resampled to 16 kHz internally. The fixed 10 ms hop differs from RAPT/SWIPE,
##' which honour the \code{windowShift} parameter.
##'
##' @inheritParams trk_pitch_rapt
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
##'   \describe{
##'     \item{\code{f0}}{REAL32, fundamental frequency in Hz, n_frames × 1.
##'       Zero indicates unvoiced frames.}
##'     \item{\code{vad}}{REAL32, voiced/unvoiced decision (0 = unvoiced, 1 = voiced),
##'       n_frames × 1.}
##'   }
##'   Frame rate: 100 Hz (fixed 10 ms hop).
##'   If \code{toFile = TRUE}: integer count of files written, returned invisibly.
##'
##' @export
##' @examples
##' \dontrun{
##' # Extract F0 using SRH
##' trk_pitch_srh("recording.wav")
##'
##' # Process with custom F0 range
##' trk_pitch_srh("speech.wav", minF = 80, maxF = 300)
##' }
##'
##' @references
##' \insertCite{Drugman2011SRH}{superassp}
trk_pitch_srh <- function(listOfFiles,
                    beginTime = 0.0,
                    endTime = 0.0,
                    minF = 50,
                    maxF = 500,
                    toFile = TRUE,
                    explicitExt = "srh",
                    outputDirectory = NULL,
                    verbose = TRUE) {

  if (is.null(listOfFiles) || length(listOfFiles) == 0) {
    cli::cli_abort("No input files specified in {.arg listOfFiles}")
  }

  listOfFiles <- fast_strip_file_protocol(listOfFiles)
  listOfFiles <- normalizePath(path.expand(listOfFiles), mustWork = FALSE)

  files_exist <- file.exists(listOfFiles)
  if (!all(files_exist)) {
    missing_files <- listOfFiles[!files_exist]
    cli::cli_abort(c(
      "!" = "Some files do not exist:",
      "x" = "{.file {fast_basename(missing_files)}}"
    ))
  }

  n_files <- length(listOfFiles)

  beginTime <- if (is.null(beginTime)) 0.0 else beginTime
  endTime <- if (is.null(endTime)) 0.0 else endTime

  if (length(beginTime) == 1) beginTime <- rep(beginTime, n_files)
  if (length(endTime) == 1) endTime <- rep(endTime, n_files)

  makeOutputDirectory(outputDirectory, FALSE, "trk_pitch_srh")

  if (verbose) format_apply_msg("trk_pitch_srh", n_files, beginTime, endTime)

  results <- vector("list", n_files)

  if (verbose && n_files > 1) {
    cli::cli_progress_bar(
      "Processing files",
      total = n_files,
      format = "{cli::pb_spin} {cli::pb_current}/{cli::pb_total} | ETA: {cli::pb_eta}"
    )
  }

  # SRH works at 16 kHz internally
  target_sr <- 16000L

  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]
    et <- endTime[i]

    results[[i]] <- tryCatch({
      # Load audio via av, resample to 16 kHz
      invisible(utils::capture.output(
        audio_data <- av::read_audio_bin(
          audio = file_path,
          start_time = if (bt > 0) bt else NULL,
          end_time = if (et > 0) et else NULL,
          channels = 1,
          sample_rate = target_sr
        ),
        type = "message"
      ))

      audio_vec <- as.numeric(audio_data)
      edge <- as.integer(round(c(minF, maxF)))
      srh_result <- srh_variant_cpp(audio_vec, target_sr, edge)

      n_frames <- length(srh_result$f0)
      if (n_frames == 0) {
        cli::cli_warn("SRH returned empty result for {.file {basename(file_path)}}")
        if (toFile) FALSE else NULL
      } else {
        out_obj <- list(
          f0  = matrix(srh_result$f0, ncol = 1),
          vad = matrix(as.numeric(srh_result$vad), ncol = 1)
        )
        attr(out_obj, "trackFormats") <- c("REAL32", "REAL32")
        attr(out_obj, "sampleRate")   <- 100
        attr(out_obj, "origFreq")     <- as.numeric(target_sr)
        attr(out_obj, "startTime")    <- 0.0
        attr(out_obj, "startRecord")  <- 1L
        attr(out_obj, "endRecord")    <- as.integer(n_frames)
        attr(out_obj, "fileInfo")     <- c(20L, 2L)
        class(out_obj) <- "AsspDataObj"

        if (toFile) {
          out_file <- generate_output_path(file_path, explicitExt, outputDirectory)
          write.AsspDataObj(out_obj, out_file)
          TRUE
        } else {
          out_obj
        }
      }
    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(file_path)}}: {conditionMessage(e)}")
      if (toFile) FALSE else NULL
    })

    if (verbose && n_files > 1) {
      cli::cli_progress_update()
    }
  }

  if (verbose && n_files > 1) {
    cli::cli_progress_done()
  }

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

attr(trk_pitch_srh, "ext") <- "srh"
attr(trk_pitch_srh, "tracks") <- c("f0", "vad")
attr(trk_pitch_srh, "outputType") <- "SSFF"
attr(trk_pitch_srh, "nativeFiletypes") <- c("wav", "flac", "mp3", "mp4", "mkv", "avi")
attr(trk_pitch_srh, "suggestCaching") <- FALSE

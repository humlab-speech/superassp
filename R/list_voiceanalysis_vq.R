##' Per-GCI voice-quality summary via voiceanalysis
##'
##' Extracts the canonical voice-source parameters
##' \insertCite{KaneGobl2013}{superassp} from an audio file using
##' \pkg{voiceanalysis}: NAQ (Normalised Amplitude Quotient), QOQ
##' (Quasi-Open Quotient), H1H2 spectral tilt, and HRF (Harmonic Richness
##' Factor). Internally runs SE-VQ for GCIs, IAIF for the glottal flow
##' derivative, then \code{vat_voice_quality()}.
##'
##' Returns one row per GCI. To get continuous-track equivalents, interp
##' to a 10 ms grid via \code{stats::approx}.
##'
##' @param listOfFiles Character vector of audio file paths.
##' @param beginTime,endTime Analysis window (seconds). Defaults: full file.
##' @param toFile Logical. If TRUE, persist results as JSTF and return paths.
##'   Default FALSE.
##' @param explicitExt File extension. Default "vqv".
##' @param outputDirectory Output directory. Default NULL = next to input.
##' @param verbose Logical. Default TRUE.
##'
##' @return If \code{toFile = FALSE} and \code{length(listOfFiles) == 1}:
##'   a named list with vectors \code{gci_time} (seconds),
##'   \code{NAQ}, \code{QOQ}, \code{H1H2}, \code{HRF}. For multiple files,
##'   a list of such lists. If \code{toFile = TRUE}: invisible vector of
##'   output paths.
##'
##' @references
##' \insertCite{KaneGobl2013}{superassp}
##' @seealso \code{\link{lst_covarep_vq}}
##' @export
lst_vq_vat <- function(listOfFiles,
                       beginTime = 0.0,
                       endTime = 0.0,
                       toFile = FALSE,
                       explicitExt = "vqv",
                       outputDirectory = NULL,
                       verbose = TRUE) {

  if (!requireNamespace("voiceanalysis", quietly = TRUE))
    cli::cli_abort(c("Package {.pkg voiceanalysis} is required.",
                     "i" = "Install via {.code pak::pkg_install('jckane/Voice_Analysis_Toolkit/voiceanalysis')}"))
  if (is.null(listOfFiles) || length(listOfFiles) == 0)
    cli::cli_abort("No input files specified in {.arg listOfFiles}")
  listOfFiles <- normalizePath(path.expand(fast_strip_file_protocol(listOfFiles)), mustWork = FALSE)
  if (!all(file.exists(listOfFiles)))
    cli::cli_abort("Some files do not exist.")

  n_files <- length(listOfFiles)
  beginTime <- if (length(beginTime) == 1) rep(beginTime, n_files) else beginTime
  endTime   <- if (length(endTime)   == 1) rep(endTime,   n_files) else endTime
  makeOutputDirectory(outputDirectory, FALSE, "lst_vq_vat")
  if (verbose) format_apply_msg("lst_vq_vat", n_files, beginTime, endTime)

  results <- vector("list", n_files)
  for (i in seq_len(n_files)) {
    file_path <- listOfFiles[i]
    bt <- beginTime[i]; et <- endTime[i]
    tryCatch({
      invisible(utils::capture.output(
        audio_data <- av::read_audio_bin(
          audio = file_path,
          start_time = if (bt > 0) bt else NULL,
          end_time   = if (et > 0) et else NULL,
          channels = 1
        ), type = "message"))
      fs <- attr(audio_data, "sample_rate")
      wave <- as.numeric(audio_data)
      mx <- max(abs(wave)); if (mx > 1) wave <- wave / mx

      se  <- voiceanalysis::vat_se_vq(wave, fs)
      if (length(se$GCI) < 3) {
        cli::cli_warn("Too few GCIs for VQ on {.file {basename(file_path)}}")
        results[[i]] <- if (toFile) NA_character_ else NULL; next
      }
      iaif <- voiceanalysis::vat_iaif(wave, fs, GCI = se$GCI)
      vq   <- voiceanalysis::vat_voice_quality(iaif$dg, fs, se$GCI)

      out <- list(
        gci_time = (se$GCI - 1) / fs,
        NAQ  = as.numeric(vq$NAQ),
        QOQ  = as.numeric(vq$QOQ),
        H1H2 = as.numeric(vq$H1H2),
        HRF  = as.numeric(vq$HRF),
        meta = list(F0mean = se$F0mean, fs = fs,
                     file = basename(file_path),
                     generator = "lst_vq_vat",
                     package = "voiceanalysis")
      )

      if (toFile) {
        out_file <- generate_output_path(file_path, explicitExt, outputDirectory)
        # JSTF-style write via jsonlite (lightweight; downstream tools should
        # be able to read this as JSON or via jsonlite::read_json)
        jsonlite::write_json(out, out_file, auto_unbox = TRUE, digits = 6)
        results[[i]] <- out_file
      } else {
        results[[i]] <- out
      }
    }, error = function(e) {
      warning("Error processing ", basename(file_path), ": ", e$message, call. = FALSE)
      results[i] <- list(NULL)
    })
  }

  if (toFile) invisible(unlist(results))
  else if (n_files == 1) results[[1]] else results
}

attr(lst_vq_vat, "ext")              <- "vqv"
attr(lst_vq_vat, "outputType")       <- "JSTF"
attr(lst_vq_vat, "nativeFiletypes")  <- c("wav", "flac", "mp3", "mp4")
attr(lst_vq_vat, "suggestCaching")   <- FALSE

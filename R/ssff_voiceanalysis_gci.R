##' Detect glottal closure instants (GCIs) using SE-VQ via voiceanalysis
##'
##' Detects GCIs using the improved SEDREAMS / SE-VQ algorithm of Kane &
##' Gobl (2013) \insertCite{KaneGobl2013}{superassp} as implemented by the
##' \pkg{voiceanalysis} package (bit-faithful Rcpp port of the MATLAB Voice
##' Analysis Toolkit). Supports both fixed-F0 and the variable-F0 variant
##' optimised for highly expressive speech.
##'
##' @inheritParams trk_pitch_rapt
##' @param var_f0 Logical. If \code{TRUE}, use the variable-F0 SE-VQ variant
##'   (recommended for expressive speech). Default \code{FALSE}.
##' @param f0_min Minimum F0 in Hz (default 20).
##' @param f0_max Maximum F0 in Hz (default 500).
##' @param use_creak Logical. If \code{TRUE}, run \pkg{voiceanalysis}'s creak
##'   detector and feed its decisions into the SE-VQ creaky post-processing
##'   step. Default \code{FALSE}.
##'
##' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
##'   \describe{
##'     \item{\code{gci_sample}}{INT32, GCI sample indices (1-based) in the
##'       resampled signal, n_gci × 1.}
##'     \item{\code{residual}}{REAL32, LP residual signal (normalised),
##'       n_samples × 1.}
##'   }
##'   The \code{gci_sample} track is a sparse event list; the \code{residual}
##'   track is per-sample at the audio rate.
##'   If \code{toFile = TRUE}: invisibly returns the count of files written.
##'
##' @references
##' \insertCite{KaneGobl2013}{superassp}
##' @seealso \code{\link{trk_covarep_vq_gci}}
##' @export
trk_gci_vat <- function(listOfFiles,
                        beginTime = 0.0,
                        endTime = 0.0,
                        var_f0 = FALSE,
                        f0_min = 20,
                        f0_max = 500,
                        use_creak = FALSE,
                        toFile = TRUE,
                        explicitExt = "gciv",
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

  makeOutputDirectory(outputDirectory, FALSE, "trk_gci_vat")
  if (verbose) format_apply_msg("trk_gci_vat", n_files, beginTime, endTime)

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

      creak_arg <- NULL
      if (use_creak) {
        cd <- voiceanalysis::vat_creak_detect(wave, fs)
        # Up-sample decision vector to per-sample
        n_samples <- length(wave); hop <- max(1, round(0.010 * fs))
        creak_per_sample <- integer(n_samples)
        for (k in seq_along(cd$decision)) {
          a <- (k - 1) * hop + 1; b <- min(n_samples, k * hop)
          creak_per_sample[a:b] <- cd$decision[k]
        }
        creak_arg <- creak_per_sample
      }

      vat_res <- voiceanalysis::vat_se_vq(
        wave, fs, f0 = NULL, VUV = NULL,
        creak = creak_arg, var_f0 = var_f0
      )

      n_gci    <- length(vat_res$GCI)
      n_samples <- length(vat_res$res)
      if (n_gci == 0) {
        cli::cli_warn("SE-VQ returned no GCIs for {.file {basename(file_path)}}")
        results[[i]] <- if (toFile) FALSE else NULL; next
      }

      out_obj <- list(
        gci_sample = matrix(as.integer(vat_res$GCI), ncol = 1),
        residual   = matrix(as.numeric(vat_res$res), ncol = 1)
      )
      # Two heterogeneous tracks at different rates is not strictly SSFF —
      # we stash GCI as a separate attribute and expose 'residual' as the
      # canonical per-sample track. Downstream code should read $gci_sample
      # via attributes when needed.
      attr(out_obj, "trackFormats") <- c("INT32", "REAL32")
      attr(out_obj, "sampleRate")   <- as.numeric(fs)
      attr(out_obj, "origFreq")     <- as.numeric(fs)
      attr(out_obj, "startTime")    <- as.numeric(bt)
      attr(out_obj, "startRecord")  <- 1L
      attr(out_obj, "endRecord")    <- as.integer(n_samples)
      attr(out_obj, "fileInfo")     <- c(20L, 2L)
      attr(out_obj, "n_gci")        <- n_gci
      attr(out_obj, "F0mean")       <- vat_res$F0mean
      class(out_obj) <- "AsspDataObj"

      if (toFile) {
        out_file <- generate_output_path(file_path, explicitExt, outputDirectory)
        write.AsspDataObj(out_obj, out_file)
        results[[i]] <- TRUE
      } else {
        results[[i]] <- out_obj
      }
    }, error = function(e) {
      warning("Error processing ", basename(file_path), ": ", e$message, call. = FALSE)
      results[i] <- list(NULL)
    })
  }

  if (toFile) invisible(sum(unlist(results), na.rm = TRUE))
  else if (n_files == 1) results[[1]] else results
}

attr(trk_gci_vat, "ext")              <- "gciv"
attr(trk_gci_vat, "tracks")           <- c("gci_sample", "residual")
attr(trk_gci_vat, "outputType")       <- "SSFF"
attr(trk_gci_vat, "nativeFiletypes")  <- c("wav", "flac", "mp3", "mp4")
attr(trk_gci_vat, "suggestCaching")   <- FALSE

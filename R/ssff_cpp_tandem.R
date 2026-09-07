#' Track pitch and voiced speech using the TANDEM-STRAIGHT algorithm
#'
#' Estimates F0 and per-frame voicing probability using a gammatone filterbank
#' combined with neural network-based pitch tracking (Hu & Wang 2010), which
#' simultaneously segregates voiced speech from noise. TANDEM is robust to noise
#' and reverberation and can track multiple simultaneous pitch sources.
#'
#' @note The core processing is currently a placeholder; full TANDEM C++ integration
#'   is under development. Results reflect the algorithm framework but may not match
#'   the published TANDEM-STRAIGHT output.
#'
#' @param listOfFiles Character vector of audio file paths. Any format supported by
#'   \pkg{av} is accepted; audio is resampled to \code{target_sample_rate} Hz internally.
#' @param minF Numeric. Minimum F0 in Hz. Default 50 Hz.
#' @param maxF Numeric. Maximum F0 in Hz. Default 500 Hz.
#' @param target_sample_rate Numeric. Internal processing sample rate in Hz.
#'   TANDEM requires 20000 Hz. Default 20000.
#' @param return_mask Logical. Return time-frequency voiced mask (currently unused).
#'   Default \code{FALSE}.
#' @param toFile Logical. If \code{TRUE}, write SSFF output files and return the
#'   paths written. If \code{FALSE}, return an \code{AsspDataObj}.
#'   Default \code{FALSE}.
#' @param explicitExt Character. Output file extension. Default \code{"tnd"}.
#' @param outputDirectory Character. Directory for output files. \code{NULL} (default)
#'   writes alongside the input file.
#' @param verbose Logical. Print per-file progress. Default \code{TRUE}.
#' @param ... Additional arguments (currently unused).
#'
#' @return If \code{toFile = FALSE}: an \code{AsspDataObj} with tracks:
#'   \describe{
#'     \item{\code{pitch}}{REAL64, fundamental frequency in Hz, n_frames × 1.
#'       Zero indicates unvoiced frames.}
#'     \item{\code{voicing_prob}}{REAL64, voicing probability, 0–1, n_frames × 1.}
#'   }
#'   Frame rate: 100 Hz (fixed 10 ms hop).
#'   If \code{toFile = TRUE}: character vector of output file paths.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic pitch tracking
#' result <- trk_tandem("speech.wav")
#' plot(result$pitch, type = "l", main = "TANDEM Pitch Track")
#'
#' # With noisy speech
#' result <- trk_tandem("noisy_speech.wav", minF = 80, maxF = 400)
#'
#' # Batch processing
#' files <- c("speaker1.wav", "speaker2.wav", "speaker3.wav")
#' results <- trk_tandem(files, verbose = TRUE)
#'
#' # Save to files
#' trk_tandem("speech.wav", toFile = TRUE, outputDirectory = "output/")
#' }
#'
#' @references
#' \insertRef{hu2010tandem}{superassp}
#'
#' \insertRef{hu2011unvoiced}{superassp}
#'
#' @seealso \code{\link{trk_pitch_rapt}}, \code{\link{trk_pitch_swipe}}, \code{\link{trk_pitch_yin}}
#'   for other pitch tracking methods
trk_tandem <- function(
  listOfFiles,
  minF = 50,
  maxF = 500,
  target_sample_rate = 20000,
  return_mask = FALSE,
  toFile = FALSE,
  explicitExt = "tnd",
  outputDirectory = NULL,
  verbose = TRUE,
  ...
) {
  # Validate inputs
  if (!is.character(listOfFiles)) {
    cli::cli_abort("listOfFiles must be a character vector")
  }
  
  # Check file existence
  missing_files <- listOfFiles[!file.exists(listOfFiles)]
  if (length(missing_files) > 0) {
    cli::cli_abort(c("{length(missing_files)} file{?s} not found:",
                     stats::setNames(missing_files, rep("*", length(missing_files)))))
  }
  
  n_files <- length(listOfFiles)
  .warn_if_lossy_input(listOfFiles)
  results <- vector("list", n_files)
  
  if (verbose) format_apply_msg("trk_tandem", n_files)
  
  for (i in seq_along(listOfFiles)) {
    results[[i]] <- tryCatch({
      if (verbose && n_files > 1) {
        message("  [", i, "/", n_files, "] ", basename(listOfFiles[i]))
      }

      # Load audio via av package
      tryCatch({
        invisible(utils::capture.output(
          audio_data <- av::read_audio_bin(
            listOfFiles[i],
            channels = 1  # TANDEM requires mono
          ),
          type = "message"
        ))
      }, error = function(e) {
        cli::cli_abort("Failed to load audio file {.file {basename(listOfFiles[i])}}: {e$message}")
      })
    
      orig_sr <- attr(audio_data, "sample_rate")
      audio_vec <- as.numeric(audio_data)

      # Resample to 20 kHz if needed (using av package)
      if (orig_sr != target_sample_rate) {
        if (verbose) {
          cli::cli_inform("Resampling {basename(listOfFiles[i])} {orig_sr} -> {target_sample_rate} Hz")
        }
        temp_wav <- tempfile(fileext = ".wav")
        on.exit(unlink(temp_wav), add = TRUE)

        tryCatch({
          invisible(utils::capture.output(
            invisible(utils::capture.output(
              av::av_audio_convert(
                listOfFiles[i],
                temp_wav,
                format = "wav",
                sample_rate = target_sample_rate,
                channels = 1
              ),
              type = "message"
            ))
          ))

          invisible(utils::capture.output(
            audio_data <- av::read_audio_bin(temp_wav, channels = 1),
            type = "message"
          ))
          audio_vec <- as.numeric(audio_data)
        }, error = function(e) {
          cli::cli_abort("Resampling failed for {.file {basename(listOfFiles[i])}}: {e$message}")
        })
      }
    
    # The vendored voicedMask constructor (src/tandem/tandem_64/, submodule --
    # not ours to patch) hard-codes reading its network weights from
    # "net/MLP*.64.dat" relative to the working directory; the net_path
    # argument passed to tandem_pitch_cpp() below is not consulted by that
    # constructor. Stage symlinks there for this call and remove exactly what
    # we added on exit -- tracking per-file (not just "did the dir exist
    # before") so a leftover net/ from an earlier interrupted run doesn't
    # suppress cleanup on every subsequent call.
    net_dir <- file.path(getwd(), "net")
    net_dir_created <- !dir.exists(net_dir)
    if (net_dir_created) dir.create(net_dir)

    net_source <- system.file("tandem_net", package = "superassp")
    net_files_created <- character(0)
    for (net_file in c("MLP1.64.dat", "MLP2.64.dat", "MLP3.64.dat")) {
      src <- file.path(net_source, net_file)
      dst <- file.path(net_dir, net_file)
      if (!file.exists(dst) && file.exists(src)) {
        ok <- tryCatch(file.symlink(src, dst), error = function(e) FALSE)
        if (!isTRUE(ok)) ok <- file.copy(src, dst)
        if (isTRUE(ok)) net_files_created <- c(net_files_created, dst)
      }
    }

    on.exit({
      unlink(net_files_created)
      if (net_dir_created && dir.exists(net_dir) &&
          length(list.files(net_dir, all.files = TRUE, no.. = TRUE)) == 0) {
        unlink(net_dir, recursive = TRUE)
      }
    }, add = TRUE)

    # Call TANDEM C++ wrapper (suppress C-level stdout/stderr)
    invisible(utils::capture.output(
      invisible(utils::capture.output(
        tandem_result <- tandem_pitch_cpp(
          audio_signal = audio_vec,
          sample_rate = target_sample_rate,
          min_pitch = minF,
          max_pitch = maxF,
          net_path = system.file("tandem_net", package = "superassp")
        ),
        type = "message"
      ))
    ))
    
    # Check status
    if ("status" %in% names(tandem_result) && tandem_result$status == "placeholder") {
      if (verbose && i == 1) {
        message("    Status: Using placeholder implementation (full TANDEM integration pending)")
      }
    }
    
    # Convert to AsspDataObj (tracks must be single-column matrices)
    assp_obj <- list(
      pitch = matrix(tandem_result$pitch, ncol = 1),
      voicing_prob = matrix(tandem_result$voicing_prob, ncol = 1)
    )
    
      # Validate result
      if (length(tandem_result$pitch) == 0) {
        cli::cli_abort("TANDEM returned empty pitch track for {.file {basename(listOfFiles[i])}}.")
      }

    # Set attributes
    attr(assp_obj, "sampleRate") <- 100  # Analysis rate (100 Hz frames)
    attr(assp_obj, "startTime") <- 0
    attr(assp_obj, "startRecord") <- 1L
    attr(assp_obj, "endRecord") <- as.integer(length(tandem_result$pitch))
    attr(assp_obj, "trackFormats") <- c("REAL64", "REAL64")
    attr(assp_obj, "fileInfo") <- as.integer(c(20L, 2L))  # SSFF format
    class(assp_obj) <- "AsspDataObj"
    
    # Write to file if requested
    if (toFile) {
      base_name <- tools::file_path_sans_ext(basename(listOfFiles[i]))
      out_dir <- if (is.null(outputDirectory)) {
        dirname(listOfFiles[i])
      } else {
        outputDirectory
      }
      
      # Create output directory if needed
      if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE)
      }
      
      output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))

      tryCatch(
        write.AsspDataObj(assp_obj, output_path),
        error = function(e) {
          cli::cli_abort("Failed to write output file {.file {output_path}}: {e$message}")
        }
      )
      output_path

    } else {
      assp_obj
    }
    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(listOfFiles[i])}}: {conditionMessage(e)}")
      if (toFile) FALSE else NULL
    })
  }
  
  # Simplify output for single file
  if (n_files == 1) {
    return(results[[1]])
  } else {
    return(results)
  }
}

# Set function attributes
attr(trk_tandem, "ext") <- "tnd"
attr(trk_tandem, "tracks") <- c("pitch", "voicing_prob")
attr(trk_tandem, "outputType") <- "SSFF"
attr(trk_tandem, "nativeFiletypes") <- c("wav")

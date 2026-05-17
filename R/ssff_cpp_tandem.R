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
    stop("listOfFiles must be a character vector")
  }
  
  # Check file existence
  missing_files <- listOfFiles[!file.exists(listOfFiles)]
  if (length(missing_files) > 0) {
    stop("File(s) not found: ", paste(missing_files, collapse = ", "))
  }
  
  n_files <- length(listOfFiles)
  results <- vector("list", n_files)
  
  if (verbose) format_apply_msg("trk_tandem", n_files)
  
  for (i in seq_along(listOfFiles)) {
    tryCatch({
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
        stop("Failed to load audio file: ", basename(listOfFiles[i]), " — ", e$message)
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
          stop("Resampling failed for ", basename(listOfFiles[i]), " — ", e$message)
        })
      }
    
    # TANDEM requires neural network files in "net/" subdirectory
    # Create temporary net/ directory with symlinks
    net_dir <- file.path(getwd(), "net")
    if (!dir.exists(net_dir)) {
      dir.create(net_dir)
      created_net_dir <- TRUE
    } else {
      created_net_dir <- FALSE
    }
    
    # Symlink or copy network files
    net_source <- system.file("tandem_net", package = "superassp")
    for (net_file in c("MLP1.64.dat", "MLP2.64.dat", "MLP3.64.dat")) {
      src <- file.path(net_source, net_file)
      dst <- file.path(net_dir, net_file)
      if (!file.exists(dst) && file.exists(src)) {
        # Try symlink first, fall back to copy
        tryCatch({
          file.symlink(src, dst)
        }, error = function(e) {
          file.copy(src, dst)
        })
      }
    }
    
    # Ensure cleanup on exit
    on.exit({
      if (created_net_dir && dir.exists(net_dir)) {
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
        stop("TANDEM returned empty pitch track for: ", basename(listOfFiles[i]))
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
      
      tryCatch({
        wrassp::write.AsspDataObj(assp_obj, output_path)
        results[[i]] <- output_path
      }, error = function(e) {
        stop("Failed to write output file: ", output_path, " — ", e$message)
      })
      
    } else {
      results[[i]] <- assp_obj
    }
    }, error = function(e) {
      cli::cli_warn("Error processing {.file {basename(listOfFiles[i])}}: {conditionMessage(e)}")
      results[[i]] <<- if (toFile) FALSE else NULL
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

#' Voice Quality Parameter Extraction via IAIF
#'
#' Extract comprehensive voice quality measures from speech signals using
#' IAIF-based glottal source analysis. Computes parameters related to glottal
#' source characteristics, spectral properties, and voice quality.
#'
#' The function automatically performs IAIF (Iterative Adaptive Inverse Filtering)
#' to extract the glottal source, then computes various voice quality parameters
#' from the glottal flow and its derivative.
#'
#' Implemented in native C++ (no Python dependency).
#'
#' @param listOfFiles Character vector of audio file paths
#' @param beginTime Numeric vector of start times in seconds (default: 0.0)
#' @param endTime Numeric vector of end times in seconds (default: 0.0, full file)
#' @param f0 Optional F0 estimate. Can be:
#'   \itemize{
#'     \item Scalar: Single F0 value for entire signal
#'     \item Vector: F0 contour (median of voiced frames will be used)
#'     \item NULL: H1-H2 will not be computed (returns NA)
#'   }
#' @param gci Optional glottal closure instants. Can be:
#'   \itemize{
#'     \item Numeric vector of sample indices
#'     \item Numeric vector of times in seconds (will be converted to samples)
#'     \item NULL: NAQ and QOQ will not be computed (return NA)
#'   }
#' @param gci_in_samples Logical; if TRUE, gci is in sample indices; if FALSE,
#'   gci is in seconds (default: FALSE)
#' @param verbose Logical; show progress messages (default: TRUE)
#' @param toFile Logical. If TRUE, write results to JSTF file. Default FALSE.
#' @param explicitExt Character. File extension for output. Default "cvq".
#' @param outputDirectory Character. Output directory path. Default NULL (use input directory).
#'
#' @return If \code{toFile=FALSE} (default), for single file: Named list with voice quality parameters.
#'   For multiple files: List of named lists. If \code{toFile=TRUE}, invisibly returns the path(s) to the written JSTF file(s).
#'
#'   Each parameter list contains:
#'   \describe{
#'     \item{\code{glottal_flow_max}}{Peak glottal flow amplitude}
#'     \item{\code{glottal_flow_min}}{Minimum glottal flow value}
#'     \item{\code{glottal_derivative_peak}}{Maximum flow derivative (MFDR)}
#'     \item{\code{NAQ}}{Normalized Amplitude Quotient (NA if no GCI)}
#'     \item{\code{QOQ}}{Quasi-Open Quotient (NA if no GCI)}
#'     \item{\code{H1_H2}}{First two harmonics difference in dB (NA if no F0)}
#'     \item{\code{HRF}}{Harmonic Richness Factor}
#'     \item{\code{PSP}}{Parabolic Spectral Parameter}
#'   }
#'
#' @references
#' \insertCite{Alku1992}{superassp}
#'
#' \insertCite{Kane2013}{superassp}
#'
#' @seealso
#' \code{\link{trk_covarep_iaif}} for glottal waveforms,
#' \code{\link{trk_covarep_srh}} for F0 estimation
#'
#' @examples
#' \dontrun{
#' # Basic usage (HRF and PSP only, no F0 or GCI)
#' vq <- lst_covarep_vq("vowel.wav")
#' print(vq$HRF)
#' print(vq$PSP)
#'
#' # With F0 for H1-H2 computation
#' vq <- lst_covarep_vq("vowel.wav", f0 = 150)
#' print(vq$H1_H2)
#'
#' # Batch processing
#' files <- c("a.wav", "e.wav", "i.wav", "o.wav", "u.wav")
#' vq_all <- lst_covarep_vq(files)
#' }
#'
#' @export
lst_covarep_vq <- function(listOfFiles,
                           beginTime = 0.0,
                           endTime = 0.0,
                           f0 = NULL,
                           gci = NULL,
                           gci_in_samples = FALSE,
                           verbose = TRUE,
                           toFile = FALSE,
                           explicitExt = "cvq",
                           outputDirectory = NULL) {

  # Validate JSTF parameters
  validate_jstf_parameters(toFile, explicitExt, outputDirectory, "lst_covarep_vq")

  # Normalize time parameters
  n_files <- length(listOfFiles)
  beginTime <- fast_recycle_times(beginTime, n_files)
  endTime <- fast_recycle_times(endTime, n_files)

  # Validate time window parameters
  validate_time_window(beginTime, endTime, n_files, "lst_covarep_vq")

  # Validate F0 if provided
  if (!is.null(f0)) {
    if (!is.numeric(f0)) stop("f0 must be numeric (scalar or vector)", call. = FALSE)
    if (any(f0[!is.na(f0)] < 0)) stop("f0 values must be non-negative", call. = FALSE)
  }

  # Validate GCI if provided
  if (!is.null(gci)) {
    if (!is.numeric(gci)) stop("gci must be numeric vector", call. = FALSE)
  }

  # Initialize results
  results <- vector("list", n_files)

  # Progress bar
  if (verbose && n_files > 1) {
    cli::cli_alert_info("Extracting voice quality from {n_files} file{?s}")
    pb <- cli::cli_progress_bar("Voice quality extraction", total = n_files)
  }

  # Process files
  for (i in seq_along(listOfFiles)) {
    file_path <- listOfFiles[i]

    if (!file.exists(file_path)) {
      warning("File not found: ", file_path, call. = FALSE)
      results[i] <- list(NULL)
      if (verbose && n_files > 1) cli::cli_progress_update()
      next
    }

    tryCatch({
      # Load audio via av
      audio_bin <- av::read_audio_bin(
        audio = file_path,
        start_time = if (beginTime[i] > 0) beginTime[i] else NULL,
        end_time = if (endTime[i] > 0) endTime[i] else NULL,
        channels = 1
      )
      sample_rate <- attr(audio_bin, "sample_rate")
      samples <- as.numeric(audio_bin) / 2147483647.0

      # Run IAIF via C++
      iaif_result <- iaif_cpp(samples, as.double(sample_rate))
      glottal_flow <- iaif_result$glottal_flow
      glottal_derivative <- iaif_result$glottal_derivative

      if (length(glottal_flow) == 0) {
        warning(format_processing_warning(file_path, "IAIF returned empty result", "IAIF glottal flow"),
                call. = FALSE)
        results[i] <- list(NULL)
        if (verbose && n_files > 1) cli::cli_progress_update()
        next
      }

      # Prepare F0 parameter
      cpp_f0 <- -1.0
      if (!is.null(f0)) {
        if (length(f0) == 1) {
          cpp_f0 <- as.numeric(f0)
        } else {
          f0_voiced <- f0[f0 > 0 & !is.na(f0)]
          if (length(f0_voiced) > 0) cpp_f0 <- median(f0_voiced)
        }
      }

      # Prepare GCI parameter
      cpp_gci <- NULL
      if (!is.null(gci)) {
        if (gci_in_samples) {
          cpp_gci <- as.integer(gci)
        } else {
          cpp_gci <- as.integer(gci * sample_rate)
        }
        cpp_gci <- cpp_gci[cpp_gci >= 0L & cpp_gci < length(glottal_flow)]
      }

      # Extract VQ params via C++
      results[[i]] <- extract_vq_params_cpp(
        glottal_flow, glottal_derivative,
        as.double(sample_rate), cpp_f0, cpp_gci
      )

    }, error = function(e) {
      warning(format_processing_error(file_path, safe_error_message(e), "voice quality extraction"),
              call. = FALSE)
      results[i] <- list(NULL)
    })

    if (verbose && n_files > 1) cli::cli_progress_update()
  }

  if (verbose && n_files > 1) cli::cli_progress_done()

  # Handle JSTF file writing
  if (toFile) {
    output_paths <- write_lst_results_to_jstf(
      results = results,
      file_paths = listOfFiles,
      beginTime = beginTime,
      endTime = endTime,
      function_name = "lst_covarep_vq",
      parameters = list(f0_provided = !is.null(f0), gci_provided = !is.null(gci)),
      explicitExt = explicitExt,
      outputDirectory = outputDirectory,
      verbose = verbose
    )
    return(invisible(output_paths))
  }

  if (n_files == 1) return(results[[1]]) else return(results)
}

# Set function attributes
attr(lst_covarep_vq, "ext") <- "cvq"
attr(lst_covarep_vq, "outputType") <- "JSTF"
attr(lst_covarep_vq, "format") <- "JSON"

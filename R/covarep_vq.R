#' Voice Quality Parameter Extraction via COVAREP
#'
#' Extract comprehensive voice quality measures from speech signals using
#' COVAREP algorithms. Computes parameters related to glottal source
#' characteristics, spectral properties, and voice quality.
#'
#' The function automatically performs IAIF (Iterative Adaptive Inverse Filtering)
#' to extract the glottal source, then computes various voice quality parameters
#' from the glottal flow and its derivative.
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
#'   F0 can be obtained from \code{trk_covarep_srh()}.
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
#' @details
#' **Voice Quality Parameters:**
#'
#' \describe{
#'   \item{\bold{NAQ (Normalized Amplitude Quotient)}}{
#'     Measures vocal fold closure characteristics.
#'     - Normal range: 0.05-0.20
#'     - Higher values indicate more breathy voice
#'     - Requires GCI detection
#'     - Clinical use: Dysphonia assessment
#'   }
#'   \item{\bold{QOQ (Quasi-Open Quotient)}}{
#'     Ratio of open phase to total glottal cycle duration.
#'     - Normal range: 0.4-0.7
#'     - Related to vocal effort and voice quality
#'     - Requires GCI detection
#'   }
#'   \item{\bold{H1-H2}}{
#'     Spectral tilt measure (dB difference between first two harmonics).
#'     - Positive values: breathy voice
#'     - Negative values: pressed phonation
#'     - Typical range: -10 to +10 dB
#'     - Requires F0 estimate
#'   }
#'   \item{\bold{HRF (Harmonic Richness Factor)}}{
#'     Ratio of high-frequency (>2kHz) to low-frequency energy.
#'     - Indicates spectral richness
#'     - Always computed (no dependencies)
#'   }
#'   \item{\bold{PSP (Parabolic Spectral Parameter)}}{
#'     Spectral envelope curvature coefficient.
#'     - Related to overall spectral shape
#'     - Always computed (no dependencies)
#'   }
#' }
#'
#' **Performance:** With Numba optimization, processes typical sustained vowel
#' (3s) in ~30-50ms total (IAIF + VQ extraction).
#'
#' **Optimization:** Uses Numba JIT for NAQ/QOQ computation (2-3x speedup).
#' Check status with \code{covarep_info()}.
#'
#' @references
#' Alku, P. (1992). "Glottal wave analysis with pitch synchronous iterative
#' adaptive inverse filtering". Speech Communication, 11(2-3), 109-118.
#'
#' Kane, J., & Gobl, C. (2013). "Evaluation of glottal closure instant
#' detection in a range of voice qualities". Speech Communication, 55(2), 295-314.
#'
#' @seealso
#' \code{\link{trk_covarep_srh}} for F0 estimation,
#' \code{\link{trk_covarep_iaif}} for glottal waveforms,
#' \code{\link{install_covarep}} for installation,
#' \code{\link{covarep_info}} for optimization status
#'
#' @examples
#' \dontrun{
#' # Basic usage (HRF and PSP only, no F0 or GCI)
#' vq <- lst_covarep_vq("vowel.wav")
#' print(vq$HRF)
#' print(vq$PSP)
#'
#' # With F0 for H1-H2 computation
#' f0_data <- trk_covarep_srh("vowel.wav", toFile = FALSE)
#' f0_mean <- median(f0_data$`F0[Hz]`[f0_data$VUV == 1, 1])
#' vq <- lst_covarep_vq("vowel.wav", f0 = f0_mean)
#' print(vq$H1_H2)  # Now computed
#'
#' # With F0 contour (median of voiced frames used)
#' vq <- lst_covarep_vq("vowel.wav", f0 = f0_data$`F0[Hz]`[, 1])
#'
#' # Batch processing
#' files <- c("a.wav", "e.wav", "i.wav", "o.wav", "u.wav")
#' vq_all <- lst_covarep_vq(files)
#'
#' # Extract specific parameters from batch
#' h1h2_values <- sapply(vq_all, function(x) x$H1_H2)
#' hrf_values <- sapply(vq_all, function(x) x$HRF)
#'
#' # With time windowing
#' vq_window <- lst_covarep_vq("speech.wav", beginTime = 1.0, endTime = 2.5)
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


  # Check COVAREP availability
  if (!covarep_available()) {
    stop("COVAREP Python module not available.\n",
         "Install with: install_covarep()\n",
         "Check status with: covarep_info()",
         call. = FALSE)
  }

  # Normalize time parameters
  n_files <- length(listOfFiles)
  beginTime <- fast_recycle_times(beginTime, n_files)
  endTime <- fast_recycle_times(endTime, n_files)


  # Validate time window parameters
  validate_time_window(beginTime, endTime, n_files, "lst_covarep_vq")

  # Validate F0 if provided
  if (!is.null(f0)) {
    if (!is.numeric(f0)) {
      stop("f0 must be numeric (scalar or vector)", call. = FALSE)
    }
    if (any(f0[!is.na(f0)] < 0)) {
      stop("f0 values must be non-negative", call. = FALSE)
    }
  }

  # Validate GCI if provided
  if (!is.null(gci)) {
    if (!is.numeric(gci)) {
      stop("gci must be numeric vector", call. = FALSE)
    }
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

    # Validate file exists
    if (!file.exists(file_path)) {
      warning("File not found: ", file_path, call. = FALSE)
      results[[i]] <- NULL
      if (verbose && n_files > 1) cli::cli_progress_update()
      next
    }

    tryCatch({
      # Load audio
      audio_data <- av_load_for_python(
        file_path,
        start_time = beginTime[i],
        end_time = endTime[i]
      )

      # Compute IAIF to get glottal flow
      iaif_result <- covarep_module$glottal$iaif_optimized$iaif_optimized(
        x = audio_data$samples,
        fs = as.integer(audio_data$sample_rate)
      )

      glottal_flow <- as.numeric(iaif_result[[1]])
      glottal_derivative <- as.numeric(iaif_result[[2]])

      # Check IAIF success
      if (length(glottal_flow) == 0) {

        warning(format_processing_warning(file_path, "IAIF computation returned empty result", "COVAREP glottal flow"),
                call. = FALSE)

        results[[i]] <- NULL
        if (verbose && n_files > 1) cli::cli_progress_update()
        next
      }

      # Prepare F0 parameter
      py_f0 <- NULL
      if (!is.null(f0)) {
        if (length(f0) == 1) {
          # Scalar F0
          py_f0 <- as.numeric(f0)
        } else {
          # Vector F0 - use median of positive values
          f0_voiced <- f0[f0 > 0 & !is.na(f0)]
          if (length(f0_voiced) > 0) {
            py_f0 <- as.numeric(median(f0_voiced))
          }
        }
      }

      # Prepare GCI parameter
      py_gci <- NULL
      if (!is.null(gci)) {
        if (gci_in_samples) {
          # Already in sample indices
          py_gci <- as.integer(gci)
        } else {
          # Convert from seconds to sample indices
          py_gci <- as.integer(gci * audio_data$sample_rate)
        }
        # Filter to valid range
        py_gci <- py_gci[py_gci >= 0 & py_gci < length(glottal_flow)]
      }

      # Extract voice quality parameters
      vq_params <- covarep_module$glottal$vq_optimized$extract_vq_params_optimized(
        glottal_flow = glottal_flow,
        glottal_derivative = glottal_derivative,
        fs = as.integer(audio_data$sample_rate),
        f0 = py_f0,
        gci = py_gci
      )

      # Convert Python dict to R list
      param_list <- list()
      param_names <- names(vq_params)

      for (name in param_names) {
        value <- vq_params[[name]]
        # Convert numpy scalars to R numerics
        param_list[[name]] <- as.numeric(value)
      }

      results[[i]] <- param_list

    }, error = function(e) {

      warning(format_processing_error(file_path, safe_error_message(e), "COVAREP voice quality extraction"),
              call. = FALSE)

      results[[i]] <- NULL
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
      parameters = list(
        f0_provided = !is.null(f0),
        gci_provided = !is.null(gci)
      ),
      explicitExt = explicitExt,
      outputDirectory = outputDirectory,
      verbose = verbose
    )
    return(invisible(output_paths))
  }


  # Return results
  if (n_files == 1) {
    return(results[[1]])
  } else {
    return(results)
  }
}


# Set function attributes
attr(lst_covarep_vq, "ext") <- "cvq"
attr(lst_covarep_vq, "outputType") <- "JSTF"
attr(lst_covarep_vq, "format") <- "JSON"


#' SRH F0 Tracking via COVAREP Python
#'
#' Estimate fundamental frequency using the Summation of Residual
#' Harmonics (SRH) algorithm from COVAREP.
#'
#' This function implements the SRH F0 estimation algorithm which:
#' \enumerate{
#'   \item Computes LP residual signal
#'   \item Analyzes harmonic structure via spectrogram
#'   \item Selects F0 based on harmonic summation criterion
#'   \item Provides voice/unvoiced decision and confidence measure
#' }
#'
#' Performance: With Numba optimization, processes 10s audio in ~67ms
#' (143x real-time, 7x speedup over unoptimized).
#'
#' @param listOfFiles Character vector of file paths to audio files
#' @param beginTime Numeric vector of start times in seconds (default: 0.0)
#' @param endTime Numeric vector of end times in seconds (default: 0.0, full file)
#' @param windowShift Frame shift in milliseconds (default: 5.0)
#' @param minF Minimum F0 in Hz (default: 50)
#' @param maxF Maximum F0 in Hz (default: 500)
#' @param toFile Logical; if TRUE, write SSFF files; if FALSE, return AsspDataObj (default: TRUE)
#' @param explicitExt File extension for output files (default: "f0")
#' @param outputDirectory Optional output directory for files
#' @param verbose Logical; show progress messages (default: TRUE)
#' @param ... Additional arguments (parallel, n_cores for batch processing)
#'
#' @return If toFile=FALSE: AsspDataObj or list of AsspDataObj with tracks:
#'   \itemize{
#'     \item \code{F0[Hz]}: Estimated fundamental frequency
#'     \item \code{VUV}: Voice/unvoiced decision (0=unvoiced, 1=voiced)
#'     \item \code{SRH}: SRH confidence values
#'   }
#'   If toFile=TRUE: Character vector of output file paths
#'
#' @details
#' The SRH algorithm analyzes the harmonic structure of the LP residual signal.
#' It computes the sum of harmonic energies minus subharmonic energies across
#' candidate F0 values. The candidate with maximum SRH value is selected as F0.
#'
#' **Optimization:** This implementation uses NumPy vectorization and optional
#' Numba JIT compilation for 5-10x speedup. Check optimization status with
#' \code{covarep_info()}.
#'
#' **Reference:**
#' T. Drugman and A. Alwan (2011). "Joint Robust Voicing Detection and Pitch
#' Estimation Based on Residual Harmonics". Interspeech 2011.
#'
#' @seealso
#' \code{\link{trk_covarep_iaif}} for glottal analysis,
#' \code{\link{install_covarep}} for installation,
#' \code{\link{covarep_info}} for optimization status
#'
#' @examples
#' \dontrun{
#' # Single file
#' f0 <- trk_covarep_srh("audio.wav", toFile = FALSE)
#' plot(f0$F0, type = "l", ylab = "F0 (Hz)")
#'
#' # Batch processing with time windowing
#' files <- c("a.wav", "b.wav", "c.wav")
#' result <- trk_covarep_srh(files,
#'                           beginTime = c(0.5, 1.0, 0.0),
#'                           endTime = c(3.0, 4.0, 2.0),
#'                           minF = 75, maxF = 300)
#'
#' # Custom F0 range for female speaker
#' f0_female <- trk_covarep_srh("female.wav",
#'                              minF = 150, maxF = 400,
#'                              toFile = FALSE)
#' }
#'
#' @export
trk_covarep_srh <- function(listOfFiles,
                            beginTime = 0.0,
                            endTime = 0.0,
                            windowShift = 5.0,
                            minF = 50,
                            maxF = 500,
                            toFile = TRUE,
                            explicitExt = "f0",
                            outputDirectory = NULL,
                            verbose = TRUE,
                            ...) {

  # Check Python module availability
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

  # Validate parameters
  if (windowShift <= 0) {
    stop("windowShift must be positive", call. = FALSE)
  }
  if (minF <= 0 || maxF <= minF) {
    stop("F0 range invalid: minF must be positive and less than maxF", call. = FALSE)
  }

  # Initialize results
  results <- vector("list", n_files)

  # Progress bar
  if (verbose && n_files > 1) {
    cli::cli_alert_info("Processing {n_files} file{?s} with COVAREP SRH")
    pb <- cli::cli_progress_bar("SRH F0 tracking", total = n_files)
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
      # Load audio to Python-compatible format
      audio_data <- av_load_for_python(
        file_path,
        start_time = beginTime[i],
        end_time = endTime[i]
      )

      # Call Python function (optimized version)
      py_result <- covarep_module$f0$f0_optimized$pitch_srh_vectorized(
        wave = audio_data$samples,
        fs = as.integer(audio_data$sample_rate),
        f0min = minF,
        f0max = maxF,
        hopsize = windowShift
      )

      # Extract results (Python returns tuple: f0, vuv, srh_values, times)
      f0_values <- as.numeric(py_result[[1]])
      vuv_values <- as.integer(py_result[[2]])
      srh_values <- as.numeric(py_result[[3]])
      times <- as.numeric(py_result[[4]])

      # Check for empty results
      if (length(f0_values) == 0) {
        warning("No F0 detected in ", basename(file_path), call. = FALSE)
        results[[i]] <- NULL
        if (verbose && n_files > 1) cli::cli_progress_update()
        next
      }

      # Convert to AsspDataObj
      obj <- list()
      obj$`F0[Hz]` <- matrix(f0_values, ncol = 1)
      obj$VUV <- matrix(vuv_values, ncol = 1)
      obj$SRH <- matrix(srh_values, ncol = 1)

      # Set attributes
      # Frame rate from windowShift (ms)
      frame_rate <- 1000.0 / windowShift
      attr(obj, "sampleRate") <- frame_rate
      attr(obj, "startTime") <- beginTime[i]
      attr(obj, "startRecord") <- 1L
      attr(obj, "endRecord") <- as.integer(length(f0_values))
      attr(obj, "trackFormats") <- c("REAL64", "INT16", "REAL64")
      attr(obj, "origFreq") <- audio_data$sample_rate
      class(obj) <- "AsspDataObj"

      # Write to file if requested
      if (toFile) {
        # Construct output path
        if (is.null(outputDirectory)) {
          out_dir <- dirname(file_path)
        } else {
          out_dir <- outputDirectory
          if (!dir.exists(out_dir)) {
            dir.create(out_dir, recursive = TRUE)
          }
        }

        base_name <- tools::file_path_sans_ext(basename(file_path))
        out_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))

        wrassp::write.AsspDataObj(obj, out_path)
        results[[i]] <- out_path
      } else {
        results[[i]] <- obj
      }

    }, error = function(e) {
      warning("Error processing ", basename(file_path), ": ",
              e$message, call. = FALSE)
      results[[i]] <- NULL
    })

    if (verbose && n_files > 1) cli::cli_progress_update()
  }

  if (verbose && n_files > 1) cli::cli_progress_done()

  # Return results
  if (n_files == 1) {
    return(results[[1]])
  } else {
    return(results)
  }
}

# Set function attributes for consistency with other DSP functions
attr(trk_covarep_srh, "ext") <- "f0"
attr(trk_covarep_srh, "tracks") <- c("F0[Hz]", "VUV", "SRH")
attr(trk_covarep_srh, "outputType") <- "SSFF"

#' Glottal Source Analysis via IAIF
#'
#' Extract glottal flow waveform and derivative using Iterative Adaptive
#' Inverse Filtering (IAIF).
#'
#' IAIF iteratively estimates the vocal tract filter and glottal source:
#' \enumerate{
#'   \item Rough glottal source estimation (high-order LPC)
#'   \item Vocal tract filter estimation from glottal estimate
#'   \item Refined glottal source extraction via inverse filtering
#' }
#'
#' Implemented in native C++ (no Python dependency).
#'
#' @param listOfFiles Character vector of file paths to audio files
#' @param beginTime Numeric vector of start times in seconds (default: 0.0)
#' @param endTime Numeric vector of end times in seconds (default: 0.0, full file)
#' @param order_vt Vocal tract LPC order (default: NULL = auto: 2*round(fs/2000)+4)
#' @param order_gl Glottal source LPC order (default: NULL = auto: 2*round(fs/4000))
#' @param leaky_coef Leaky integration coefficient for lip radiation compensation (default: 0.99)
#' @param hpfilt Apply high-pass filter (default: TRUE)
#' @param toFile Logical; if TRUE, write SSFF files; if FALSE, return AsspDataObj (default: TRUE)
#' @param explicitExt File extension for output files (default: "glf")
#' @param outputDirectory Optional output directory for files
#' @param verbose Logical; show progress messages (default: TRUE)
#' @param ... Additional arguments (parallel, n_cores for batch processing)
#'
#' @return If toFile=FALSE: AsspDataObj or list of AsspDataObj with tracks:
#'   \itemize{
#'     \item \code{glottal_flow}: Estimated glottal flow waveform
#'     \item \code{glottal_derivative}: Glottal flow derivative (MFDR)
#'   }
#'   If toFile=TRUE: Character vector of output file paths
#'
#' @details
#' The IAIF algorithm separates the glottal source contribution from the
#' vocal tract filtering effect in speech signals. The output glottal flow
#' waveform can be used for voice quality analysis, dysphonia assessment,
#' and speaker characterization.
#'
#' **Filter Orders:**
#' The default orders follow the original COVAREP/MATLAB implementation:
#' \itemize{
#'   \item Vocal tract: \code{2 * round(fs/2000) + 4} (typically 12-20)
#'   \item Glottal source: \code{2 * round(fs/4000)} (typically 2-4)
#' }
#'
#' @references
#' \insertCite{Alku1992}{superassp}
#'
#' @seealso
#' \code{\link{lst_covarep_vq}} for voice quality parameters,
#' \code{\link{trk_covarep_srh}} for F0 estimation
#'
#' @examples
#' \dontrun{
#' # Single file - extract glottal flow
#' glottal <- trk_covarep_iaif("vowel.wav", toFile = FALSE)
#' plot(glottal$glottal_flow, type = "l", ylab = "Glottal Flow")
#' plot(glottal$glottal_derivative, type = "l", ylab = "Glottal Derivative")
#'
#' # Batch processing
#' files <- c("a.wav", "e.wav", "i.wav")
#' result <- trk_covarep_iaif(files, toFile = TRUE)
#'
#' # Custom filter orders
#' glottal_custom <- trk_covarep_iaif("audio.wav",
#'                                    order_vt = 16,
#'                                    order_gl = 3,
#'                                    toFile = FALSE)
#' }
#'
#' @export
trk_covarep_iaif <- function(listOfFiles,
                             beginTime = 0.0,
                             endTime = 0.0,
                             order_vt = NULL,
                             order_gl = NULL,
                             leaky_coef = 0.99,
                             hpfilt = TRUE,
                             toFile = TRUE,
                             explicitExt = "glf",
                             outputDirectory = NULL,
                             verbose = TRUE,
                             ...) {

  # Normalize time parameters
  n_files <- length(listOfFiles)
  beginTime <- fast_recycle_times(beginTime, n_files)
  endTime <- fast_recycle_times(endTime, n_files)

  # Validate parameters
  if (leaky_coef <= 0 || leaky_coef >= 1) {
    stop("leaky_coef must be between 0 and 1 (typically 0.95-0.99)", call. = FALSE)
  }

  # Initialize results
  results <- vector("list", n_files)

  # Progress bar
  if (verbose && n_files > 1) {
    cli::cli_alert_info("Processing {n_files} file{?s} with IAIF (C++)")
    pb <- cli::cli_progress_bar("Glottal analysis", total = n_files)
  }

  # Process files
  for (i in seq_along(listOfFiles)) {
    file_path <- listOfFiles[i]

    # Validate file exists
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

      # Convert INT32 to float64 [-1, 1]
      samples <- as.numeric(audio_bin) / 2147483647.0

      # C++ parameters (-1 means auto)
      cpp_order_vt <- if (is.null(order_vt)) -1L else as.integer(order_vt)
      cpp_order_gl <- if (is.null(order_gl)) -1L else as.integer(order_gl)

      # Call C++ IAIF
      cpp_result <- iaif_cpp(samples, as.double(sample_rate),
                             cpp_order_vt, cpp_order_gl,
                             leaky_coef, hpfilt)

      glottal_flow <- cpp_result$glottal_flow
      glottal_derivative <- cpp_result$glottal_derivative

      # Check for empty results
      if (length(glottal_flow) == 0) {
        warning("IAIF failed for ", basename(file_path), call. = FALSE)
        results[i] <- list(NULL)
        if (verbose && n_files > 1) cli::cli_progress_update()
        next
      }

      # Convert to AsspDataObj
      obj <- list()
      obj$glottal_flow <- matrix(glottal_flow, ncol = 1)
      obj$glottal_derivative <- matrix(glottal_derivative, ncol = 1)

      attr(obj, "sampleRate") <- as.numeric(sample_rate)
      attr(obj, "startTime") <- as.numeric(beginTime[i])
      attr(obj, "startRecord") <- 1L
      attr(obj, "endRecord") <- as.integer(length(glottal_flow))
      attr(obj, "trackFormats") <- c("REAL64", "REAL64")
      attr(obj, "origFreq") <- as.numeric(sample_rate)
      class(obj) <- "AsspDataObj"
      AsspFileFormat(obj) <- "SSFF"
      AsspDataFormat(obj) <- as.integer(2)

      # Write to file if requested
      if (toFile) {
        out_dir <- if (is.null(outputDirectory)) dirname(file_path) else outputDirectory
        if (!is.null(outputDirectory) && !dir.exists(out_dir)) {
          dir.create(out_dir, recursive = TRUE)
        }
        base_name <- tools::file_path_sans_ext(basename(file_path))
        out_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))
        write.AsspDataObj(obj, out_path)
        results[[i]] <- out_path
      } else {
        results[[i]] <- obj
      }

    }, error = function(e) {
      warning("Error processing ", basename(file_path), ": ",
              e$message, call. = FALSE)
      results[i] <- list(NULL)
    })

    if (verbose && n_files > 1) cli::cli_progress_update()
  }

  if (verbose && n_files > 1) cli::cli_progress_done()

  if (n_files == 1) return(results[[1]]) else return(results)
}

# Set function attributes
attr(trk_covarep_iaif, "ext") <- "glf"
attr(trk_covarep_iaif, "tracks") <- c("glottal_flow", "glottal_derivative")
attr(trk_covarep_iaif, "outputType") <- "SSFF"

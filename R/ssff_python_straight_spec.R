#' STRAIGHT Spectral Analysis (Legacy Vocoder)
#'
#' Extracts pitch-adaptive spectral envelope using the legacy STRAIGHT algorithm.
#' Provides 99% correlation with the original MATLAB implementation.
#'
#' @param listOfFiles Character vector of audio file paths, or AVAudio object
#' @param beginTime Numeric; start time in seconds (default: 0, start of file)
#' @param endTime Numeric; end time in seconds (default: 0, end of file)
#' @param f0_floor Numeric; minimum F0 in Hz (default: 40)
#' @param f0_ceil Numeric; maximum F0 in Hz (default: 800)
#' @param fft_size Numeric; FFT size for spectral analysis (default: 2048)
#' @param frame_shift Numeric; frame shift in milliseconds (default: 1.0)
#' @param toFile Logical; write output to SSFF file (default: FALSE for in-memory)
#' @param explicitExt Character; file extension for output (default: "strspec")
#' @param outputDirectory Character; output directory (default: NULL, same as input)
#' @param verbose Logical; print progress messages (default: TRUE)
#' @param ... Additional arguments (for S7 dispatch compatibility)
#'
#' @return If `toFile = FALSE`: AsspDataObj (single file) or list of AsspDataObj
#'   (multiple files) with spectral tracks.
#'   If `toFile = TRUE`: Character vector of output file paths (returned invisibly).
#'
#' @details
#' The STRAIGHT spectral analysis algorithm performs pitch-adaptive smoothing
#' of the time-frequency representation, providing a high-quality spectral
#' envelope suitable for vocoding and voice conversion.
#'
#' **Features**:
#' - Pitch-adaptive smoothing (no voicing artifacts)
#' - High-resolution spectral representation
#' - Compatible with STRAIGHT synthesis
#'
#' **Accuracy**: 99% correlation with MATLAB STRAIGHT spectral output
#'
#' **Audio Loading**: Uses `av` package for universal media support.
#'
#' @references
#' \insertRef{Kawahara.1999.SpeechCommunication}{superassp}
#'
#' @examples
#' \dontrun{
#' # Check availability
#' if (!straight_available()) {
#'   install_legacy_straight()
#' }
#'
#' # Single file spectral analysis
#' wav_file <- system.file("samples", "sustained", "a1.wav", package = "superassp")
#' spec_data <- trk_straight_spec(wav_file, toFile = FALSE)
#'
#' # Plot spectrogram
#' image(t(spec_data$spec), col = heat.colors(256))
#'
#' # Batch processing
#' files <- list.files("audio/", pattern = "\\.wav$", full.names = TRUE)
#' spec_results <- trk_straight_spec(files, toFile = FALSE)
#' }
#'
#' @seealso
#' \code{\link{trk_straight_f0}}, \code{\link{straight_synth}},
#' \code{\link{install_legacy_straight}}
#'
#' @export
trk_straight_spec <- function(listOfFiles, beginTime = 0.0, endTime = 0.0,
                              f0_floor = 40, f0_ceil = 800,
                              fft_size = 2048, frame_shift = 1.0,
                              toFile = FALSE, explicitExt = "strspec",
                              outputDirectory = NULL, verbose = TRUE, ...) {
  
  .Deprecated("trk_d4c", package = "superassp",
    msg = "trk_straight_spec() is deprecated. Use trk_d4c() (WORLD C++) for spectral analysis.")

  # Check if STRAIGHT is available
  if (!straight_available()) {
    stop("Legacy STRAIGHT not available. Install with: install_legacy_straight()",
         call. = FALSE)
  }
  
  # Handle single file
  if (is.character(listOfFiles) && length(listOfFiles) == 1) {
    return(.straight_spec_single(
      file_path = listOfFiles,
      beginTime = beginTime,
      endTime = endTime,
      f0_floor = f0_floor,
      f0_ceil = f0_ceil,
      fft_size = fft_size,
      frame_shift = frame_shift,
      toFile = toFile,
      explicitExt = explicitExt,
      outputDirectory = outputDirectory,
      verbose = verbose
    ))
  }
  
  # Batch processing
  n_files <- length(listOfFiles)
  
  if (verbose) {
    message(sprintf("Processing %d files with STRAIGHT spectral analysis...", n_files))
  }
  
  results <- vector("list", n_files)
  
  if (verbose) {
    pb <- utils::txtProgressBar(min = 0, max = n_files, style = 3)
  }
  
  for (i in seq_along(listOfFiles)) {
    results[[i]] <- .straight_spec_single(
      file_path = listOfFiles[i],
      beginTime = if (length(beginTime) == n_files) beginTime[i] else beginTime[1],
      endTime = if (length(endTime) == n_files) endTime[i] else endTime[1],
      f0_floor = f0_floor,
      f0_ceil = f0_ceil,
      fft_size = fft_size,
      frame_shift = frame_shift,
      toFile = toFile,
      explicitExt = explicitExt,
      outputDirectory = outputDirectory,
      verbose = FALSE
    )
    
    if (verbose) {
      utils::setTxtProgressBar(pb, i)
    }
  }
  
  if (verbose) {
    close(pb)
    message("\nProcessing complete!")
  }
  
  if (toFile) {
    return(invisible(unlist(results)))
  } else {
    return(results)
  }
}


#' Internal function for single-file STRAIGHT spectral analysis
#' @keywords internal
#' @noRd
.straight_spec_single <- function(file_path, beginTime, endTime, f0_floor,
                                  f0_ceil, fft_size, frame_shift, toFile,
                                  explicitExt, outputDirectory, verbose) {
  
  # Load audio (av primary, read_audio fallback for non-av formats)
  audio_data <- tryCatch(
    av::read_audio_bin(
      audio = file_path,
      start_time = if (beginTime > 0) beginTime else NULL,
      end_time = if (endTime > 0) endTime else NULL,
      channels = 1
    ),
    error = function(e) {
      assp_obj  <- read_audio(file_path, begin = beginTime, end = endTime)
      raw_int32 <- as.integer(assp_obj$audio[, 1]) * 65536L
      attr(raw_int32, "sample_rate") <- attr(assp_obj, "sampleRate")
      attr(raw_int32, "channels")    <- 1L
      raw_int32
    }
  )
  
  sample_rate <- attr(audio_data, "sample_rate")
  audio_float <- as.numeric(audio_data) / 2147483647.0
  
  # Import Python modules
  .setup_straight_path <- get(".setup_straight_path", envir = asNamespace("superassp"))
  .setup_straight_path()
  
  py_spec <- reticulate::import("legacy_STRAIGHT.spectral")
  py_f0 <- reticulate::import("legacy_STRAIGHT.f0_extraction")
  np <- reticulate::import("numpy", convert = FALSE)
  
  audio_np <- np$array(audio_float, dtype = "float32")
  
  if (verbose) {
    message("Running STRAIGHT spectral analysis...")
  }
  
  # First extract F0 (required for pitch-adaptive smoothing)
  f0_result <- py_f0$MulticueF0v14(
    x = audio_np,
    fs = as.integer(sample_rate),
    f0_floor = f0_floor,
    f0_ceil = f0_ceil,
    frame_shift = frame_shift / 1000.0
  )
  
  f0_values <- reticulate::py_to_r(f0_result[[1]])
  
  # Run spectral analysis
  spec_result <- py_spec$exstraightspec(
    x = audio_np,
    fs = as.integer(sample_rate),
    f0 = np$array(f0_values, dtype = "float64"),
    fft_size = as.integer(fft_size),
    frame_shift = frame_shift / 1000.0
  )
  
  # Extract spectral data
  spec_data <- reticulate::py_to_r(spec_result$spectrogram)
  time_axis <- reticulate::py_to_r(spec_result$time_axis)
  freq_axis <- reticulate::py_to_r(spec_result$freq_axis)
  
  n_frames <- ncol(spec_data)
  n_freqs <- nrow(spec_data)
  
  # Create AsspDataObj
  assp_obj <- list()
  assp_obj$spec <- t(spec_data)  # Transpose to [frames x freqs]
  
  attr(assp_obj, "sampleRate") <- sample_rate
  attr(assp_obj, "origFreq") <- sample_rate
  attr(assp_obj, "startTime") <- beginTime
  attr(assp_obj, "startRecord") <- 1L
  attr(assp_obj, "endRecord") <- as.integer(n_frames)
  attr(assp_obj, "trackFormats") <- "REAL32"
  
  class(assp_obj) <- "AsspDataObj"
  
  # Write to file if requested
  if (toFile) {
    base_name <- tools::file_path_sans_ext(basename(file_path))
    out_dir <- if (!is.null(outputDirectory)) outputDirectory else dirname(file_path)
    output_path <- file.path(out_dir, paste0(base_name, ".", explicitExt))
    
    wrassp::write.AsspDataObj(assp_obj, output_path)
    
    if (verbose) {
      message("Written: ", output_path)
    }
    
    return(invisible(output_path))
  }
  
  return(assp_obj)
}


# Set function attributes
attr(trk_straight_spec, "ext") <- "strspec"
attr(trk_straight_spec, "tracks") <- c("spec")
attr(trk_straight_spec, "outputType") <- "SSFF"
attr(trk_straight_spec, "nativeFiletypes") <- c("wav")

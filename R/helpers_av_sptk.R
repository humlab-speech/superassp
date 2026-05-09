##' SPTK Helper Functions
##'
##' @description Internal helper functions shared across SPTK-based DSP functions
##' @name sptk_helpers
##' @keywords internal
NULL


##' Convert SPTK C++ pitch result to AsspDataObj
##'
##' @param pitch_result List returned from rapt_cpp, swipe_cpp, reaper_cpp, or dio_cpp
##' @param windowShift Frame shift in milliseconds
##' @return AsspDataObj with f0 track
##' @keywords internal
create_f0_asspobj <- function(pitch_result, windowShift) {
  # Extract F0 matrix and metadata
  f0_matrix <- pitch_result$f0
  sample_rate <- as.numeric(pitch_result$sample_rate)  # Original sample rate
  n_frames <- as.integer(pitch_result$n_frames)  # Ensure integer
  
  # Calculate effective frame rate
  frame_rate <- 1000.0 / windowShift  # windowShift is in ms

  # Create AsspDataObj
  out_obj <- list(
    f0 = f0_matrix
  )

  # Set attributes matching wrassp format
  attr(out_obj, "trackFormats") <- "REAL32"
  attr(out_obj, "sampleRate") <- frame_rate  # Frames per second
  attr(out_obj, "origFreq") <- sample_rate  # Original audio sample rate
  attr(out_obj, "startTime") <- 0.0
  attr(out_obj, "startRecord") <- 1L
  attr(out_obj, "endRecord") <- n_frames
  attr(out_obj, "fileInfo") <- c(20L, 2L)  # SSFF format

  class(out_obj) <- "AsspDataObj"

  return(out_obj)
}


##' Convert SPTK REAPER epoch times to pitch mark AsspDataObj
##'
##' @param epoch_times Vector of epoch times in seconds (irregular intervals)
##' @param sample_rate Original audio sample rate
##' @param windowShift Frame shift in milliseconds
##' @return AsspDataObj with pm (pitch mark) track as binary indicator
##' @keywords internal
create_pitchmark_asspobj <- function(epoch_times, sample_rate, windowShift) {
  # Calculate frame parameters
  frame_rate <- 1000.0 / windowShift  # Frames per second
  frame_shift_sec <- windowShift / 1000.0  # Frame shift in seconds

  # Determine number of frames needed
  if (length(epoch_times) == 0) {
    # No epochs - return empty track
    n_frames <- 1L
    pm_values <- matrix(0L, nrow = 1, ncol = 1)
  } else {
    # Calculate time range
    max_time <- max(epoch_times)
    n_frames <- as.integer(ceiling(max_time / frame_shift_sec)) + 1L

    # Create binary indicator grid at regular intervals
    pm_values <- matrix(0L, nrow = n_frames, ncol = 1)

    # Mark frames that contain epochs
    for (epoch_time in epoch_times) {
      frame_idx <- as.integer(floor(epoch_time / frame_shift_sec)) + 1L
      if (frame_idx > 0 && frame_idx <= n_frames) {
        pm_values[frame_idx, 1] <- 1L  # Mark this frame as containing pitch mark
      }
    }
  }

  # Create AsspDataObj
  out_obj <- list(
    pm = pm_values
  )

  # Set attributes matching reaper_pm format
  attr(out_obj, "trackFormats") <- "INT16"
  attr(out_obj, "sampleRate") <- frame_rate  # Frames per second
  attr(out_obj, "origFreq") <- sample_rate  # Original audio sample rate
  attr(out_obj, "startTime") <- 0.0
  attr(out_obj, "startRecord") <- 1L
  attr(out_obj, "endRecord") <- n_frames
  attr(out_obj, "fileInfo") <- c(20L, 2L)  # SSFF format

  class(out_obj) <- "AsspDataObj"

  # Store raw epoch times as attribute for advanced users
  attr(out_obj, "epoch_times") <- epoch_times
  attr(out_obj, "n_epochs") <- length(epoch_times)

  return(out_obj)
}


##' Generate output file path for SSFF files
##'
##' @param input_file Input file path
##' @param ext Output extension
##' @param output_dir Output directory (optional)
##' @return Output file path
##' @keywords internal
generate_output_path <- function(input_file, ext, output_dir = NULL) {
  if (is.null(output_dir)) {
    out_dir <- dirname(input_file)
  } else {
    out_dir <- output_dir
  }

  base_name <- fast_file_path_sans_ext(c(basename(input_file)))[1]
  output_file <- file.path(out_dir, paste0(base_name, ".", ext))

  return(output_file)
}


##' Convert SPTK C++ aperiodicity result to AsspDataObj
##'
##' @param ap_result List returned from d4c_cpp
##' @param windowShift Frame shift in milliseconds
##' @return AsspDataObj with aperiodicity track
##' @keywords internal
create_aperiodicity_asspobj <- function(ap_result, windowShift) {
  # Extract aperiodicity matrix and metadata
  ap_matrix <- ap_result$aperiodicity
  sample_rate <- as.numeric(ap_result$sample_rate)  # Original sample rate
  n_frames <- as.integer(ap_result$n_frames)  # Ensure integer
  
  # Calculate effective frame rate
  frame_rate <- 1000.0 / windowShift  # windowShift is in ms

  # Create AsspDataObj
  out_obj <- list(
    aperiodicity = ap_matrix
  )

  # Set attributes matching wrassp format
  attr(out_obj, "trackFormats") <- "REAL32"
  attr(out_obj, "sampleRate") <- frame_rate  # Frames per second
  attr(out_obj, "origFreq") <- sample_rate  # Original audio sample rate
  attr(out_obj, "startTime") <- 0.0
  attr(out_obj, "startRecord") <- 1L
  attr(out_obj, "endRecord") <- n_frames
  attr(out_obj, "fileInfo") <- c(20L, 2L)  # SSFF format

  class(out_obj) <- "AsspDataObj"

  return(out_obj)
}

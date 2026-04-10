##' SPTK Helper Functions
##'
##' @description Internal helper functions shared across SPTK-based DSP functions
##' @name sptk_helpers
##' @keywords internal
NULL


##' Emit a consistent "Applying <fun>()" progress message
##'
##' Prints either:
##'   "Applying `fun()` to N recording(s)"
##' or, when a time window is active:
##'   "Applying `fun()` to a X.X second long portion of N recording(s)"
##'
##' @param fun_name Character. Public function name, e.g. "trk_rapt".
##' @param n_files  Integer. Number of files being processed.
##' @param beginTime Numeric vector of begin times (seconds).
##' @param endTime   Numeric vector of end times (seconds; 0 = end of file).
##' @keywords internal
format_apply_msg <- function(fun_name, n_files, beginTime = NULL, endTime = NULL) {
  has_window <- !is.null(endTime) && any(endTime > 0, na.rm = TRUE)
  if (has_window) {
    durs <- ifelse(endTime > 0, endTime - beginTime, NA_real_)
    durs <- durs[!is.na(durs)]
    if (length(durs) > 0) {
      d <- round(mean(durs), 1)
      cli::cli_inform(
        "Applying {.fun {fun_name}} to a {d} second long portion of {cli::no(n_files)} recording{?s}"
      )
      return(invisible(NULL))
    }
  }
  cli::cli_inform("Applying {.fun {fun_name}} to {cli::no(n_files)} recording{?s}")
}


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


##' Convert Snack pitch result (4 tracks) to AsspDataObj
##'
##' @param res List from snackp_cpp (f0, voicing, rms, acpeak, sample_rate, n_frames)
##' @param windowShift Frame shift in milliseconds
##' @return AsspDataObj with f0, voicing, rms, acpeak tracks
##' @keywords internal
create_snackp_asspobj <- function(res, windowShift) {
  n_frames   <- as.integer(res$n_frames)
  frame_rate <- 1000.0 / windowShift

  out_obj <- list(
    f0      = res$f0,
    voicing = res$voicing,
    rms     = res$rms,
    acpeak  = res$acpeak
  )

  attr(out_obj, "trackFormats")  <- c("REAL32", "REAL32", "REAL32", "REAL32")
  attr(out_obj, "sampleRate")    <- frame_rate
  attr(out_obj, "origFreq")      <- as.numeric(res$sample_rate)
  attr(out_obj, "startTime")     <- 0.0
  attr(out_obj, "startRecord")   <- 1L
  attr(out_obj, "endRecord")     <- n_frames
  attr(out_obj, "fileInfo")      <- c(20L, 2L)
  class(out_obj) <- "AsspDataObj"
  out_obj
}


##' Convert Snack formant result to AsspDataObj
##'
##' @param res List from snackf_cpp (fm, bw, sample_rate, n_frames)
##' @param windowShift Frame shift in milliseconds
##' @param numFormants Number of formants
##' @return AsspDataObj with fm and bw tracks
##' @keywords internal
create_formant_asspobj <- function(res, windowShift, numFormants) {
  n_frames   <- as.integer(res$n_frames)
  frame_rate <- 1000.0 / windowShift

  out_obj <- list(
    fm = res$fm,
    bw = res$bw
  )

  attr(out_obj, "trackFormats")  <- c("REAL32", "REAL32")
  attr(out_obj, "sampleRate")    <- frame_rate
  attr(out_obj, "origFreq")      <- as.numeric(res$sample_rate)
  attr(out_obj, "startTime")     <- 0.0
  attr(out_obj, "startRecord")   <- 1L
  attr(out_obj, "endRecord")     <- n_frames
  attr(out_obj, "fileInfo")      <- c(20L, 2L)
  class(out_obj) <- "AsspDataObj"
  out_obj
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

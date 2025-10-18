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
  sample_rate <- pitch_result$sample_rate
  n_frames <- pitch_result$n_frames

  # Create AsspDataObj
  out_obj <- list(
    f0 = f0_matrix
  )

  # Set attributes
  attr(out_obj, "sampleRate") <- sample_rate
  attr(out_obj, "startTime") <- 0.0
  attr(out_obj, "startRecord") <- 1L
  attr(out_obj, "endRecord") <- n_frames
  attr(out_obj, "windowShift") <- windowShift / 1000.0  # Convert ms to seconds

  class(out_obj) <- c("AsspDataObj", "list")

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

#' S7 Method System for DSP Functions
#'
#' NOTE: S7 runtime method dispatch is disabled due to incompatibility with
#' existing function signatures (dispatch arguments cannot have default values).
#'
#' Instead, users can:
#' 1. Use AVAudio class for in-memory audio processing
#' 2. Convert AVAudio to temporary files with avaudio_to_tempfile()
#' 3. Pass temp files to existing DSP functions
#'
#' Example:
#' ```r
#' audio <- read_avaudio("speech.wav")
#' temp_file <- avaudio_to_tempfile(audio)
#' result <- trk_rapt(temp_file, toFile = FALSE)
#' unlink(temp_file)
#' ```
#'
#' Future versions may add explicit AVAudio support to individual functions.
#'
#' @name s7-methods
NULL

#' Setup S7 Method Dispatch for DSP Functions
#'
#' Currently disabled - S7 dispatch incompatible with existing function signatures.
#' This function is a placeholder for future S7 integration.
#'
#' @return NULL (called for side effects)
#' @keywords internal
.setup_s7_methods <- function() {
  # Disabled - incompatible with existing function signatures
  # Future: Add explicit AVAudio methods to individual functions
  invisible(NULL)
}

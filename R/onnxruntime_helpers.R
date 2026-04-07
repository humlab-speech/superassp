# onnxruntime_helpers.R — R-level helpers for ONNX Runtime session management
#
# These wrap the C++ Rcpp exports with nicer error messages and
# provide the high-level session API used by trk_crepe, trk_brouhaha, etc.


#' Create an ONNX Runtime inference session
#'
#' @param model_path Path to an ONNX model file (.onnx).
#' @param num_threads Number of intra-op threads. 0 = auto (ORT default).
#'
#' @return An external pointer to the ORT session (class "ort_session").
#'
#' @details
#' The session holds the loaded model and is reusable across multiple
#' inference calls. It is automatically released when garbage collected.
#'
#' @examples
#' \dontrun{
#' sess <- ort_session("model.onnx")
#' result <- ort_run(sess, list(input = rnorm(1024)), list(c(1L, 1024L)))
#' }
#'
#' @keywords internal
ort_session <- function(model_path, num_threads = 0L) {
  .ort_ensure()

  if (!file.exists(model_path)) {
    cli::cli_abort("Model file not found: {model_path}")
  }

  ort_create_session_cpp(normalizePath(model_path, mustWork = TRUE),
                          as.integer(num_threads))
}


#' Run ONNX Runtime inference
#'
#' @param session An ORT session from \code{ort_session()}.
#' @param inputs Named list of input data. Each element is a numeric vector.
#' @param shapes List of integer vectors specifying the shape of each input.
#'   Must be in the same order as \code{inputs}.
#' @param output_names Character vector of output tensor names to fetch.
#'   NULL = fetch all outputs.
#'
#' @return Named list of output tensors (numeric vectors with "shape" attribute).
#'
#' @keywords internal
ort_run <- function(session, inputs, shapes, output_names = NULL) {
  if (!inherits(session, "ort_session")) {
    cli::cli_abort("session must be created by ort_session()")
  }

  input_names <- names(inputs)
  if (is.null(input_names) || any(input_names == "")) {
    cli::cli_abort("All inputs must be named.")
  }

  # Convert shapes to list of integer vectors
  shapes <- lapply(shapes, as.integer)

  ort_run_cpp(session,
              input_names,
              unname(inputs),
              shapes,
              output_names)
}


#' Get ORT session input/output metadata
#'
#' @param session An ORT session from \code{ort_session()}.
#' @return Named list of tensor info (name, shape, type).
#' @keywords internal
ort_input_info <- function(session) {
  if (!inherits(session, "ort_session")) {
    cli::cli_abort("session must be created by ort_session()")
  }
  ort_session_input_info_cpp(session)
}

#' @rdname ort_input_info
#' @keywords internal
ort_output_info <- function(session) {
  if (!inherits(session, "ort_session")) {
    cli::cli_abort("session must be created by ort_session()")
  }
  ort_session_output_info_cpp(session)
}


#' Ensure ONNX Runtime is available, stop if not
#' @keywords internal
.ort_ensure <- function() {
  .ort_set_cached_path()
  if (!ort_available_cpp()) {
    cli::cli_abort(c(
      "ONNX Runtime is not installed.",
      "i" = "Install it with: {.code install_onnxruntime()}",
      "i" = "This is required for CREPE, Brouhaha, and Swift-F0 C++ inference."
    ))
  }
  invisible(TRUE)
}

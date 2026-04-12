##' @keywords internal
.onAttach <- function(libname, pkgname) {
}

##' @keywords internal
.onUnload <- function(libpath) {
  # Release ONNX Runtime environment and shared library before unload
  # to prevent use-after-free in XPtr finalizers
  tryCatch(ort_cleanup_cpp(), error = function(e) NULL)
  library.dynam.unload("superassp", libpath)
}

##' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Fix S3 method dispatch for base generics (print, summary).
  # R's namespace loader sees these names in our namespace (auto-created lazy
  # bindings) and treats them as local generics, so the methods never reach
  # base's S3 methods table.  Re-register them explicitly.
  ns <- asNamespace(pkgname)
  registerS3method("print", "AsspDataObj", ns$print.AsspDataObj, envir = asNamespace("base"))
  registerS3method("print", "JsonTrackObj", ns$print.JsonTrackObj, envir = asNamespace("base"))
  registerS3method("summary", "JsonTrackObj", ns$summary.JsonTrackObj, envir = asNamespace("base"))

  # Setup S7 method dispatch for DSP functions (lst_*, trk_*)
  # This enables AVAudio object support while maintaining backward compatibility
  tryCatch({
    .setup_s7_methods()
  }, error = function(e) {
    warning("Failed to setup S7 method dispatch: ", e$message, call. = FALSE)
  })

  # Register psychoacoustic units (Bark scale, etc.) with units package
  tryCatch({
    .onLoad_psychoacoustic_units()
  }, error = function(e) {
    # Silent - units package may not be available
  })
}

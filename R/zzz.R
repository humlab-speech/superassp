##' @keywords internal
.onAttach <- function(libname, pkgname) {}

# Column names referenced inside dplyr/tidyr verbs, with() and other
# non-standard-evaluation calls. codetools (R CMD check) cannot see these
# bindings, so they are declared here.
utils::globalVariables(c(
  "Hz", "audio", "extension", "frame_time", "i1", "mediaFile",
  "output", "times_norm", "times_orig", "times_rel"
))

##' @keywords internal
.onUnload <- function(libpath) {
  # Release ONNX Runtime environment and shared library before unload
  # to prevent use-after-free in XPtr finalizers
  tryCatch(ort_cleanup_cpp(), error = function(e) NULL)
  library.dynam.unload("superassp", libpath)
}

##' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Default thread count for OpenMP-accelerated DSP kernels. Honoured by
  # any future #pragma omp parallel blocks that read getOption(). Users
  # can override via options(superassp.threads = N). The OMP_NUM_THREADS
  # environment variable always wins at the OpenMP runtime level.
  if (is.null(getOption("superassp.threads"))) {
    options(superassp.threads = max(1L, parallel::detectCores(logical = FALSE) - 1L))
  }

  # Fix S3 method dispatch for base generics (print, summary).
  # R's namespace loader sees these names in our namespace (auto-created lazy
  # bindings) and treats them as local generics, so the methods never reach
  # base's S3 methods table.  Re-register them explicitly.
  ns <- asNamespace(pkgname)
  registerS3method("print", "AsspDataObj", ns$print.AsspDataObj, envir = asNamespace("base"))
  registerS3method("print", "JsonTrackObj", ns$print.JsonTrackObj, envir = asNamespace("base"))
  registerS3method("summary", "JsonTrackObj", ns$summary.JsonTrackObj, envir = asNamespace("base"))

  # S7 classes (AVAudio) register their print/summary/etc methods for base
  # generics via S7::method<- at parse time; those only reach base's S3
  # methods table once S7::methods_register() runs at load time (see
  # ?S7::methods_register). Without this, an installed+library()-loaded
  # package silently falls back to S7's default print.S7_object.
  S7::methods_register()

  # Setup S7 method dispatch for DSP functions (lst_*, trk_*)
  # This enables AVAudio object support while maintaining backward compatibility
  tryCatch({
    .setup_s7_methods()
  }, error = function(e) {
    cli::cli_warn("Failed to setup S7 method dispatch: {e$message}")
  })

  # Register psychoacoustic units (Bark scale, etc.) with units package
  tryCatch({
    .onLoad_psychoacoustic_units()
  }, error = function(e) {
    # Silent - units package may not be available
  })
}

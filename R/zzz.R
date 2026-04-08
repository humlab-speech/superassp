##' @keywords internal
.onAttach <- function(libname, pkgname) {
}

##' @keywords internal
.onLoad <- function(libname, pkgname) {
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

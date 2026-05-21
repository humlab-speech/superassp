#' F0 and voicing detection
#' Computes F0 and voiced/unvoiced decisions. The "srh" method is a MATLAB-
#' faithful port of Drugman & Alwan (2011) Summation of Residual Harmonics.
#' The "praat" method delegates to pladdrr cross-correlation pitch tracking.
#' @param x Speech signal (numeric vector, samples)
#' @param fs Sampling frequency (Hz)
#' @param f0_min Minimum F0 in Hz (default 20)
#' @param f0_max Maximum F0 in Hz (default 500)
#' @param method "srh" (default, MATLAB-faithful) or "praat" (cross-check)
#' @return List with f0, VUV, SRHVal (srh only), time, fs
#' @export
.vat_pitch <- function(x, fs, f0_min = 20, f0_max = 500,
                       method = c("srh", "praat")) {
  method <- match.arg(method)
  if (method == "srh") {
    return(vat_srh_pitch_cpp(x, fs, f0_min, f0_max))
  }
  if (!requireNamespace("pladdrr", quietly = TRUE)) {
    stop("Package 'pladdrr' is required for method='praat'. ",
         "Install it or use method='srh'.")
  }

  # Write to temp WAV and load in pladdrr
  tmp <- tempfile(fileext = ".wav")
  on.exit(unlink(tmp), add = TRUE)

  if (requireNamespace("tuneR", quietly = TRUE)) {
    x_scaled <- x / (max(abs(x)) + 1e-10)
    w <- tuneR::Wave(left = as.integer(x_scaled * 32767),
                     samp.rate = as.integer(fs), bit = 16L)
    tuneR::writeWave(w, tmp)
  } else {
    stop("'tuneR' is required to write temp WAV for pladdrr pitch tracking.")
  }

  snd <- pladdrr::Sound(tmp)
  pitch_obj <- snd$to_pitch_cc(time_step  = 0.01,
                                pitch_floor   = max(f0_min, 50),
                                pitch_ceiling = f0_max)

  n_samples <- length(x)
  n_frames  <- pitch_obj$get_number_of_frames()
  xmin      <- pitch_obj$get_xmin()
  xmax      <- pitch_obj$get_xmax()

  pitch_vals  <- pitch_obj$get_values_vector()  # length n_frames, 0 = unvoiced
  pitch_times <- seq(xmin, xmax, length.out = n_frames)

  # Interpolate to per-sample
  times_samples <- (seq_len(n_samples) - 1) / fs
  if (n_frames >= 2) {
    f0_interp <- stats::approx(pitch_times, pitch_vals,
                               xout  = times_samples,
                               method = "linear",
                               rule   = 2)$y
  } else {
    f0_interp <- rep(0, n_samples)
  }
  f0_interp[is.na(f0_interp)] <- 0

  VUV <- as.integer(f0_interp > f0_min)
  f0_interp[VUV == 0] <- 0

  list(f0 = f0_interp, VUV = VUV)
}

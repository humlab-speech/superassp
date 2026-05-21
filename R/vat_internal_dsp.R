#' DSP primitives (MATLAB-faithful)
#' Wrappers over the Rcpp/Armadillo DSP kernels used internally by the
#' Voice Analysis Toolkit. Exposed for unit testing against MATLAB reference
#' outputs.
#' @name vat_dsp
NULL

#' @rdname vat_dsp
#' @param N window length
#' @param symmetric logical, TRUE for `hamming(N)` style, FALSE for `hamming(N,'periodic')`
#' @export
.vat_hamming <- function(N, symmetric = TRUE) vat_hamming_cpp(as.integer(N), symmetric)

#' @rdname vat_dsp
#' @export
.vat_hanning <- function(N, symmetric = TRUE) vat_hanning_cpp(as.integer(N), symmetric)

#' @rdname vat_dsp
#' @param beta Kaiser beta parameter
#' @export
.vat_kaiser <- function(N, beta) vat_kaiser_cpp(as.integer(N), beta)

#' @rdname vat_dsp
#' @param b numerator coefficients
#' @param a denominator coefficients
#' @param x input signal
#' @param zi initial state (optional)
#' @export
.vat_filter <- function(b, a, x, zi = NULL) vat_filter_cpp(b, a, x, zi)

#' @rdname vat_dsp
#' @export
.vat_filtfilt <- function(b, a, x) vat_filtfilt_cpp(b, a, x)

#' @rdname vat_dsp
#' @param n FIR order (length = n+1)
#' @param Wn cutoff(s) in [0,1] where 1 = Nyquist
#' @param type "low", "high", "bandpass", "stop"
#' @export
.vat_fir1 <- function(n, Wn, type = "low") vat_fir1_cpp(as.integer(n), Wn, type)

#' @rdname vat_dsp
#' @export
.vat_butter <- function(n, Wn, type = "low") vat_butter_cpp(as.integer(n), Wn, type)

#' @rdname vat_dsp
#' @export
.vat_medfilt1 <- function(x, n) vat_medfilt1_cpp(x, as.integer(n))

#' @rdname vat_dsp
#' @param xq query points
#' @param y data values
#' @param method "linear" or "spline"
#' @export
.vat_interp1 <- function(x, y, xq, method = c("linear", "spline")) {
  method <- match.arg(method)
  if (method == "linear") vat_interp1_linear_cpp(x, y, xq)
  else vat_interp1_spline_cpp(x, y, xq)
}

#' @rdname vat_dsp
#' @param min_peak_height minimum height to accept
#' @param min_peak_distance minimum sample distance between peaks
#' @export
.vat_findpeaks <- function(x, min_peak_height = -Inf, min_peak_distance = 1) {
  mph <- if (is.infinite(min_peak_height) && min_peak_height < 0) -1e300 else min_peak_height
  vat_findpeaks_cpp(x, mph, as.integer(min_peak_distance))
}

#' @rdname vat_dsp
#' @param p upsample factor
#' @param q downsample factor
#' @export
.vat_resample <- function(x, p, q, beta = 5.0) vat_resample_cpp(x, as.integer(p), as.integer(q), beta)

#' @rdname vat_dsp
#' @param nfft FFT length (default = length(x))
#' @export
.vat_fft <- function(x, nfft = -1L) vat_fft_cpp(x, as.integer(nfft))

# ---------- LPC ----------

#' LPC primitives
#' @name vat_lpc
NULL

#' @rdname vat_lpc
#' @param s input signal
#' @param p model order
#' @export
.vat_lpcauto <- function(s, p) vat_lpcauto_cpp(s, as.integer(p))

#' @rdname vat_lpc
#' @export
.vat_burg <- function(s, p) vat_burg_cpp(s, as.integer(p))

#' @rdname vat_lpc
#' @param ar AR coefficients (ar(0) = 1)
#' @export
.vat_lpcar2rf <- function(ar) vat_lpcar2rf_cpp(ar)

#' @rdname vat_lpc
#' @export
.vat_lpcar2ra <- function(ar) vat_lpcar2ra_cpp(ar)

#' @rdname vat_lpc
#' @param k reflection coefficients
#' @export
.vat_lpcrf2rr <- function(k) vat_lpcrf2rr_cpp(k)

#' @rdname vat_lpc
#' @param ar1,ar2 AR coefficient vectors
#' @param symmetric symmetric Itakura distance if TRUE
#' @export
.vat_distitar <- function(ar1, ar2, symmetric = FALSE) vat_distitar_cpp(ar1, ar2, symmetric)

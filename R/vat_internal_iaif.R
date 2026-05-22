#' Iterative Adaptive Inverse Filtering (IAIF)
#' Estimates the glottal flow derivative from a speech signal using the IAIF
#' algorithm (Alku et al. 1992). Translates \code{IAIF.m} (John Kane, 2012).
#' @param x Speech signal (numeric vector)
#' @param fs Sampling frequency (Hz)
#' @param GCI Glottal closure instants (integer sample positions). If NULL,
#'   computed via \code{.vat_se_vq()}.
#' @param p LPC order. Default: \code{round(fs/1000) + 2}.
#' @return List with:
#'   \describe{
#'     \item{g_iaif}{Glottal flow derivative estimate (same length as x)}
#'     \item{ar_lpc}{LPC coefficient matrix (p+1 x nGCI)}
#'     \item{e_lpc}{Residual energy vector (nGCI)}
#'   }
#' @references \insertCite{Alku1992IAIF}{superassp}
#' @param backend "cpp" (default, full Rcpp pipeline) or "r" (legacy).
#' @keywords internal
#' @noRd
.vat_iaif <- function(x, fs, GCI = NULL, p = NULL, backend = c("cpp", "r")) {
  backend <- match.arg(backend)
  x <- as.numeric(x)
  if (is.null(p)) p <- round(fs / 1000) + 2

  if (is.null(GCI)) {
    se_res <- .vat_se_vq(x, fs)
    GCI <- se_res$GCI
  }
  GCI <- as.integer(GCI)

  if (backend == "cpp") {
    res <- vat_iaif_full_cpp(x, fs, GCI, as.integer(p))
    return(list(g       = res$g,
                g_iaif  = res$dg,
                dg      = res$dg,
                ar_lpc  = res$ar_lpc,
                e_lpc   = res$e_lpc))
  }

  # --- Part 1: High-pass FIR filter (50 Hz) to remove DC ---
  nc  <- 704L
  fc  <- 50
  nfc <- fc / (fs / 2)
  fir_hp <- signal::fir1(nc, nfc, type = "high")  # returns Ma object
  x_hp   <- as.numeric(signal::filter(fir_hp, x))
  # Linear-phase delay compensation: advance by nc/2 samples
  half_nc <- nc %/% 2L
  n_total  <- length(x)
  x_filt <- c(x_hp[(half_nc + 1L):n_total], x[(n_total - half_nc + 1L):n_total])

  # --- Parts 2-3: Pre-emphasis (LPC order 1) ---
  res1 <- calc_residual_r(x_filt, x_filt, 1L, GCI)
  x_emph <- res1$vector_res

  # --- Parts 4-5: First glottal estimate (LPC order p) ---
  res2 <- calc_residual_r(x_filt, x_emph, p, GCI)
  ug1  <- res2$vector_res   # first glottal flow derivative

  # --- Parts 7-8: Remove glottal source effect (LPC order 4) ---
  res3 <- calc_residual_r(x_filt, ug1, 4L, GCI)
  vt_signal <- res3$vector_res

  # --- Parts 9-10: Second glottal estimate (LPC order p) ---
  res4 <- calc_residual_r(x_filt, vt_signal, p, GCI)

  list(g_iaif  = res4$vector_res,
       ar_lpc  = res4$ar_lpc,
       e_lpc   = res4$e_lpc)
}

# ── Internal: per-GCI LPC + inverse filter ──────────────────────────────────

#' Per-GCI LPC analysis and inverse filtering
#' R wrapper around the Rcpp implementation of \code{calc_residual.m}.
#' If GCI vector is too short, falls back to pure R.
#' @keywords internal
calc_residual_r <- function(x, x_lpc, ord_lpc, GCI) {
  vat_calc_residual_cpp(as.numeric(x), as.numeric(x_lpc),
                    as.integer(ord_lpc), as.integer(GCI))
}

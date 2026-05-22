#' Maxima Dispersion Quotient (MDQ)
#' Computes the Maxima Dispersion Quotient (Kane & Gobl 2013) from the LP
#' residual signal. Useful for breathy-to-tense voice discrimination.
#' Translates \code{MDQ/get_MDQ.m} (recovered from git commit bb9b314).
#' @param res LP residual signal (e.g. from \code{.vat_se_vq()$res})
#' @param fs Sampling frequency (Hz)
#' @param GCI Glottal closure instants (integer sample positions)
#' @return Numeric vector of MDQ values (one per GCI).
#' @references \insertCite{Kane2013MDQ}{superassp}
#' @keywords internal
#' @noRd
.vat_mdq <- function(res, fs, GCI) {
  res <- as.numeric(res)
  GCI <- as.integer(GCI)
  vat_mdq_cpp(res, fs, GCI)
}

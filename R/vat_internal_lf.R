#' LF model: Newton-Raphson area-balance solver
#' Direct port of \code{lf_Area_newton.c} (Kane 2012). Solves for the LF
#' model parameters alpha (open-phase decay) and epsilon (return-phase
#' decay) given Tc, fs, Tp, Te, Ta, EE.
#' @param Tc closure time
#' @param fs sampling frequency
#' @param Tp time of peak glottal flow
#' @param Te time of negative excitation
#' @param Ta return-phase time constant
#' @param EE excitation strength
#' @return list with `alpha` and `epsi`
#' @export
.vat_lf_area_newton <- function(Tc, fs, Tp, Te, Ta, EE) {
  vat_lf_area_newton_cpp(Tc, fs, Tp, Te, Ta, EE)
}

#' Convert Rd parameter to LF Ra/Rk/Rg.
#' @param Rd shape descriptor
#' @param EE excitation strength
#' @param F0 fundamental frequency (Hz)
#' @return list with Ra, Rk, Rg
#' @export
.vat_rd2r <- function(Rd, EE, F0) vat_rd2r_cpp(Rd, EE, F0)

#' Generate one LF model pulse.
#' @param F0 fundamental frequency (Hz)
#' @param fs sampling frequency (Hz)
#' @param Ra,Rk,Rg LF R-parameters
#' @param EE excitation strength
#' @return numeric vector (one glottal-flow-derivative pulse)
#' @export
.vat_lf_cont <- function(F0, fs, Ra, Rk, Rg, EE = 1.0) {
  vat_lf_cont_cpp(F0, fs, Ra, Rk, Rg, EE)
}

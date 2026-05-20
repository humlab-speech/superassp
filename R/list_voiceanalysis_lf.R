##' Synthesise an LF model glottal pulse via voiceanalysis
##'
##' Convenience wrapper around \code{voiceanalysis::vat_rd2r()} and
##' \code{voiceanalysis::vat_lf_cont()} for generating Liljencrants-Fant
##' \insertCite{Fant1985LF}{superassp} glottal flow derivative pulses from
##' a Rd shape descriptor. Useful for source-modelling experiments and
##' synthesis demos.
##'
##' This is a pure synthesis utility — it does not analyse an input file.
##' For per-cycle Rd estimation on real speech, a future
##' \code{lst_lf_vat_fit()} will land once \code{dyProg_LF} is ported.
##'
##' @param Rd Numeric. LF shape descriptor (typical range 0.3–2.5;
##'   smaller = tenser).
##' @param F0 Numeric. Fundamental frequency in Hz.
##' @param fs Numeric. Sampling frequency in Hz.
##' @param EE Numeric. Excitation strength (default 1.0).
##'
##' @return Named list:
##'   \describe{
##'     \item{\code{pulse}}{Numeric vector — LF glottal flow derivative
##'       sampled at \code{fs}.}
##'     \item{\code{Rd}}{The Rd that was used.}
##'     \item{\code{Ra}, \code{Rk}, \code{Rg}}{Derived LF R-parameters.}
##'     \item{\code{F0}, \code{fs}, \code{EE}}{Inputs (for reproducibility).}
##'   }
##'
##' @references
##' \insertCite{Fant1985LF}{superassp}
##' @export
lst_lf_vat_synthesis <- function(Rd, F0, fs, EE = 1.0) {
  if (!requireNamespace("voiceanalysis", quietly = TRUE))
    cli::cli_abort(c("Package {.pkg voiceanalysis} is required.",
                     "i" = "Install via {.code pak::pkg_install('jckane/Voice_Analysis_Toolkit/voiceanalysis')}"))
  r <- voiceanalysis::vat_rd2r(Rd, EE, F0)
  pulse <- voiceanalysis::vat_lf_cont(F0, fs, r$Ra, r$Rk, r$Rg, EE)
  list(pulse = as.numeric(pulse),
       Rd = Rd, Ra = r$Ra, Rk = r$Rk, Rg = r$Rg,
       F0 = F0, fs = fs, EE = EE)
}

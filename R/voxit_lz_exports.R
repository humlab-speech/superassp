#' Lempel-Ziv complexity
#'
#' Lempel-Ziv (1976) complexity of a binary sequence, as computed by Voxit
#' (`100 * calc_lz_complexity(S, 'exhaustive', 1)`).
#'
#' @param S Logical or integer vector (non-zero / NA count as TRUE)
#' @param type Only "exhaustive" (LZ76 exhaustive production history)
#' @param normalize If TRUE, return `100 * C / (n / log2(n))` as Voxit does
#'
#' @return Numeric complexity score
#' @keywords internal
#' @noRd
lz_complexity_cpp <- function(S, type = "exhaustive", normalize = TRUE) {
  .Call(`_superassp_lz_complexity_cpp`, as.logical(S), type, normalize)
}

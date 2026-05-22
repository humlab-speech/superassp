#' Lempel-Ziv complexity
#'
#' Compute binary sequence complexity via Lempel-Ziv factorization.
#'
#' @param S Logical or integer vector (converted to binary)
#' @param type Complexity measure: "exhaustive" (LZ76) or "primitive"
#' @param normalize If TRUE, normalize to 0–100 range
#'
#' @return Numeric complexity score
#' @export
lz_complexity_cpp <- function(S, type = "exhaustive", normalize = TRUE) {
  .Call(`_superassp_lz_complexity_cpp`, as.logical(S), type, normalize)
}

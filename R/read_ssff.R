#' Read an SSFF or audio file into an AsspDataObj
#'
#' A user-facing wrapper around the internal ASSP C-level reader.
#' Interface is identical to the legacy \code{read.AsspDataObj}.
#'
#' @param fname Path to an SSFF or native ASSP audio file (WAV, AU, NIST, etc.).
#' @param begin Start of region to read (seconds, or samples if \code{samples=TRUE}). Default 0 = file start.
#' @param end   End of region to read (seconds, or samples if \code{samples=TRUE}). Default 0 = file end.
#' @param samples Logical. If \code{TRUE}, \code{begin}/\code{end} are in samples; otherwise in seconds.
#' @return An \code{AsspDataObj}.
#' @seealso \code{\link{read_audio}} for universal format support including MP3/MP4.
#' @export
read_ssff <- function(fname, begin = 0, end = 0, samples = FALSE) {
  fname <- prepareFiles(fname)
  if (inherits(begin, "integer")) begin <- as.numeric(begin)
  if (inherits(end, "integer"))   end   <- as.numeric(end)
  .External("getDObj2", fname, begin = begin, end = end, samples = samples,
            PACKAGE = "superassp")
}

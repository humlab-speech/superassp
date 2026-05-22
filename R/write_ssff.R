#' Write an AsspDataObj to an SSFF file
#'
#' User-facing wrapper around the ASSP C-level writer.
#' Interface mirrors the legacy \code{write.AsspDataObj}.
#'
#' @param dobj An \code{AsspDataObj}.
#' @param file Output file path. Defaults to the \code{filePath} attribute of \code{dobj}.
#' @return Invisibly, the resolved output file path (after \code{path.expand}).
#'   Use this for chained pipelines. The file written is an SSFF container with
#'   the tracks stored in \code{dobj} encoded according to
#'   \code{attr(dobj, "trackFormats")}.
#' @seealso \code{\link{write_jstf}} for JSTF (JSON) output.
#' @examples
#' \dontrun{
#' wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
#' f0  <- trk_pitch_rapt(wav, toFile = FALSE)
#'
#' out <- tempfile(fileext = ".f0")
#' write_ssff(f0, file = out)
#'
#' f0_back <- read_ssff(out)
#' identical(names(f0), names(f0_back))
#' }
#' @export
write_ssff <- function(dobj, file = attr(dobj, "filePath")) {
  if (is.null(file))
    stop("File path not set internally. Please specify!")
  file <- path.expand(file)
  .Call("writeDObj_", dobj, file, PACKAGE = "superassp")
  invisible(file)
}

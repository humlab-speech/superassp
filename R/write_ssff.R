#' Write an AsspDataObj to an SSFF file
#'
#' User-facing wrapper around the ASSP C-level writer.
#' Interface mirrors the legacy \code{write.AsspDataObj}.
#'
#' @param dobj An \code{AsspDataObj}.
#' @param file Output file path. Defaults to the \code{filePath} attribute of \code{dobj}.
#' @return Invisibly, \code{file}.
#' @seealso \code{\link{write_json_track}} for JSTF (JSON) output.
#' @export
write_ssff <- function(dobj, file = attr(dobj, "filePath")) {
  if (is.null(file))
    stop("File path not set internally. Please specify!")
  file <- path.expand(file)
  .Call("writeDObj_", dobj, file, PACKAGE = "superassp")
  invisible(file)
}

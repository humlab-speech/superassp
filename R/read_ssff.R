#' Read an SSFF or audio file into an AsspDataObj
#'
#' A user-facing wrapper around the internal ASSP C-level reader.
#' Interface is identical to the legacy \code{read.AsspDataObj}.
#'
#' @param fname Path to an SSFF or native ASSP audio file (WAV, AU, NIST, etc.).
#' @param begin Start of region to read (seconds, or samples if \code{samples=TRUE}). Default 0 = file start.
#' @param end   End of region to read (seconds, or samples if \code{samples=TRUE}). Default 0 = file end.
#' @param samples Logical. If \code{TRUE}, \code{begin}/\code{end} are in samples; otherwise in seconds.
#' @return An \code{AsspDataObj}. For audio files, contains an \code{audio}
#'   track (n_samples x n_channels). For SSFF tracks, contains one matrix per
#'   stored track (e.g. \code{F0}, \code{fm}, \code{bw}, \code{rms}) at the
#'   analysis frame rate. Standard attributes include \code{sampleRate},
#'   \code{startTime}, \code{startRecord}, \code{endRecord}, \code{trackFormats}
#'   and \code{filePath}.
#' @seealso \code{\link{read_audio}} for universal format support including MP3/MP4.
#' @examples
#' \dontrun{
#' # Read an audio file from the bundled samples
#' wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
#' au  <- read_ssff(wav)
#' names(au)               # "audio"
#' attr(au, "sampleRate")  # native sample rate
#'
#' # Read an SSFF parameter track produced earlier by a trk_* function
#' f0_path <- tempfile(fileext = ".f0")
#' trk_pitch_rapt(wav, toFile = TRUE, outputDirectory = dirname(f0_path),
#'                explicitExt = "f0")
#' f0_obj <- read_ssff(file.path(dirname(f0_path),
#'                               paste0(tools::file_path_sans_ext(basename(wav)), ".f0")))
#' names(f0_obj)
#' }
#' @export
read_ssff <- function(fname, begin = 0, end = 0, samples = FALSE) {
  fname <- prepareFiles(fname)
  if (inherits(begin, "integer")) begin <- as.numeric(begin)
  if (inherits(end, "integer"))   end   <- as.numeric(end)
  .External("getDObj2", fname, begin = begin, end = end, samples = samples,
            PACKAGE = "superassp")
}

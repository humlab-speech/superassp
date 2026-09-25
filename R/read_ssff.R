#' Peek at an SSFF file's header (dataRate/startRecord/numRecords) cheaply
#'
#' Opens the file and reads only its header (no records), so callers can
#' learn the analysis-frame rate without paying for a full read.
#' @keywords internal
#' @noRd
ssff_header_peek <- function(fname) {
  .Call("getSSFFHeader", fname, PACKAGE = "superassp")
}

#' Round a time point to the nearest analysis-frame boundary
#' @keywords internal
#' @noRd
snap_to_nearest_frame <- function(t, dataRate) {
  round(t * dataRate) / dataRate
}

#' Read an SSFF or audio file into an AsspDataObj
#'
#' A user-facing wrapper around the internal ASSP C-level reader.
#' Interface is identical to the legacy \code{read.AsspDataObj}.
#'
#' @param fname Path to an SSFF or native ASSP audio file (WAV, AU, NIST, etc.).
#' @param begin Start of region to read (seconds, or samples if \code{samples=TRUE}). Default 0 = file start.
#' @param end   End of region to read (seconds, or samples if \code{samples=TRUE}). Default 0 = file end.
#' @param samples Logical. If \code{TRUE}, \code{begin}/\code{end} are in samples; otherwise in seconds.
#' @param zero_to_na Logical. If \code{TRUE}, stored values that are exactly
#'   \code{0} are returned as \code{NA} for every track that is not sampled
#'   audio (SSFF has no NA encoding; \code{0} is its substitute). Default
#'   \code{FALSE}, i.e. the values exactly as stored. \code{\link{read_track}}
#'   passes \code{TRUE} for SSFF files.
#' @param tracks Optional character vector of track names to read. \code{NULL}
#'   (default) reads every track; other tracks are skipped without being
#'   converted, which is considerably faster for files that store more than
#'   one track.
#' @param threads Number of threads used to convert large files (default 1,
#'   serial). Values above 1 need a build with OpenMP support; the results are
#'   identical either way.
#' @param snap One of \code{"none"} (default) or \code{"nearest"}. A single
#'   time point (\code{begin == end}, both non-zero) that does not fall
#'   exactly on an analysis-frame boundary errors by default \emph{— this
#'   matches \code{wrassp::read.AsspDataObj()} exactly, since faithfulness to
#'   the reference implementation is the priority}. Pass \code{"nearest"} to
#'   instead round such a request to the nearest frame and return it. No
#'   effect when \code{samples = TRUE}, when reading a range
#'   (\code{begin != end}), or for \code{begin == end == 0} (whole file).
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
read_ssff <- function(fname, begin = 0, end = 0, samples = FALSE,
                      zero_to_na = FALSE, tracks = NULL, threads = 1L,
                      snap = c("none", "nearest")) {
  fname <- prepareFiles(fname)
  if (inherits(begin, "integer")) begin <- as.numeric(begin)
  if (inherits(end, "integer"))   end   <- as.numeric(end)
  if (!is.null(tracks) && !is.character(tracks)) {
    cli::cli_abort("{.arg tracks} must be a character vector of track names, or {.code NULL}.")
  }
  if (!is.logical(zero_to_na) || length(zero_to_na) != 1L || is.na(zero_to_na)) {
    cli::cli_abort("{.arg zero_to_na} must be {.code TRUE} or {.code FALSE}.")
  }
  threads <- suppressWarnings(as.integer(threads)[1])
  if (is.na(threads) || threads < 1L) threads <- 1L
  snap <- match.arg(snap)

  single_time_point <- !isTRUE(samples) && begin == end && begin > 0

  if (snap == "nearest" && single_time_point) {
    header <- ssff_header_peek(fname)
    t_snapped <- snap_to_nearest_frame(begin, header[["dataRate"]])
    begin <- t_snapped
    end   <- t_snapped
  }

  tryCatch(
    .External("getDObj2", fname, begin = begin, end = end, samples = samples,
              zero_to_na = zero_to_na, tracks = tracks, threads = threads,
              PACKAGE = "superassp"),
    error = function(e) {
      if (snap == "none" && single_time_point) {
        cli::cli_abort(c(
          "Cannot read a single time point at {.val {begin}} s: it does not fall on an exact analysis-frame boundary.",
          "i" = "Pass {.code snap = \"nearest\"} to round to the nearest frame, or supply a {.arg begin}/{.arg end} pair spanning at least one frame."
        ), parent = e, call = rlang::caller_env())
      } else {
        stop(e)
      }
    }
  )
}

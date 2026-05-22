#' Memoized media-info lookup
#'
#' Internal helper that wraps [av::av_media_info()] with a per-session cache
#' keyed by canonical path, file size and modification time. Repeated calls to
#' the same unchanged file return from the cache without spawning ffprobe.
#'
#' The cache lives in a package-private environment created at load time
#' (see [.media_info_cache]). It is automatically invalidated when the file's
#' size or mtime changes between calls; manual invalidation is also available
#' via [media_info_cache_clear()].
#'
#' @param fname A character path to a media file.
#' @param force If `TRUE`, bypass any cached entry and re-probe the file.
#' @return The list returned by [av::av_media_info()] (with `audio`, `video`,
#'   `duration`, `bit_rate`, etc.).
#' @keywords internal
#' @noRd
media_info <- function(fname, force = FALSE) {
  if (length(fname) != 1L || !is.character(fname)) {
    stop("media_info() expects a single character path.")
  }
  cache <- .media_info_cache
  if (!file.exists(fname)) {
    return(av::av_media_info(fname))
  }
  fi  <- file.info(fname, extra_cols = FALSE)
  key <- normalizePath(fname, winslash = "/", mustWork = FALSE)
  hit <- cache[[key]]
  if (!force && !is.null(hit) &&
      isTRUE(hit$size == fi$size) &&
      isTRUE(hit$mtime == fi$mtime)) {
    return(hit$info)
  }
  info <- av::av_media_info(fname)
  cache[[key]] <- list(size = fi$size, mtime = fi$mtime, info = info)
  info
}

#' Clear the media-info cache (internal).
#' @keywords internal
#' @noRd
media_info_cache_clear <- function() {
  rm(list = ls(.media_info_cache, all.names = TRUE), envir = .media_info_cache)
  invisible(NULL)
}

# Package-private cache environment. Created at load time via .onLoad in zzz.R
# (or lazily here if zzz.R has not initialised it yet).
.media_info_cache <- new.env(parent = emptyenv())

# hf_model_cache.R — generic download-and-cache helper for ONNX models hosted
# on the Hugging Face Hub, via the huggingfaceR package (Suggests, guarded).
#
# Mirrors the caching convention already used for the ONNX Runtime binary
# itself (see onnxruntime.R): check the R user cache directory first, only
# hit the network on a genuine cache miss, one-time cli::cli_inform() notice.
#
# huggingfaceR::hf_hub_download() has no persistent cache of its own (a NULL
# `dest` writes to a tempfile that vanishes with the session), so that part
# is implemented here instead.


#' Get a cached path to a Hugging Face Hub model file, downloading if needed
#'
#' Downloads \code{filename} from the Hugging Face Hub repo \code{repo_id} on
#' first use and caches it under the R user cache directory
#' (\code{tools::R_user_dir("superassp", "cache")}); subsequent calls return
#' the cached path without any network access.
#'
#' @param repo_id Hugging Face Hub repo id, e.g. \code{"username/model-name"}.
#' @param filename Path to the file inside the repo, e.g. \code{"model.onnx"}.
#' @param subdir Cache subdirectory name for this model (keeps different
#'   models from colliding when their filenames share a basename), e.g.
#'   \code{"swift-f0"}.
#' @param revision Git revision (tag or commit SHA) to pin. Required, and
#'   deliberately has no default — pinning to \code{"main"} would let the
#'   model change under the package without a version bump, undermining
#'   DSP output faithfulness.
#' @param repo_type One of \code{"model"}, \code{"dataset"}, \code{"space"}.
#'   Default \code{"model"}.
#'
#' @return Character path to the cached local file.
#' @keywords internal
.hf_get_cached_model <- function(repo_id, filename, subdir, revision,
                                 repo_type = "model") {

  if (missing(revision) || !nzchar(revision)) {
    cli::cli_abort(c(
      "x" = "{.arg revision} must be a pinned tag or commit SHA.",
      "i" = "Do not pass {.val main} — it would let the model change under superassp without a version bump."
    ))
  }

  cache_dir  <- file.path(tools::R_user_dir("superassp", "cache"), "onnx", subdir)
  model_path <- file.path(cache_dir, basename(filename))

  if (file.exists(model_path)) {
    return(model_path)
  }

  if (!requireNamespace("huggingfaceR", quietly = TRUE)) {
    cli::cli_abort(c(
      "x" = "Package {.pkg huggingfaceR} is required to download the {.val {subdir}} model.",
      "i" = "Install with {.code install.packages(\"huggingfaceR\")}."
    ))
  }

  cli::cli_inform(c(
    "i" = "Downloading {.val {subdir}} model from Hugging Face Hub ({.val {repo_id}}).",
    "i" = "This is a one-time download, cached at {.path {cache_dir}}."
  ))

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

  tryCatch(
    huggingfaceR::hf_hub_download(
      repo_id   = repo_id,
      filename  = filename,
      repo_type = repo_type,
      revision  = revision,
      dest      = cache_dir
    ),
    error = function(e) {
      cli::cli_abort(c(
        "x" = "Failed to download {.val {subdir}} model from Hugging Face Hub: {conditionMessage(e)}",
        "i" = "Check your internet connection and try again."
      ))
    }
  )

  if (!file.exists(model_path)) {
    cli::cli_abort("Download reported success but {.path {model_path}} is missing.")
  }

  model_path
}

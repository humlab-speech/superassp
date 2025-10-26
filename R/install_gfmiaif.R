#' Install Python dependencies for GFM-IAIF
#'
#' This function installs the required Python packages for the GFM-IAIF
#' source-filter separation algorithm. It installs numpy and scipy (required),
#' and optionally numba for 5-10x speedup.
#'
#' @param install_numba If TRUE (default), also install numba for JIT compilation.
#'   This provides 5-10x speedup for the Levinson-Durbin algorithm but adds ~50MB
#'   to the installation size.
#' @param method Installation method passed to \code{reticulate::py_install()}.
#'   Default is "auto".
#' @param conda Path to conda executable (if using conda). Default is "auto".
#'
#' @details
#' \strong{Requirements:}
#' \itemize{
#'   \item \strong{numpy} (>= 1.19): Numerical computing
#'   \item \strong{scipy} (>= 1.5): Scientific computing (signal processing, FFT)
#'   \item \strong{numba} (>= 0.54, optional): JIT compilation for performance
#' }
#'
#' \strong{Installation size:}
#' \itemize{
#'   \item Without numba: ~50-100 MB
#'   \item With numba: ~100-150 MB
#' }
#'
#' \strong{Performance:}
#' \itemize{
#'   \item Without numba: ~1.2 ms per frame (512 samples)
#'   \item With numba: ~0.25 ms per frame (4-5x faster)
#' }
#'
#' The function will use pip by default. If you prefer conda, set the conda
#' parameter or use reticulate::use_condaenv() before calling this function.
#'
#' @return Invisible NULL. Called for side effects (installing packages).
#'
#' @examples
#' \dontrun{
#' # Basic installation (with numba for best performance)
#' install_gfmiaif()
#'
#' # Install without numba (smaller installation, slower processing)
#' install_gfmiaif(install_numba = FALSE)
#'
#' # Use conda instead of pip
#' install_gfmiaif(method = "conda")
#' }
#'
#' @seealso \code{\link{trk_gfmiaif}} for using GFM-IAIF after installation
#'
#' @export
install_gfmiaif <- function(install_numba = TRUE,
                           method = "auto",
                           conda = "auto") {

  message("Installing Python dependencies for GFM-IAIF...")
  message("This may take a few minutes.\n")

  # Core dependencies
  packages <- c("numpy>=1.19", "scipy>=1.5")

  # Add numba if requested
  if (install_numba) {
    packages <- c(packages, "numba>=0.54")
    message("Installing: numpy, scipy, numba (with JIT compilation)")
    message("Installation size: ~100-150 MB")
  } else {
    message("Installing: numpy, scipy (without numba)")
    message("Installation size: ~50-100 MB")
    message("Note: Install numba later for 5-10x speedup")
  }

  # Install packages
  tryCatch({
    reticulate::py_install(
      packages = packages,
      method = method,
      conda = conda,
      pip = TRUE
    )

    message("\n✓ GFM-IAIF dependencies installed successfully!")

    # Check if numba is available
    if (install_numba && reticulate::py_module_available("numba")) {
      message("✓ Numba JIT compilation available (maximum performance)")
    } else if (install_numba) {
      warning("Numba installation may have failed. GFM-IAIF will work but may be slower.",
              call. = FALSE)
    } else {
      message("ℹ Running without numba (install with install_numba=TRUE for 5-10x speedup)")
    }

    message("\nYou can now use trk_gfmiaif() to process audio files.")

  }, error = function(e) {
    stop(
      "Failed to install GFM-IAIF dependencies.\n\n",
      "Error: ", conditionMessage(e), "\n\n",
      "Try installing manually:\n",
      "  reticulate::py_install(c('numpy', 'scipy', 'numba'))\n",
      call. = FALSE
    )
  })

  invisible(NULL)
}

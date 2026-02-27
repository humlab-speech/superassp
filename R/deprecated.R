#' Deprecated EGG Functions - Migrated to eggstract Package
#'
#' @description
#' These functions have been migrated to the \strong{eggstract} package as part
#' of a package restructuring effort. The eggstract package provides a unified
#' API for electroglottographic (EGG) signal analysis with support for multiple
#' output formats (dataframe, SSFF, suggestion).
#'
#' **Migration Timeline:**
#' - These wrapper functions will remain available until superassp v0.10.0
#' - Expected removal: 6-12 months from 2025-10-30
#' - Users should update their code to use eggstract functions directly
#'
#' @name deprecated-egg-functions
#' @keywords internal
NULL

#' @rdname deprecated-egg-functions
#' @export
trk_egg_f0_deprecated <- function(...) {
  .Deprecated(
    new = "eggstract::egg_f0",
    package = "superassp",
    msg = paste(
      "trk_egg_f0() has moved to the eggstract package.",
      "Use: eggstract::egg_f0(..., output_format = 'ssff')",
      "\nThis wrapper will be removed in superassp v0.10.0."
    )
  )

  if (!requireNamespace("eggstract", quietly = TRUE)) {
    cli::cli_abort(c(
      "trk_egg_f0() has moved to the eggstract package",
      "x" = "eggstract package not installed",
      "i" = "Install with: {.code remotes::install_github('humlab-speech/eggstract')}"
    ))
  }

  # Forward to eggstract with SSFF output format
  eggstract::egg_f0(..., output_format = "ssff")
}

#' Helper Function: EGG Migration Information
#'
#' @description
#' Displays information about the migration of EGG functions from superassp to eggstract.
#'
egg_migration_info_superassp <- function() {
  cli::cli_h1("EGG Function Migration to eggstract")

  cli::cli_alert_info("The following function has moved to the eggstract package:")
  cli::cli_ul(c(
    "trk_egg_f0() → eggstract::egg_f0(..., output_format = 'ssff')"
  ))

  cli::cli_h2("Why the migration?")
  cli::cli_ul(c(
    "Unified API: All EGG functions now have consistent parameters",
    "Multiple workflows: Support for dataframe, SSFF, and Suggestion outputs",
    "Focused packages: EGG analysis consolidated in one dedicated package",
    "Better maintenance: Centralized EGG documentation and development"
  ))

  cli::cli_h2("How to migrate your code")
  cli::cli_text("Old code:")
  cli::cli_code("result <- trk_egg_f0(files, toFile = FALSE)")

  cli::cli_text("\nNew code:")
  cli::cli_code("result <- eggstract::egg_f0(files, output_format = 'ssff')")

  cli::cli_h2("Installation")
  cli::cli_code("remotes::install_github('humlab-speech/eggstract')")

  # Check if eggstract is installed
  if (requireNamespace("eggstract", quietly = TRUE)) {
    version <- utils::packageVersion("eggstract")
    cli::cli_alert_success("eggstract v{version} is installed")
  } else {
    cli::cli_alert_warning("eggstract is not yet installed")
  }

  cli::cli_h2("Timeline")
  cli::cli_ul(c(
    "Now: Deprecated wrappers available (with warnings)",
    "superassp v0.10.0: Wrappers will be removed (6-12 months)"
  ))

  cli::cli_text("\nFor more information, see: {.url https://github.com/humlab-speech/eggstract}")

  invisible(NULL)
}

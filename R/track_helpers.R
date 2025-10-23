# Track Name Helper Functions
#
# Internal helper functions for track name manipulation, template expansion,
# and unit assignment following Titze 2015 / Nylén 2024 standards.

#' Detect if track name has placeholder 'i'
#'
#' Checks if a track name template contains the placeholder 'i' using the
#' uniform pattern: uppercase letter + 'i' + (bracket or end-of-string).
#'
#' @param name Character. Track name template to check.
#' @return Logical. TRUE if name contains placeholder, FALSE otherwise.
#'
#' @details
#' The pattern `[A-Z]i(\\[|$)` matches:
#' - Fi\[Hz\] → TRUE (formant frequency template)
#' - Bi\[Hz\] → TRUE (bandwidth template)
#' - LPCi → TRUE (LP coefficient template)
#' - ARFi → TRUE (ARF coefficient template)
#' - Hi\[dB\] → TRUE (harmonic template)
#' - Ai\[dB\] → TRUE (amplitude template)
#'
#' But excludes:
#' - foi\[Hz\] → FALSE (lowercase before i)
#' - intensity → FALSE (lowercase, i not at end)
#' - pitch\[Hz\] → FALSE (no i before bracket)
#' - gain\[dB\] → FALSE (no i)
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' .has_placeholder("Fi[Hz]")     # TRUE
#' .has_placeholder("LPCi")       # TRUE
#' .has_placeholder("fo[Hz]")     # FALSE
#' .has_placeholder("intensity")  # FALSE
#' }
.has_placeholder <- function(name) {
  # Pattern: uppercase letter + 'i' + (bracket or end-of-string)
  grepl("[A-Z]i(\\[|$)", name)
}


#' Expand track template to column names
#'
#' Expands a track name template containing placeholder 'i' to a vector of
#' numbered column names.
#'
#' @param template Character. Track name template (e.g., "Fi\[Hz\]", "LPCi").
#' @param n_cols Integer. Number of columns to generate.
#' @return Character vector of expanded column names.
#'
#' @details
#' Uses the uniform placeholder pattern where the last 'i' after an uppercase
#' letter and before '[' or end-of-string is replaced with sequential numbers.
#'
#' For templates without placeholders, falls back to appending "_1", "_2", etc.
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' .expand_track_template("Fi[Hz]", 4)
#' # [1] "F1[Hz]" "F2[Hz]" "F3[Hz]" "F4[Hz]"
#'
#' .expand_track_template("LPCi", 12)
#' # [1] "LPC1" "LPC2" "LPC3" ... "LPC12"
#'
#' .expand_track_template("Hi[dB]", 3)
#' # [1] "H1[dB]" "H2[dB]" "H3[dB]"
#'
#' # Non-template (fallback)
#' .expand_track_template("unknown", 3)
#' # [1] "unknown_1" "unknown_2" "unknown_3"
#' }
.expand_track_template <- function(template, n_cols) {

  # Check if template has placeholder
  if (!.has_placeholder(template)) {
    # No placeholder - use fallback numbering
    return(paste0(template, "_", seq_len(n_cols)))
  }

  # Expand by substituting last 'i' before [ or end with numbers
  col_names <- character(n_cols)

  for (j in seq_len(n_cols)) {
    # Replace: (uppercase)i(bracket or end) → (uppercase)(number)(bracket or end)
    col_names[j] <- sub("([A-Z])i(\\[|$)", paste0("\\1", j, "\\2"), template)
  }

  col_names
}


#' Clean track names for R data frames
#'
#' Converts SSFF-style track names (with brackets and hyphens) to R-friendly
#' column names (with underscores).
#'
#' @param names Character vector. Track names to clean.
#' @return Character vector of cleaned names.
#'
#' @details
#' Applies the following transformations:
#' 1. `\[Hz\]` → `_Hz` (brackets to underscores)
#' 2. `H1-H2` → `H1_H2` (hyphens to underscores)
#' 3. `(local)` → `_local` (parentheses to underscores)
#' 4. Multiple spaces → single underscore
#'
#' This follows the hybrid naming strategy where SSFF files use scientific
#' notation with brackets, but R data frames use clean underscore notation.
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' .clean_track_names(c("fo[Hz]", "F1[Hz]", "H1-H2c[dB]"))
#' # [1] "fo_Hz" "F1_Hz" "H1_H2c_dB"
#'
#' .clean_track_names(c("Jitter (local)[%]", "Shimmer (dda)[%]"))
#' # [1] "Jitter_local_pct" "Shimmer_dda_pct"
#' }
.clean_track_names <- function(names) {
  # Convert brackets to underscores: [Hz] → _Hz
  names <- gsub("\\[([^]]+)\\]", "_\\1", names)

  # Convert hyphens to underscores: H1-H2 → H1_H2
  names <- gsub("-", "_", names)

  # Convert parentheses to underscores: (local) → _local
  names <- gsub("\\s*\\(([^)]+)\\)", "_\\1", names)

  # Remove multiple spaces, replace with underscore
  names <- gsub("\\s+", "_", names)

  names
}


#' Parse unit from column name
#'
#' Extracts the unit suffix from a cleaned column name.
#'
#' @param col_name Character. Column name to parse.
#' @return Character. Unit name (e.g., "Hz", "dB") or NA if no unit found.
#'
#' @details
#' Expects cleaned column names in format: `param_unit`
#' - `fo_Hz` → "Hz"
#' - `F1_Hz` → "Hz"
#' - `H1_H2c_dB` → "dB"
#' - `Jitter_local_pct` → "pct" (percent)
#' - `intensity` → NA (no unit)
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' .parse_unit_from_colname("fo_Hz")        # "Hz"
#' .parse_unit_from_colname("CPP_dB")      # "dB"
#' .parse_unit_from_colname("frame_time")  # NA
#' }
.parse_unit_from_colname <- function(col_name) {
  # Pattern: _<unit> at end of string
  # Common units: Hz, dB, pct (%), us (microseconds), ms, s
  if (grepl("_([A-Za-z]+)$", col_name)) {
    unit <- sub(".*_([A-Za-z]+)$", "\\1", col_name)
    return(unit)
  }

  NA_character_
}


#' Assign units to data frame columns
#'
#' Assigns R units package units to columns based on their unit suffix.
#'
#' @param df Data frame. Data frame with cleaned column names.
#' @return Data frame with units assigned to appropriate columns.
#'
#' @details
#' Maps common unit suffixes to R units package units:
#' - `_Hz` → Hz
#' - `_dB` → dB (dimensionless, but labeled)
#' - `_pct` → percent (dimensionless)
#' - `_us` → microseconds
#' - `_ms` → milliseconds
#' - `_s` → seconds
#' - `_Bark` → Bark (dimensionless)
#' - `_mel` → mel (dimensionless)
#' - `_ERB` → ERB (dimensionless)
#'
#' Requires the 'units' package to be installed.
#'
#' @keywords internal
#' @examples
#' \dontrun{
#' df <- data.frame(frame_time = 1:10, fo_Hz = rnorm(10, 120, 10))
#' df <- .assign_track_units(df)
#' class(df$fo_Hz)  # "units"
#' }
.assign_track_units <- function(df) {
  if (!requireNamespace("units", quietly = TRUE)) {
    warning("Package 'units' not available. Skipping unit assignment.")
    return(df)
  }

  # Map unit suffixes to R units
  unit_map <- list(
    Hz = "Hz",
    kHz = "kHz",
    dB = "dB",
    pct = "1",      # Percent (dimensionless)
    us = "us",      # Microseconds
    ms = "ms",      # Milliseconds
    s = "s",        # Seconds
    Bark = "1",     # Bark scale (dimensionless)
    mel = "1",      # Mel scale (dimensionless)
    ERB = "1",      # ERB scale (dimensionless)
    semitone = "1"  # Semitones (dimensionless)
  )

  for (col in names(df)) {
    unit_suffix <- .parse_unit_from_colname(col)

    if (!is.na(unit_suffix) && unit_suffix %in% names(unit_map)) {
      unit_str <- unit_map[[unit_suffix]]

      # Try to assign unit
      tryCatch({
        # Special handling for dB (not a standard unit)
        if (unit_suffix == "dB") {
          # Store as attribute since dB is not in udunits
          attr(df[[col]], "unit") <- "dB"
          class(df[[col]]) <- c("dB_vector", "numeric")
        } else {
          df[[col]] <- units::set_units(df[[col]], unit_str, mode = "standard")
        }
      }, error = function(e) {
        # If unit assignment fails, continue without it
        warning(sprintf("Could not assign unit '%s' to column '%s': %s",
                       unit_str, col, e$message), call. = FALSE)
      })
    }
  }

  df
}


#' Get track label mapping
#'
#' Returns a list mapping column names to display labels.
#'
#' @return Named list with 'short' and 'full' label mappings.
#'
#' @details
#' Provides two levels of labels:
#' - **short**: For plot axes (e.g., "fo \[Hz\]", "F1 \[Hz\]")
#' - **full**: For documentation/papers (e.g., "Frequency of oscillation \[Hz\]")
#'
#' @keywords internal
.get_track_label_mapping <- function() {
  list(
    # Short labels (for plots)
    short = list(
      fo_Hz = "fo [Hz]",
      F1_Hz = "F1 [Hz]",
      F2_Hz = "F2 [Hz]",
      F3_Hz = "F3 [Hz]",
      F4_Hz = "F4 [Hz]",
      F5_Hz = "F5 [Hz]",
      F6_Hz = "F6 [Hz]",
      F7_Hz = "F7 [Hz]",
      F8_Hz = "F8 [Hz]",
      B1_Hz = "B1 [Hz]",
      B2_Hz = "B2 [Hz]",
      B3_Hz = "B3 [Hz]",
      B4_Hz = "B4 [Hz]",
      B5_Hz = "B5 [Hz]",
      B6_Hz = "B6 [Hz]",
      B7_Hz = "B7 [Hz]",
      B8_Hz = "B8 [Hz]",
      H1_dB = "H1 [dB]",
      H2_dB = "H2 [dB]",
      H4_dB = "H4 [dB]",
      A1_dB = "A1 [dB]",
      A2_dB = "A2 [dB]",
      A3_dB = "A3 [dB]",
      H2k_dB = "H2k [dB]",
      H5k_dB = "H5k [dB]",
      H1_H2_dB = "H1-H2 [dB]",
      H1_H2c_dB = "H1-H2c [dB]",
      H1_H2u_dB = "H1-H2u [dB]",
      H2_H4_dB = "H2-H4 [dB]",
      H2_H4c_dB = "H2-H4c [dB]",
      H1_A1_dB = "H1-A1 [dB]",
      H1_A1c_dB = "H1-A1c [dB]",
      H1_A2_dB = "H1-A2 [dB]",
      H1_A2c_dB = "H1-A2c [dB]",
      H1_A3_dB = "H1-A3 [dB]",
      H1_A3c_dB = "H1-A3c [dB]",
      H4_H2k_dB = "H4-H2k [dB]",
      H4_H2kc_dB = "H4-H2kc [dB]",
      H2k_H5k_dB = "H2k-H5k [dB]",
      H2k_H5kc_dB = "H2k-H5kc [dB]",
      CPP_dB = "CPP [dB]",
      HNR05_dB = "HNR05 [dB]",
      HNR15_dB = "HNR15 [dB]",
      HNR25_dB = "HNR25 [dB]",
      HNR35_dB = "HNR35 [dB]",
      SHR_dB = "SHR [dB]",
      Energy_dB = "Energy [dB]",
      frame_time = "Time [s]"
    ),

    # Full descriptive labels (for documentation/papers)
    full = list(
      fo_Hz = "Frequency of oscillation [Hz]",
      F1_Hz = "First formant frequency [Hz]",
      F2_Hz = "Second formant frequency [Hz]",
      F3_Hz = "Third formant frequency [Hz]",
      F4_Hz = "Fourth formant frequency [Hz]",
      F5_Hz = "Fifth formant frequency [Hz]",
      F6_Hz = "Sixth formant frequency [Hz]",
      F7_Hz = "Seventh formant frequency [Hz]",
      F8_Hz = "Eighth formant frequency [Hz]",
      B1_Hz = "First formant bandwidth [Hz]",
      B2_Hz = "Second formant bandwidth [Hz]",
      B3_Hz = "Third formant bandwidth [Hz]",
      B4_Hz = "Fourth formant bandwidth [Hz]",
      B5_Hz = "Fifth formant bandwidth [Hz]",
      B6_Hz = "Sixth formant bandwidth [Hz]",
      B7_Hz = "Seventh formant bandwidth [Hz]",
      B8_Hz = "Eighth formant bandwidth [Hz]",
      H1_dB = "First harmonic amplitude [dB]",
      H2_dB = "Second harmonic amplitude [dB]",
      H4_dB = "Fourth harmonic amplitude [dB]",
      A1_dB = "Amplitude of harmonic nearest F1 [dB]",
      A2_dB = "Amplitude of harmonic nearest F2 [dB]",
      A3_dB = "Amplitude of harmonic nearest F3 [dB]",
      H2k_dB = "Harmonic nearest 2 kHz [dB]",
      H5k_dB = "Harmonic nearest 5 kHz [dB]",
      H1_H2_dB = "H1-H2 amplitude difference [dB]",
      H1_H2c_dB = "H1-H2 corrected for formants [dB]",
      H1_H2u_dB = "H1-H2 uncorrected [dB]",
      H2_H4_dB = "H2-H4 amplitude difference [dB]",
      H2_H4c_dB = "H2-H4 corrected for formants [dB]",
      H1_A1_dB = "H1-A1 spectral tilt [dB]",
      H1_A1c_dB = "H1-A1 corrected for formants [dB]",
      H1_A2_dB = "H1-A2 spectral tilt [dB]",
      H1_A2c_dB = "H1-A2 corrected for formants [dB]",
      H1_A3_dB = "H1-A3 spectral tilt [dB]",
      H1_A3c_dB = "H1-A3 corrected for formants [dB]",
      H4_H2k_dB = "H4-H2k mid-frequency tilt [dB]",
      H4_H2kc_dB = "H4-H2k corrected for formants [dB]",
      H2k_H5k_dB = "H2k-H5k high-frequency tilt [dB]",
      H2k_H5kc_dB = "H2k-H5k corrected for formants [dB]",
      CPP_dB = "Cepstral Peak Prominence [dB]",
      HNR05_dB = "Harmonics-to-Noise Ratio 0-500 Hz [dB]",
      HNR15_dB = "Harmonics-to-Noise Ratio 0-1500 Hz [dB]",
      HNR25_dB = "Harmonics-to-Noise Ratio 0-2500 Hz [dB]",
      HNR35_dB = "Harmonics-to-Noise Ratio 0-3500 Hz [dB]",
      SHR_dB = "Subharmonic-to-Harmonic Ratio [dB]",
      Energy_dB = "Energy [dB]",
      frame_time = "Frame time [s]"
    )
  )
}


#' Infer track label from column name
#'
#' Intelligently infers a display label from a cleaned column name.
#'
#' @param col Character. Column name.
#' @param full Logical. If TRUE, generate full descriptive label.
#' @return Character. Display label.
#'
#' @details
#' First checks the predefined mapping from `.get_track_label_mapping()`.
#' If not found, applies intelligent inference rules:
#' - Formants: `F1_Hz` → "F1 \[Hz\]" or "First formant frequency \[Hz\]"
#' - Bandwidths: `B1_Hz` → "B1 \[Hz\]" or "First formant bandwidth \[Hz\]"
#' - Generic: `param_unit` → "param [unit]"
#'
#' @keywords internal
.infer_track_label <- function(col, full = FALSE) {
  # Try to get from mapping
  mapping <- .get_track_label_mapping()
  labels <- if (full) mapping$full else mapping$short

  if (!is.null(labels[[col]])) {
    return(labels[[col]])
  }

  # Intelligent inference for formants
  if (grepl("^F([0-9]+)_Hz$", col)) {
    n <- as.integer(gsub("^F([0-9]+)_Hz$", "\\1", col))
    if (full && n <= 8) {
      ord <- c("First", "Second", "Third", "Fourth", "Fifth", "Sixth", "Seventh", "Eighth")
      return(paste0(ord[n], " formant frequency [Hz]"))
    } else {
      return(paste0("F", n, " [Hz]"))
    }
  }

  # Intelligent inference for bandwidths
  if (grepl("^B([0-9]+)_Hz$", col)) {
    n <- as.integer(gsub("^B([0-9]+)_Hz$", "\\1", col))
    if (full && n <= 8) {
      ord <- c("First", "Second", "Third", "Fourth", "Fifth", "Sixth", "Seventh", "Eighth")
      return(paste0(ord[n], " formant bandwidth [Hz]"))
    } else {
      return(paste0("B", n, " [Hz]"))
    }
  }

  # Generic: param_unit → "param [unit]"
  if (grepl("_([A-Za-z]+)$", col)) {
    return(gsub("_([A-Za-z]+)$", " [\\1]", col))
  }

  # No conversion possible
  col
}


#' Generate track labels for data frame
#'
#' Creates a named list mapping column names to short display labels.
#'
#' @param col_names Character vector. Column names.
#' @return Named list mapping column names to labels.
#'
#' @keywords internal
.generate_track_labels <- function(col_names) {
  labels <- list()

  for (col in col_names) {
    labels[[col]] <- .infer_track_label(col, full = FALSE)
  }

  labels
}


#' Generate track descriptions for data frame
#'
#' Creates a named list mapping column names to full descriptive labels.
#'
#' @param col_names Character vector. Column names.
#' @return Named list mapping column names to descriptions.
#'
#' @keywords internal
.generate_track_descriptions <- function(col_names) {
  descriptions <- list()

  for (col in col_names) {
    descriptions[[col]] <- .infer_track_label(col, full = TRUE)
  }

  descriptions
}

#!/usr/bin/env Rscript
#
# Validate Track Names Against Standards
#
# This script validates all track names in the codebase against:
# - Titze 2015 standards
# - Nylén 2024 extensions
# - superassp conventions
#
# Usage: Rscript scripts/validate_track_names.R

library(stringr)

cat("=== Track Name Validation ===\n\n")

# Load mapping
CSV_FILE <- "TRACK_NAMES_MAPPING.csv"
if (!file.exists(CSV_FILE)) {
  stop("Mapping CSV not found: ", CSV_FILE)
}

mapping <- read.csv(CSV_FILE, stringsAsFactors = FALSE)
cat("Loaded", nrow(mapping), "track definitions\n\n")

# Validation rules based on Titze 2015 + Nylén 2024

validate_f0 <- function(name) {
  errors <- c()

  # Should be lowercase 'fo' (letter o for oscillation) not 'F0' or 'f0'
  if (grepl("^F0", name) || grepl("^f0", name)) {
    errors <- c(errors, "Should use 'fo' with letter o for oscillation (Titze 2015)")
  }

  # Should have unit [Hz]
  if (!grepl("\\[Hz\\]", name)) {
    errors <- c(errors, "Missing unit [Hz]")
  }

  errors
}

validate_formant <- function(name) {
  errors <- c()

  # Should be uppercase 'F' with number
  if (grepl("^f[0-9]", name)) {
    errors <- c(errors, "Should use uppercase 'F' for formants (Titze 2015)")
  }

  # Should have subscript number (F1, F2, not fm)
  if (name == "fm" || name == "Fi") {
    errors <- c(errors, "Should use specific formant number (F1, F2, F3...) when possible")
  }

  # Should have unit [Hz]
  if (!grepl("\\[Hz\\]", name) && name != "fm") {
    errors <- c(errors, "Missing unit [Hz]")
  }

  errors
}

validate_bandwidth <- function(name) {
  errors <- c()

  # Should be 'B' with number
  if (grepl("^b[0-9]", name)) {
    errors <- c(errors, "Should use uppercase 'B' for bandwidth (Titze 2015)")
  }

  # Should have specific number (not generic 'bw')
  if (name == "bw" || name == "Bi") {
    errors <- c(errors, "Should use specific bandwidth number (B1, B2, B3...) when possible")
  }

  # Should have unit [Hz]
  if (!grepl("\\[Hz\\]", name) && name != "bw") {
    errors <- c(errors, "Missing unit [Hz]")
  }

  errors
}

validate_harmonic <- function(name) {
  errors <- c()

  # Should have unit [dB]
  if (!grepl("\\[dB\\]", name)) {
    errors <- c(errors, "Missing unit [dB]")
  }

  # Check format: H1, H2, A1, A2, H2k (not H2K)
  if (grepl("K", name)) {
    errors <- c(errors, "Should use lowercase 'k' for kHz (H2k not H2K)")
  }

  errors
}

validate_harmonic_diff <- function(name) {
  errors <- c()

  # Should have hyphen separator
  if (grepl("H[0-9][HA]", name) && !grepl("-", name)) {
    errors <- c(errors, "Should use hyphen separator (H1-H2 not H1H2)")
  }

  # Should have unit [dB]
  if (!grepl("\\[dB\\]", name)) {
    errors <- c(errors, "Missing unit [dB]")
  }

  # Check corrected/uncorrected notation
  if (grepl("[cu]$", gsub("\\[dB\\]", "", name))) {
    # Has suffix, should be before unit
    if (!grepl("[cu]\\[dB\\]$", name)) {
      errors <- c(errors, "Correction suffix should be before unit (H1-H2c[dB])")
    }
  }

  errors
}

validate_voice_quality <- function(name) {
  errors <- c()

  # CPP, HNR, SHR should have [dB]
  if (grepl("^(CPP|HNR|SHR)", name) && !grepl("\\[dB\\]", name)) {
    errors <- c(errors, "Missing unit [dB]")
  }

  # Jitter should have [%] or [us]
  if (grepl("Jitter", name) && !grepl("\\[(%|us)\\]", name)) {
    errors <- c(errors, "Missing unit [%] or [us]")
  }

  # Shimmer should have [%] or [dB]
  if (grepl("Shimmer", name) && !grepl("\\[(%|dB)\\]", name)) {
    errors <- c(errors, "Missing unit [%] or [dB]")
  }

  # Should not use parentheses (use underscores)
  if (grepl("\\(", name)) {
    errors <- c(errors, "Should use underscores instead of parentheses (Jitter_local not 'Jitter (local)')")
  }

  # Energy should have [dB]
  if (grepl("^Energy", name) && !grepl("\\[dB\\]", name)) {
    errors <- c(errors, "Missing unit [dB]")
  }

  errors
}

validate_spectral <- function(name) {
  errors <- c()

  # RMS should have [dB]
  if (grepl("^RMS", name) && !grepl("\\[dB\\]", name)) {
    errors <- c(errors, "Missing unit [dB]")
  }

  # ZCR should have [Hz]
  if (grepl("^ZCR", name) && !grepl("\\[Hz\\]", name)) {
    errors <- c(errors, "Missing unit [Hz]")
  }

  # MFCC could be uppercase (optional)
  if (grepl("^mfcc", name)) {
    # This is stylistic, not an error
  }

  errors
}

# Run validation
results <- list()

for (i in 1:nrow(mapping)) {
  row <- mapping[i, ]
  errors <- c()

  # Validate new name based on category
  if (row$category == "f0") {
    errors <- validate_f0(row$new_name)
  } else if (row$category == "formant") {
    errors <- validate_formant(row$new_name)
  } else if (row$category == "bandwidth") {
    errors <- validate_bandwidth(row$new_name)
  } else if (row$category == "harmonic") {
    errors <- validate_harmonic(row$new_name)
  } else if (row$category == "harmonic_diff") {
    errors <- validate_harmonic_diff(row$new_name)
  } else if (row$category == "voice_quality") {
    errors <- validate_voice_quality(row$new_name)
  } else if (row$category == "spectral") {
    errors <- validate_spectral(row$new_name)
  }

  if (length(errors) > 0) {
    results[[length(results) + 1]] <- list(
      function_name = row$function_name,
      old_name = row$old_name,
      new_name = row$new_name,
      category = row$category,
      errors = errors
    )
  }
}

# Report results
cat("=== Validation Results ===\n\n")

if (length(results) == 0) {
  cat("✓ All new track names comply with standards!\n\n")
  cat("Standards checked:\n")
  cat("  - Titze 2015: Symbolic notation for fo (oscillation), formants, harmonics\n")
  cat("  - Nylén 2024: Corrected formant notation\n")
  cat("  - Explicit units for all measurements\n\n")
} else {
  cat(sprintf("✗ Found %d validation issues\n\n", length(results)))

  # Group by category
  by_category <- split(results, sapply(results, function(x) x$category))

  for (cat_name in names(by_category)) {
    cat(sprintf("## Category: %s (%d issues)\n\n", cat_name, length(by_category[[cat_name]])))

    for (item in by_category[[cat_name]]) {
      cat(sprintf("**%s**: `%s` → `%s`\n", item$function_name, item$old_name, item$new_name))
      for (error in item$errors) {
        cat(sprintf("  - %s\n", error))
      }
      cat("\n")
    }
  }
}

# Summary statistics
cat("=== Summary Statistics ===\n\n")
cat(sprintf("Total tracks validated: %d\n", nrow(mapping)))
cat(sprintf("Tracks with issues: %d\n", length(results)))
cat(sprintf("Tracks passing validation: %d (%.1f%%)\n",
            nrow(mapping) - length(results),
            100 * (nrow(mapping) - length(results)) / nrow(mapping)))

if (length(results) > 0) {
  cat("\nIssues by category:\n")
  issue_counts <- table(sapply(results, function(x) x$category))
  print(issue_counts)
}

cat("\n")

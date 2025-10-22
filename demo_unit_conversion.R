#!/usr/bin/env Rscript
#
# Demonstration of Automatic Unit Conversion in AsspDataObj
#
# This script demonstrates the new automatic unit conversion feature
# for as.data.frame() and as_tibble() methods on AsspDataObj objects.

library(superassp)
library(tibble)
library(units)

cat("=== Automatic Unit Conversion in AsspDataObj ===\n\n")

# Load test audio file
test_wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
cat("Test file:", test_wav, "\n\n")

# ---- Example 1: F0 Tracking with Unit Conversion ----
cat("## Example 1: F0 Tracking (Pitch)\n\n")

f0_obj <- fo(test_wav, toFile = FALSE, verbose = FALSE)
cat("AsspDataObj track names:", paste(names(f0_obj), collapse = ", "), "\n\n")

# Convert to data frame WITH unit conversion (default)
cat("### With unit conversion (convert_units = TRUE):\n")
df_with_units <- as.data.frame(f0_obj, convert_units = TRUE)
cat("Column 'fo[Hz]' class:", class(df_with_units[["fo[Hz]"]]), "\n")
cat("Unit:", as.character(units::deparse_unit(df_with_units[["fo[Hz]"]])), "\n")
cat("First 10 values:\n")
print(head(df_with_units[["fo[Hz]"]], 10))
cat("\n")

# Convert to data frame WITHOUT unit conversion
cat("### Without unit conversion (convert_units = FALSE):\n")
df_no_units <- as.data.frame(f0_obj, convert_units = FALSE)
cat("Column 'fo[Hz]' class:", class(df_no_units[["fo[Hz]"]]), "\n")
cat("First 10 values:\n")
print(head(df_no_units[["fo[Hz]"]], 10))
cat("\n")

# ---- Example 2: Unit Conversion with Tibbles ----
cat("## Example 2: Converting to Tibble\n\n")

tbl <- as_tibble(f0_obj, convert_units = TRUE)
cat("Tibble dimensions:", nrow(tbl), "x", ncol(tbl), "\n")
cat("Columns:", paste(names(tbl), collapse = ", "), "\n")
cat("Column 'fo[Hz]' is units?", inherits(tbl[["fo[Hz]"]], "units"), "\n")
cat("\nFirst 5 rows:\n")
print(head(tbl, 5))
cat("\n")

# ---- Example 3: Integration with Psychoacoustic Scales ----
cat("## Example 3: Integration with Psychoacoustic Scales\n\n")

# Extract Hz values and convert to Bark scale
hz_values <- df_with_units[["fo[Hz]"]]
cat("Original F0 range (Hz):\n")
cat("  Min:", min(hz_values[!is.na(hz_values)]), "\n")
cat("  Max:", max(hz_values[!is.na(hz_values)]), "\n")
cat("  Mean:", mean(hz_values[!is.na(hz_values)]), "\n\n")

# Convert non-zero Hz to Bark
nonzero_hz <- as.numeric(hz_values[hz_values > 0])
if (length(nonzero_hz) > 0) {
  bark_values <- hz_to_bark(nonzero_hz, as_units = TRUE)
  cat("Converted to Bark scale:\n")
  cat("  Min:", min(bark_values, na.rm = TRUE), "\n")
  cat("  Max:", max(bark_values, na.rm = TRUE), "\n")
  cat("  Mean:", mean(bark_values, na.rm = TRUE), "\n")
  cat("  Unit:", as.character(units::deparse_unit(bark_values)), "\n\n")
}

# ---- Example 4: Column Name Pattern Matching ----
cat("## Example 4: Unit Label Pattern Matching\n\n")

test_cases <- c(
  "fo[Hz]",         # Match: Hz
  "intensity[dB]",  # Match: dB
  "fm[kHz]",        # Match: kHz
  "pitch[Bark]",    # Match: Bark
  "frame_time",     # No match: no unit label
  "fo[Hz]_1",       # No match: unit not at end
  "data"            # No match: no brackets
)

cat("Testing unit label extraction:\n")
for (test_case in test_cases) {
  unit <- .parse_unit_from_colname(test_case)
  cat(sprintf("  %-20s -> %s\n", paste0("\"", test_case, "\""),
              ifelse(is.na(unit), "(no unit)", unit)))
}
cat("\n")

# ---- Example 5: Unit Conversion with NA Handling ----
cat("## Example 5: NA Handling (na.zeros = TRUE)\n\n")

tbl_with_na <- as_tibble(f0_obj, convert_units = TRUE, na.zeros = TRUE)
tbl_without_na <- as_tibble(f0_obj, convert_units = TRUE, na.zeros = FALSE)

count_na_with <- sum(is.na(tbl_with_na[["fo[Hz]"]]))
count_na_without <- sum(is.na(tbl_without_na[["fo[Hz]"]]))
count_zero_without <- sum(as.numeric(tbl_without_na[["fo[Hz]"]]) == 0, na.rm = TRUE)

cat("With na.zeros = TRUE:\n")
cat("  NA count:", count_na_with, "\n\n")

cat("With na.zeros = FALSE:\n")
cat("  NA count:", count_na_without, "\n")
cat("  Zero count:", count_zero_without, "\n\n")

cat("Note: Zero values in pitch tracks typically indicate\n")
cat("      unvoiced regions. Setting na.zeros=TRUE converts\n")
cat("      them to NA for easier statistical analysis.\n\n")

# ---- Summary ----
cat("## Summary\n\n")
cat("Benefits of automatic unit conversion:\n")
cat("  - Explicit units make code more readable\n")
cat("  - Units package provides automatic unit checking\n")
cat("  - Easy integration with psychoacoustic scales\n")
cat("  - Prevents unit confusion in data analysis\n")
cat("  - Opt-out available with convert_units = FALSE\n\n")

cat("Track naming convention:\n")
cat("  - Unit labels must end with [<unit>]\n")
cat("  - Examples: fo[Hz], intensity[dB], pitch[Bark]\n")
cat("  - Without unit labels, columns remain numeric\n\n")

cat("=== Demonstration Complete ===\n")

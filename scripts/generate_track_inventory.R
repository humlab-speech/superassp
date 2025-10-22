#!/usr/bin/env Rscript
#
# Generate Comprehensive Track Name Inventory
#
# This script analyzes all DSP functions in superassp to create a complete
# inventory of track names, categorize them, and generate statistics.
#
# Usage: Rscript scripts/generate_track_inventory.R

library(stringr)

# Configuration
BASE_DIR <- "."
R_DIR <- file.path(BASE_DIR, "R")
OUTPUT_FILE <- file.path(BASE_DIR, "TRACK_INVENTORY.md")
CSV_FILE <- file.path(BASE_DIR, "TRACK_NAMES_MAPPING.csv")

cat("=== Track Name Inventory Generator ===\n\n")
cat("Scanning directory:", R_DIR, "\n")

# Read mapping CSV
if (file.exists(CSV_FILE)) {
  mapping <- read.csv(CSV_FILE, stringsAsFactors = FALSE)
  cat("Loaded mapping with", nrow(mapping), "entries\n\n")
} else {
  stop("Mapping CSV not found: ", CSV_FILE)
}

# Generate statistics
stats <- list()

# 1. Count by category
stats$by_category <- table(mapping$category)

# 2. Count by function
stats$by_function <- table(mapping$function_name)

# 3. Count requiring changes
stats$needs_change <- sum(mapping$old_name != mapping$new_name)
stats$no_change <- sum(mapping$old_name == mapping$new_name)

# 4. Count by change type
stats$needs_unit <- sum(!grepl("\\[", mapping$old_name) & grepl("\\[", mapping$new_name))
stats$needs_case_change <- sum(grepl("^F0", mapping$old_name) & grepl("^f0", mapping$new_name))
stats$needs_hyphen <- sum(grepl("H[0-9][HA]", mapping$old_name) & grepl("-", mapping$new_name))

# 5. Identify unique old and new names
stats$unique_old <- length(unique(mapping$old_name))
stats$unique_new <- length(unique(mapping$new_name))

# Generate markdown report
md_lines <- c(
  "# Track Name Inventory",
  "",
  paste("Generated:", Sys.time()),
  "",
  "## Summary Statistics",
  "",
  paste("- **Total track definitions**:", nrow(mapping)),
  paste("- **Unique old names**:", stats$unique_old),
  paste("- **Unique new names**:", stats$unique_new),
  paste("- **Functions affected**:", length(stats$by_function)),
  paste("- **Tracks needing changes**:", stats$needs_change, sprintf("(%.1f%%)", 100 * stats$needs_change / nrow(mapping))),
  paste("- **Tracks unchanged**:", stats$no_change, sprintf("(%.1f%%)", 100 * stats$no_change / nrow(mapping))),
  "",
  "## Changes Required",
  "",
  paste("- **Add units**:", stats$needs_unit, "tracks"),
  paste("- **Case normalization**:", stats$needs_case_change, "tracks"),
  paste("- **Add hyphens**:", stats$needs_hyphen, "tracks (harmonic differences)"),
  ""
)

# Category breakdown
md_lines <- c(md_lines,
  "## Tracks by Category",
  "",
  "| Category | Count | Percentage |",
  "|----------|-------|------------|"
)

for (cat_name in names(sort(stats$by_category, decreasing = TRUE))) {
  count <- stats$by_category[cat_name]
  pct <- sprintf("%.1f%%", 100 * count / nrow(mapping))
  md_lines <- c(md_lines, sprintf("| %s | %d | %s |", cat_name, count, pct))
}

md_lines <- c(md_lines, "")

# Function breakdown (top 20)
md_lines <- c(md_lines,
  "## Functions with Most Track Definitions (Top 20)",
  "",
  "| Function | Track Count | File |",
  "|----------|-------------|------|"
)

func_counts <- sort(stats$by_function, decreasing = TRUE)[1:min(20, length(stats$by_function))]
for (func_name in names(func_counts)) {
  count <- func_counts[func_name]
  file <- unique(mapping$file[mapping$function_name == func_name])[1]
  md_lines <- c(md_lines, sprintf("| `%s` | %d | %s |", func_name, count, file))
}

md_lines <- c(md_lines, "")

# Detailed category tables
for (cat_name in c("f0", "formant", "bandwidth", "harmonic", "harmonic_diff", "voice_quality")) {
  cat_data <- mapping[mapping$category == cat_name, ]

  if (nrow(cat_data) == 0) next

  md_lines <- c(md_lines,
    sprintf("## Category: %s (%d tracks)", cat_name, nrow(cat_data)),
    ""
  )

  # Show unique mappings
  unique_mappings <- unique(cat_data[, c("old_name", "new_name", "unit", "notes")])
  unique_mappings <- unique_mappings[order(unique_mappings$old_name), ]

  md_lines <- c(md_lines,
    "| Old Name | New Name | Unit | Notes |",
    "|----------|----------|------|-------|"
  )

  for (i in 1:nrow(unique_mappings)) {
    row <- unique_mappings[i, ]
    md_lines <- c(md_lines, sprintf("| `%s` | `%s` | %s | %s |",
                                   row$old_name, row$new_name,
                                   ifelse(row$unit == "", "—", row$unit),
                                   row$notes))
  }

  md_lines <- c(md_lines, "")

  # List functions using this category
  funcs <- unique(cat_data$function_name)
  md_lines <- c(md_lines,
    sprintf("**Functions** (%d): %s", length(funcs), paste(paste0("`", funcs, "`"), collapse = ", ")),
    ""
  )
}

# Appendix: Complete mapping table
md_lines <- c(md_lines,
  "## Appendix: Complete Mapping Table",
  "",
  "| Function | File | Old Name | New Name | Category | Unit | Notes |",
  "|----------|------|----------|----------|----------|------|-------|"
)

for (i in 1:nrow(mapping)) {
  row <- mapping[i, ]
  md_lines <- c(md_lines, sprintf("| `%s` | %s | `%s` | `%s` | %s | %s | %s |",
                                 row$function_name, row$file, row$old_name, row$new_name,
                                 row$category,
                                 ifelse(row$unit == "", "—", row$unit),
                                 row$notes))
}

# Write output
writeLines(md_lines, OUTPUT_FILE)

cat("\n=== Inventory Report ===\n\n")
cat("Total tracks:", nrow(mapping), "\n")
cat("Categories:", length(stats$by_category), "\n")
cat("Functions:", length(stats$by_function), "\n")
cat("Needs changes:", stats$needs_change, sprintf("(%.1f%%)\n", 100 * stats$needs_change / nrow(mapping)))
cat("\nReport written to:", OUTPUT_FILE, "\n")
cat("\nTop 5 categories:\n")
print(head(sort(stats$by_category, decreasing = TRUE), 5))
cat("\nTop 5 functions:\n")
print(head(sort(stats$by_function, decreasing = TRUE), 5))

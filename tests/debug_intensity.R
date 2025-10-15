#!/usr/bin/env Rscript
# Debug script for praat_intensity_opt

library(superassp)

# Setup Python
Sys.unsetenv("RETICULATE_PYTHON")
reticulate::use_python("/usr/bin/python3", required = FALSE)

test_file <- "signalfiles/AVQI/input/sv1.wav"

# Source the Python script
reticulate::source_python(system.file("python", "praat_intensity.py", package = "superassp"))

# Call Python function directly
cat("Testing Python function directly...\n")
py_windowShape <- reticulate::py_eval("pm.WindowShape.GAUSSIAN1")

result_df <- reticulate::py$praat_intensity(
  normalizePath(test_file, mustWork = TRUE),
  beginTime = 0.0,
  endTime = 0.0,
  time_step = 0.0,  # Automatic
  minimal_f0_frequency = 50.0,
  subtract_mean = TRUE,
  windowShape = py_windowShape,
  relativeWidth = 1.0
)

cat("\nDataFrame structure:\n")
print(str(result_df))

cat("\nDataFrame dimensions:\n")
print(dim(result_df))

cat("\nDataFrame head:\n")
print(head(result_df))

cat("\nColumn names:\n")
print(names(result_df))

if(nrow(result_df) > 0) {
  cat("\nColumn classes:\n")
  print(sapply(result_df, class))

  cat("\nTrying to extract intensity column:\n")
  int_col <- result_df$`Intensity(dB)`
  cat("Intensity class:", class(int_col), "\n")
  cat("Intensity length:", length(int_col), "\n")
  cat("Intensity head:", head(int_col), "\n")
} else {
  cat("\n*** DataFrame is EMPTY! ***\n")
}

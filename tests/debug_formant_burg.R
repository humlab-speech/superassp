#!/usr/bin/env Rscript
# Debug script for praat_formant_burg_opt

library(superassp)

# Setup Python
Sys.unsetenv("RETICULATE_PYTHON")
reticulate::use_python("/usr/bin/python3", required = FALSE)

test_file <- "signalfiles/AVQI/input/sv1.wav"

# Source the Python script
reticulate::source_python(system.file("python", "praat_formant_burg.py", package = "superassp"))

# Call Python function directly
cat("Testing Python function directly...\n")
py_windowShape <- reticulate::py_eval("pm.WindowShape.GAUSSIAN1")
py_spec_shape <- reticulate::py_eval("pm.SpectralAnalysisWindowShape.GAUSSIAN")

result_df <- reticulate::py$praat_formant_burg(
  normalizePath(test_file, mustWork = TRUE),
  beginTime = 0.0,
  endTime = 0.0,
  timeStep = 0.005,
  number_of_formants = 5,
  maxHzFormant = 5500.0,
  windowLength = 0.025,
  pre_emphasis = 50.0,
  track_formants = FALSE,
  number_of_tracks = 3,
  reference_F1 = 550,
  reference_F2 = 1650,
  reference_F3 = 2750,
  reference_F4 = 3850,
  reference_F5 = 4950,
  frequency_cost = 1.0,
  bandwidth_cost = 1.0,
  transition_cost = 1.0,
  windowShape = py_windowShape,
  relativeWidth = 1.0,
  spectrogram_window_shape = py_spec_shape,
  spectrogram_resolution = 40.0
)

cat("\nDataFrame structure:\n")
print(str(result_df))

cat("\nDataFrame head:\n")
print(head(result_df))

cat("\nDataFrame classes:\n")
print(sapply(result_df, class))

cat("\nColumn names:\n")
print(names(result_df))

cat("\nTrying to extract F1 column:\n")
f1_col <- result_df$`F1(Hz)`
cat("F1 class:", class(f1_col), "\n")
cat("F1 type:", typeof(f1_col), "\n")
cat("F1 head:", head(f1_col), "\n")

cat("\nTrying to convert F1 to matrix:\n")
tryCatch({
  f1_matrix <- as.matrix(f1_col)
  cat("Success! Matrix dimensions:", dim(f1_matrix), "\n")
  cat("Matrix class:", class(f1_matrix), "\n")
}, error = function(e) {
  cat("Error:", conditionMessage(e), "\n")

  # Try converting to numeric first
  cat("\nTrying to convert to numeric first:\n")
  f1_numeric <- as.numeric(f1_col)
  cat("Numeric class:", class(f1_numeric), "\n")
  cat("Numeric head:", head(f1_numeric), "\n")

  f1_matrix <- as.matrix(f1_numeric)
  cat("Matrix after numeric conversion - dimensions:", dim(f1_matrix), "\n")
})

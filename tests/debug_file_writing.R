#!/usr/bin/env Rscript
# Debug file writing issue

library(superassp)

# Setup Python
Sys.unsetenv("RETICULATE_PYTHON")
reticulate::use_python("/usr/bin/python3", required = FALSE)

test_file <- "signalfiles/AVQI/input/sv1.wav"

# Test intensity which works for toFile=FALSE
cat("Testing intensity with toFile=FALSE...\n")
result <- praat_intensity_opt(
  test_file,
  toFile = FALSE
)

cat("\nAttributes of AsspDataObj:\n")
print(attributes(result))

cat("\nTrack formats:\n")
print(attr(result, "trackFormats"))

cat("\n\nNow trying to write to file...\n")
temp_dir <- tempdir()
temp_file <- file.path(temp_dir, "test.int")

tryCatch({
  wrassp::write.AsspDataObj(dobj = result, file = temp_file)
  cat("✓ Successfully written!\n")
  cat("File size:", file.size(temp_file), "bytes\n")
  unlink(temp_file)
}, error = function(e) {
  cat("✗ Error:", conditionMessage(e), "\n")
})

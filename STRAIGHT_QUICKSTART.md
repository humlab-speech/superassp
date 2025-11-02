# STRAIGHT Legacy Vocoder - Quick Start Guide

**Package**: superassp v0.8.x  
**Date**: October 29, 2025  
**Status**: ✅ Ready for testing

---

## Installation

### Step 1: Ensure superassp is loaded

```r
# Load the package (or use devtools::load_all() for development)
library(superassp)
```

### Step 2: Install Python dependencies

```r
# Install with Numba optimization (recommended - 20% faster)
install_legacy_straight(install_numba = TRUE)

# Or install without optimization
install_legacy_straight(install_numba = FALSE)
```

### Step 3: Verify installation

```r
# Check if available
straight_available()  # Should return TRUE

# Get detailed information
straight_info()
```

---

## Basic Usage

### F0 Extraction

```r
# Single file
wav_file <- "speech.wav"
f0_data <- trk_straight_f0(wav_file, toFile = FALSE)

# Plot F0 contour
plot(f0_data$f0[,1], type = "l", ylab = "F0 (Hz)", main = "STRAIGHT F0")

# Check voice/unvoiced decisions
table(f0_data$vuv[,1])
```

###Spectral Analysis

```r
# Extract spectral envelope
spec_data <- trk_straight_spec(wav_file, toFile = FALSE)

# Visualize spectrogram
image(t(spec_data$spec), col = heat.colors(256),
      xlab = "Frame", ylab = "Frequency", main = "STRAIGHT Spectrum")
```

### Speech Synthesis

```r
# Extract parameters
f0_data <- trk_straight_f0(wav_file, toFile = FALSE)
spec_data <- trk_straight_spec(wav_file, toFile = FALSE)

# Synthesize
audio_synth <- straight_synth(
  f0 = f0_data$f0[,1],
  spec = spec_data$spec,
  sample_rate = attr(f0_data, "sampleRate")
)

# Save
av::av_audio_convert(audio_synth, "resynthesis.wav",
                     sample_rate = attr(f0_data, "sampleRate"))
```

---

## Complete Pipeline

```r
# One-step analysis-synthesis
result <- straight_pipeline(
  input_file = "speech.wav",
  output_file = "resynthesis.wav",
  f0_floor = 60,
  f0_ceil = 400,
  verbose = TRUE
)

# Access components
plot(result$f0$f0[,1], type = "l")
image(t(result$spec$spec))
```

---

## Voice Modification Examples

### Pitch Shifting

```r
# Load audio
wav_file <- "speech.wav"

# Extract parameters
f0_data <- trk_straight_f0(wav_file, toFile = FALSE)
spec_data <- trk_straight_spec(wav_file, toFile = FALSE)

# Shift pitch up 50%
f0_modified <- f0_data$f0[,1] * 1.5

# Synthesize
audio_shifted <- straight_synth(
  f0 = f0_modified,
  spec = spec_data$spec,
  sample_rate = attr(f0_data, "sampleRate"),
  output_file = "pitch_shifted.wav"
)
```

### Voice Conversion (Basic)

```r
# Extract parameters from source and target
source_f0 <- trk_straight_f0("source.wav", toFile = FALSE)
target_f0 <- trk_straight_f0("target.wav", toFile = FALSE)
source_spec <- trk_straight_spec("source.wav", toFile = FALSE)

# Apply target F0 to source spectrum (simplified)
audio_converted <- straight_synth(
  f0 = target_f0$f0[,1],
  spec = source_spec$spec,
  sample_rate = attr(source_f0, "sampleRate"),
  output_file = "voice_converted.wav"
)
```

---

## Batch Processing

```r
# Process multiple files
files <- list.files("audio/", pattern = "\\.wav$", full.names = TRUE)

# Extract F0 for all files
f0_results <- trk_straight_f0(
  files,
  f0_floor = 80,
  f0_ceil = 400,
  toFile = FALSE,
  verbose = TRUE
)

# Extract spectral envelopes
spec_results <- trk_straight_spec(
  files,
  fft_size = 2048,
  toFile = FALSE,
  verbose = TRUE
)

# Access individual results
f0_file1 <- f0_results[[1]]
spec_file1 <- spec_results[[1]]
```

---

## Writing to SSFF Files (emuR Integration)

```r
# Write F0 to SSFF file
trk_straight_f0("speech.wav", toFile = TRUE, outputDirectory = "features/")
# Creates: features/speech.strf0

# Write spectral envelope to SSFF file
trk_straight_spec("speech.wav", toFile = TRUE, outputDirectory = "features/")
# Creates: features/speech.strspec

# Use in emuR database
library(reindeer)
corpus <- corpus("mydb_emuDB")
f0_data <- ask_for(corpus, "STRAIGHT_F0")
```

---

## Performance Tips

### Use Numba Optimization

```r
# Install with Numba for ~20% speedup
install_legacy_straight(install_numba = TRUE)

# Verify Numba is available
info <- straight_available(detailed = TRUE)
info$numba_available  # Should be TRUE
info$optimization     # Should say "Numba JIT (~20% faster)"
```

### First Run Compilation

```r
# First run will compile Numba functions (~0.5s overhead)
# Subsequent runs are faster

# Warm-up call (optional)
dummy_audio <- rnorm(16000)  # 1 second at 16kHz
# ... process dummy audio to compile functions
```

### Process Long Files in Chunks

```r
# For very long recordings
chunk_duration <- 30  # seconds
total_duration <- av::av_audio_info("long.wav")$duration

results <- list()
for (i in seq(0, total_duration, chunk_duration)) {
  results[[length(results) + 1]] <- trk_straight_f0(
    "long.wav",
    beginTime = i,
    endTime = min(i + chunk_duration, total_duration),
    toFile = FALSE
  )
}
```

---

## Troubleshooting

### Module Not Available

```r
# Check status
straight_info()

# Reinstall
install_legacy_straight(install_numba = TRUE)

# Check Python configuration
reticulate::py_config()
```

### Slow Performance

```r
# Install Numba if not already installed
install_legacy_straight(install_numba = TRUE)

# Check optimization status
info <- straight_available(detailed = TRUE)
if (!info$numba_available) {
  message("Numba not available - install for 20% speedup")
}
```

### Import Errors

```r
# Check Python path
library(reticulate)
py_run_string("import sys; print(sys.path)")

# Verify module location
inst_path <- system.file("python/legacy_STRAIGHT", package = "superassp")
print(inst_path)
dir.exists(inst_path)  # Should be TRUE
```

---

## Comparison with Other F0 Trackers

```r
# STRAIGHT (high quality, vocoding)
f0_straight <- trk_straight_f0(wav_file, toFile = FALSE)

# RAPT (fast, robust)
f0_rapt <- trk_rapt(wav_file, toFile = FALSE)

# SWIPE (noise robust)
f0_swipe <- trk_swipe(wav_file, toFile = FALSE)

# Compare results
plot(f0_straight$f0[,1], type = "l", col = "red", ylab = "F0 (Hz)")
lines(f0_rapt$f0[,1], col = "blue")
lines(f0_swipe$f0[,1], col = "green")
legend("topright", c("STRAIGHT", "RAPT", "SWIPE"),
       col = c("red", "blue", "green"), lty = 1)
```

---

## Advanced Usage

### Custom F0 Range

```r
# For male speech
f0_male <- trk_straight_f0(wav_file,
                           f0_floor = 50,
                           f0_ceil = 250,
                           toFile = FALSE)

# For female speech
f0_female <- trk_straight_f0(wav_file,
                             f0_floor = 120,
                             f0_ceil = 500,
                             toFile = FALSE)

# For child speech
f0_child <- trk_straight_f0(wav_file,
                            f0_floor = 150,
                            f0_ceil = 600,
                            toFile = FALSE)
```

### Time Windowing

```r
# Extract F0 from specific time window
f0_window <- trk_straight_f0(wav_file,
                             beginTime = 1.0,   # Start at 1 second
                             endTime = 3.5,     # End at 3.5 seconds
                             toFile = FALSE)
```

### Custom Frame Shift

```r
# Higher temporal resolution (more frames)
f0_highres <- trk_straight_f0(wav_file,
                              frame_shift = 0.5,  # 0.5 ms
                              toFile = FALSE)

# Lower temporal resolution (fewer frames, faster)
f0_lowres <- trk_straight_f0(wav_file,
                             frame_shift = 5.0,  # 5 ms
                             toFile = FALSE)
```

---

## Getting Help

### Documentation

```r
# Function help
?trk_straight_f0
?trk_straight_spec
?straight_synth
?straight_pipeline
?install_legacy_straight

# Package information
?superassp
```

### Module Information

```r
# Detailed status
straight_info()

# Check Python modules
library(reticulate)
py_config()
```

### Support

- **GitHub Issues**: https://github.com/humlab-speech/superassp/issues
- **Documentation**: See `STRAIGHT_INTEGRATION_SUMMARY.md`
- **Python Module**: See `inst/python/legacy_STRAIGHT/README.md`

---

## Next Steps

1. ✅ Install and verify: `install_legacy_straight(install_numba = TRUE)`
2. ✅ Try basic F0 extraction: `trk_straight_f0("speech.wav", toFile = FALSE)`
3. ✅ Experiment with synthesis: `straight_pipeline("speech.wav", "output.wav")`
4. 🔄 Integrate into your workflow
5. 🔄 Provide feedback and report issues

---

**Happy vocodinng with STRAIGHT!**

For more information, see:
- `STRAIGHT_INTEGRATION_SUMMARY.md` - Complete integration details
- `inst/python/legacy_STRAIGHT/README.md` - Python module documentation
- `R_INTEGRATION_COMPLETE.md` - Integration completion report

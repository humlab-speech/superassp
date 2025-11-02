# Voxit Quick Start Guide

## Installation

```r
library(superassp)

# Install voxit module with optimizations
install_voxit(install_numba = TRUE)

# Install SAcC pitch tracker (required)
install_sacc(install_numba = TRUE)

# Verify installation
voxit_available()
voxit_info()
```

## Basic Usage

### Single File Analysis

```r
# With pre-computed word alignments
features <- lst_voxit(
  "speech.wav",
  alignmentFiles = "speech_align.csv"
)

# View results
str(features)
# List of 11
#  $ WPM                             : int 145
#  $ pause_count                     : int 12
#  $ long_pause_count                : int 1
#  $ average_pause_length            : num 0.45
#  $ average_pause_rate              : num 0.38
#  $ rhythmic_complexity_of_pauses   : num 87.3
#  $ average_pitch                   : num 180.5
#  $ pitch_range                     : num 1.2
#  $ pitch_speed                     : num 0.15
#  $ pitch_acceleration              : num 0.08
#  $ pitch_entropy                   : num 3.45

# Access individual features
print(paste("Speaking rate:", features$WPM, "words/minute"))
print(paste("Average pitch:", round(features$average_pitch, 1), "Hz"))
print(paste("Pitch range:", round(features$pitch_range, 2), "octaves"))
```

### Time Windowing

```r
# Analyze specific segment
features <- lst_voxit(
  "speech.wav",
  alignmentFiles = "speech_align.csv",
  beginTime = 1.0,    # Start at 1 second
  endTime = 5.0,      # End at 5 seconds
  minF = 75,          # Adjust F0 range if needed
  maxF = 500
)
```

### Batch Processing

```r
# Process multiple files
audio_files <- c("audio1.wav", "audio2.wav", "audio3.wav")
align_files <- c("align1.csv", "align2.csv", "align3.csv")

results <- lst_voxit(
  audio_files,
  alignmentFiles = align_files,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = 4
)

# Results is a named list
names(results)
# [1] "audio1.wav" "audio2.wav" "audio3.wav"

# Access individual results
results$audio1.wav$WPM
results$audio2.wav$average_pitch
```

### Convert to Data Frame

```r
library(tidyverse)

# Single file
df_single <- as.data.frame(t(unlist(features)))

# Multiple files
df_batch <- results %>%
  map_dfr(~ as.data.frame(t(unlist(.))), .id = "file")

# View as tibble
df_batch %>% as_tibble()
# A tibble: 3 × 12
#   file      WPM pause_count long_pause_count average_pause_length average_pause_rate rhythmic_complexity_of_pauses average_pitch pitch_range pitch_speed pitch_acceleration pitch_entropy
#   <chr>   <dbl>       <dbl>            <dbl>                <dbl>              <dbl>                         <dbl>         <dbl>       <dbl>       <dbl>              <dbl>         <dbl>
# 1 audio1…   145          12                1                 0.45               0.38                          87.3          180.        1.2        0.15               0.08          3.45
# 2 audio2…   132          15                0                 0.52               0.42                          82.1          195.        1.4        0.18               0.10          3.68
# 3 audio3…   158           9                2                 0.38               0.31                          91.2          172.        1.1        0.12               0.06          3.22
```

## Alignment File Format

Create alignment CSV files with these columns:

```csv
word,case,start,end
hello,success,0.0,0.5
there,success,0.6,1.2
how,success,1.5,1.8
are,success,1.9,2.2
you,success,2.3,2.8
```

Required columns:
- `word`: Word text
- `start`: Start time in seconds
- `end`: End time in seconds

Optional columns:
- `case`: Category (e.g., "success", "[noise]")

## Feature Descriptions

### Temporal Features

- **WPM**: Words per minute (speaking rate)
  - Total words divided by duration, normalized to 60 seconds
  
- **pause_count**: Number of pauses between 100ms and 3000ms
  - Excludes very short gaps (<100ms) and very long pauses (>3s)
  
- **long_pause_count**: Number of pauses exceeding 3 seconds
  - Captures disfluencies and long hesitations
  
- **average_pause_length**: Mean duration of valid pauses (seconds)
  - Computed over pauses in 100-3000ms range
  
- **average_pause_rate**: Pauses per second
  - Normalized pause frequency
  
- **rhythmic_complexity_of_pauses**: Normalized Lempel-Ziv complexity (%)
  - Measures predictability of speech/pause rhythm
  - Higher = more complex/unpredictable rhythm
  - Sampled at 100 Hz (every 10ms)
  - Values typically range 50-100%

### Pitch Features

- **average_pitch**: Mean F0 in Hz
  - Computed over voiced regions only
  
- **pitch_range**: F0 range in octaves
  - Difference between highest and lowest pitch (log scale)
  
- **pitch_speed**: F0 velocity in octaves/second
  - Signed directionless: `mean(abs(velocity)) * sign(mean(velocity))`
  - Measures how fast pitch changes
  
- **pitch_acceleration**: F0 acceleration in octaves/second²
  - Signed directionless: `mean(abs(accel)) * sign(mean(accel))`
  - Measures pitch change dynamics
  - Note: Savitzky-Golay smoothed before differentiation
  
- **pitch_entropy**: Shannon entropy of F0 distribution (bits)
  - Computed over 25-bin histogram (±1 octave from mean)
  - Higher = more varied/unpredictable pitch

## Performance Tips

### 1. Use Numba for Speed

```r
# Install with numba (2-3x faster)
install_voxit(install_numba = TRUE)

# Check if available
voxit_info()$numba  # Should be TRUE
```

### 2. Use Cython for Maximum Speed

```r
# Install with cython (3-5x faster, requires C compiler)
install_voxit(compile_cython = TRUE)

# Check if available  
voxit_info()$cython  # Should be TRUE
```

### 3. Parallel Processing for Batches

```r
# Use all cores minus one
results <- lst_voxit(
  files,
  alignmentFiles = aligns,
  parallel = TRUE,
  n_cores = parallel::detectCores() - 1
)
```

## Troubleshooting

### Module Not Found

```r
if (!voxit_available()) {
  install_voxit(install_numba = TRUE)
}
```

### SAcC Not Available

```r
if (!sacc_available()) {
  install_sacc(install_numba = TRUE)
}
```

### Alignment File Issues

```r
# Check CSV format
align <- readr::read_csv("speech_align.csv")
str(align)

# Required columns: word, start, end
all(c("word", "start", "end") %in% names(align))
```

### Audio Loading Problems

```r
# Test audio loading separately
audio <- av::read_audio_bin("speech.wav")
attr(audio, "sample_rate")  # Should show sample rate

# Try different formats
av::av_media_info("speech.wav")
```

## Example Workflow

Complete analysis pipeline:

```r
library(superassp)
library(tidyverse)

# Setup
install_voxit(install_numba = TRUE)
install_sacc(install_numba = TRUE)

# Define files
audio_dir <- "path/to/audio"
audio_files <- list.files(audio_dir, pattern = "\\.wav$", full.names = TRUE)
align_files <- list.files(audio_dir, pattern = "_align\\.csv$", full.names = TRUE)

# Process all files
results <- lst_voxit(
  audio_files,
  alignmentFiles = align_files,
  minF = 75,
  maxF = 500,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = 8
)

# Convert to data frame
df <- results %>%
  map_dfr(~ as.data.frame(t(unlist(.))), .id = "file") %>%
  as_tibble()

# Add metadata
df <- df %>%
  mutate(
    speaker = str_extract(file, "S\\d+"),
    condition = str_extract(file, "cond[AB]")
  )

# Analyze
df %>%
  group_by(condition) %>%
  summarise(
    mean_WPM = mean(WPM),
    mean_pitch = mean(average_pitch),
    mean_complexity = mean(rhythmic_complexity_of_pauses)
  )

# Visualize
ggplot(df, aes(x = WPM, y = average_pitch, color = condition)) +
  geom_point() +
  theme_minimal() +
  labs(title = "Speaking Rate vs Pitch by Condition")
```

## References

- Voxit toolbox: Voice and articulation complexity measures
- SAcC: Ellis & Weiss (2010) - Pitch tracking via harmonic model
- Lempel-Ziv complexity: Measure of sequence complexity
- Savitzky-Golay filter: Smoothing before differentiation

## See Also

- `?lst_voxit` - Full function documentation
- `?install_voxit` - Installation options
- `?voxit_info` - Check module status
- `VOXIT_INTEGRATION_SUMMARY.md` - Technical details

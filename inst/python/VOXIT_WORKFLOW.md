# Voxit Complete Workflow - From Installation to Analysis

This guide demonstrates the complete workflow for using Voxit in superassp.

## 1. Installation

```r
# Install or load superassp
library(superassp)

# Install Voxit with optimizations (recommended)
install_voxit(install_numba = TRUE)

# Install SAcC for pitch tracking
install_sacc(install_numba = TRUE)

# Verify installation
voxit_available()  # Should return TRUE
sacc_available()   # Should return TRUE

# Check optimization status
info <- voxit_info()
print(info)
# $available: TRUE
# $version: "1.0.0"
# $optimized: TRUE
# $numba: TRUE
# $features: 11
```

## 2. Prepare Input Data

### Word Alignments (Required)

Create CSV files with word timing information:

```csv
word,case,start,end
hello,success,0.00,0.45
there,success,0.55,1.05
how,success,1.20,1.45
are,success,1.50,1.75
you,success,1.80,2.25
today,success,2.30,2.85
```

**Required columns:**
- `word` - Word text
- `start` - Start time (seconds)
- `end` - End time (seconds)

**Optional columns:**
- `case` - Status (e.g., "success", "[noise]")

### Audio Files

Any format supported by av package:
- WAV, MP3, MP4, M4A
- FLAC, OGG, AAC
- Video files (audio extracted)

## 3. Single File Analysis

### Basic Usage

```r
# Analyze single file
features <- lst_voxit(
  listOfFiles = "speech.wav",
  alignmentFiles = "speech_align.csv"
)

# View all features
print(features)
# $WPM: 145
# $pause_count: 12
# $long_pause_count: 1
# $average_pause_length: 0.45
# $average_pause_rate: 0.38
# $rhythmic_complexity_of_pauses: 87.3
# $average_pitch: 180.5
# $pitch_range: 1.2
# $pitch_speed: 0.15
# $pitch_acceleration: 0.08
# $pitch_entropy: 3.45
```

### With Time Windowing

```r
# Analyze specific segment (1-5 seconds)
features <- lst_voxit(
  "speech.wav",
  alignmentFiles = "speech_align.csv",
  beginTime = 1.0,
  endTime = 5.0
)
```

### With Custom F0 Range

```r
# Adjust for different voice types
# Female voice: higher F0 range
features <- lst_voxit(
  "female_speech.wav",
  alignmentFiles = "female_align.csv",
  minF = 120,   # Default: 60 Hz
  maxF = 600    # Default: 600 Hz
)

# Male voice: lower F0 range
features <- lst_voxit(
  "male_speech.wav",
  alignmentFiles = "male_align.csv",
  minF = 60,
  maxF = 300
)
```

## 4. Batch Processing

### Process Multiple Files

```r
# Define file lists
audio_dir <- "path/to/audio"
audio_files <- list.files(audio_dir, pattern = "\\.wav$", full.names = TRUE)
align_files <- list.files(audio_dir, pattern = "_align\\.csv$", full.names = TRUE)

# Process all files in parallel
results <- lst_voxit(
  listOfFiles = audio_files,
  alignmentFiles = align_files,
  verbose = TRUE,
  parallel = TRUE,
  n_cores = 4  # Use 4 cores
)

# Results is a named list
names(results)
# [1] "audio1.wav" "audio2.wav" "audio3.wav" ...
```

### Progress Monitoring

```r
# With verbose output
results <- lst_voxit(
  audio_files,
  alignmentFiles = align_files,
  verbose = TRUE,      # Show progress
  parallel = TRUE,
  n_cores = 8
)

# Output:
# ℹ Processing 20 files using 8 cores
# ✔ Processing files
```

## 5. Data Analysis

### Convert to Data Frame

```r
library(tidyverse)

# Single file - convert to data frame
df_single <- as.data.frame(t(unlist(features)))
print(df_single)

# Multiple files - create tidy data frame
df_batch <- results %>%
  map_dfr(~ as.data.frame(t(unlist(.))), .id = "file")

# View
df_batch %>% as_tibble()
# A tibble: 20 × 12
#   file      WPM pause_count long_pause_count average_pause_length ...
#   <chr>   <dbl>       <dbl>            <dbl>                <dbl> ...
# 1 audio1…   145          12                1                 0.45 ...
# 2 audio2…   132          15                0                 0.52 ...
```

### Add Metadata

```r
# Extract information from filenames
df_batch <- df_batch %>%
  mutate(
    speaker = str_extract(file, "S\\d+"),
    condition = str_extract(file, "cond[AB]"),
    trial = str_extract(file, "trial\\d+")
  )
```

### Descriptive Statistics

```r
# Overall statistics
df_batch %>%
  summarise(
    mean_WPM = mean(WPM),
    sd_WPM = sd(WPM),
    mean_pitch = mean(average_pitch),
    sd_pitch = sd(average_pitch),
    mean_complexity = mean(rhythmic_complexity_of_pauses),
    sd_complexity = sd(rhythmic_complexity_of_pauses)
  )

# By condition
df_batch %>%
  group_by(condition) %>%
  summarise(
    n = n(),
    mean_WPM = mean(WPM),
    mean_pitch = mean(average_pitch),
    mean_complexity = mean(rhythmic_complexity_of_pauses)
  )
```

## 6. Visualization

### Speaking Rate

```r
library(ggplot2)

# Speaking rate distribution
ggplot(df_batch, aes(x = WPM)) +
  geom_histogram(bins = 20, fill = "steelblue", alpha = 0.7) +
  labs(title = "Distribution of Speaking Rate",
       x = "Words Per Minute",
       y = "Count") +
  theme_minimal()

# By condition
ggplot(df_batch, aes(x = condition, y = WPM, fill = condition)) +
  geom_boxplot() +
  labs(title = "Speaking Rate by Condition") +
  theme_minimal()
```

### Pitch Dynamics

```r
# Pitch vs speaking rate
ggplot(df_batch, aes(x = WPM, y = average_pitch, color = condition)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Speaking Rate vs Average Pitch",
       x = "Words Per Minute",
       y = "Average Pitch (Hz)") +
  theme_minimal()

# Pitch range distribution
ggplot(df_batch, aes(x = pitch_range)) +
  geom_density(fill = "coral", alpha = 0.5) +
  labs(title = "Distribution of Pitch Range",
       x = "Pitch Range (octaves)") +
  theme_minimal()
```

### Rhythmic Complexity

```r
# Complexity vs pause rate
ggplot(df_batch, aes(x = average_pause_rate, 
                     y = rhythmic_complexity_of_pauses,
                     color = condition)) +
  geom_point(size = 3) +
  labs(title = "Rhythmic Complexity vs Pause Rate",
       x = "Average Pause Rate (pauses/sec)",
       y = "Rhythmic Complexity (%)") +
  theme_minimal()
```

### Multi-Feature Comparison

```r
# Normalize features for comparison
df_normalized <- df_batch %>%
  mutate(across(WPM:pitch_entropy, scale)) %>%
  select(file, condition, WPM, average_pitch, pitch_range, 
         rhythmic_complexity_of_pauses)

# Convert to long format
df_long <- df_normalized %>%
  pivot_longer(cols = WPM:rhythmic_complexity_of_pauses,
               names_to = "feature",
               values_to = "value")

# Faceted plot
ggplot(df_long, aes(x = condition, y = value, fill = condition)) +
  geom_boxplot() +
  facet_wrap(~ feature, scales = "free_y") +
  labs(title = "Normalized Features by Condition") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

## 7. Statistical Analysis

### T-tests

```r
# Compare conditions
t.test(WPM ~ condition, data = df_batch)
t.test(average_pitch ~ condition, data = df_batch)
t.test(rhythmic_complexity_of_pauses ~ condition, data = df_batch)
```

### Correlations

```r
# Feature correlations
library(corrplot)

cor_matrix <- df_batch %>%
  select(WPM:pitch_entropy) %>%
  cor(use = "pairwise.complete.obs")

corrplot(cor_matrix, method = "color", type = "upper",
         addCoef.col = "black", number.cex = 0.7)
```

### Linear Models

```r
# Predict speaking rate from pitch features
model <- lm(WPM ~ average_pitch + pitch_range + pitch_entropy, 
            data = df_batch)
summary(model)

# Predict rhythmic complexity
model2 <- lm(rhythmic_complexity_of_pauses ~ 
             average_pause_rate + long_pause_count + WPM,
             data = df_batch)
summary(model2)
```

## 8. Export Results

### Save to CSV

```r
# Save data frame
write_csv(df_batch, "voxit_results.csv")

# Save with timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
write_csv(df_batch, paste0("voxit_results_", timestamp, ".csv"))
```

### Save to RDS

```r
# Save R object
saveRDS(results, "voxit_results.rds")

# Load later
results <- readRDS("voxit_results.rds")
```

### Generate Report

```r
# Create summary report
library(knitr)
library(rmarkdown)

# Simple markdown report
report <- "
# Voxit Analysis Report

## Summary Statistics
```{r echo=FALSE}
df_batch %>% 
  summarise(across(WPM:pitch_entropy, 
                  list(mean = mean, sd = sd)))
```

## Plots
```{r echo=FALSE}
ggplot(df_batch, aes(x = WPM, y = average_pitch)) +
  geom_point() + theme_minimal()
```
"

writeLines(report, "report.Rmd")
render("report.Rmd", output_format = "html_document")
```

## 9. Troubleshooting

### Module Not Found

```r
# Check installation
if (!voxit_available()) {
  cat("Installing voxit...\n")
  install_voxit(install_numba = TRUE)
}
```

### Alignment File Errors

```r
# Validate alignment file
validate_alignment <- function(align_file) {
  df <- read_csv(align_file, show_col_types = FALSE)
  
  # Check required columns
  required <- c("word", "start", "end")
  missing <- setdiff(required, names(df))
  
  if (length(missing) > 0) {
    stop("Missing columns: ", paste(missing, collapse = ", "))
  }
  
  # Check data types
  if (!is.numeric(df$start) || !is.numeric(df$end)) {
    stop("Start and end times must be numeric")
  }
  
  # Check ordering
  if (any(df$end < df$start)) {
    stop("End time before start time detected")
  }
  
  cat("✓ Alignment file is valid\n")
  return(TRUE)
}

# Use it
validate_alignment("speech_align.csv")
```

### Performance Issues

```r
# Check optimization status
info <- voxit_info()
if (!info$numba) {
  cat("Numba not available - installing for 2-3x speedup...\n")
  install_voxit(install_numba = TRUE)
}

# Reduce parallel cores if memory issues
results <- lst_voxit(
  files,
  alignmentFiles = aligns,
  n_cores = 2  # Use fewer cores
)
```

## 10. Complete Example

```r
library(superassp)
library(tidyverse)

# Setup
install_voxit(install_numba = TRUE)
install_sacc(install_numba = TRUE)

# Process files
audio_files <- list.files("audio", pattern = "\\.wav$", full.names = TRUE)
align_files <- list.files("audio", pattern = "_align\\.csv$", full.names = TRUE)

results <- lst_voxit(
  audio_files,
  alignmentFiles = align_files,
  minF = 75,
  maxF = 500,
  parallel = TRUE,
  n_cores = 8
)

# Analyze
df <- results %>%
  map_dfr(~ as.data.frame(t(unlist(.))), .id = "file") %>%
  mutate(
    speaker = str_extract(file, "S\\d+"),
    condition = str_extract(file, "cond[AB]")
  )

# Visualize
ggplot(df, aes(x = WPM, y = average_pitch, color = condition)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(title = "Speaking Rate vs Pitch")

# Export
write_csv(df, "voxit_analysis.csv")

cat("✓ Analysis complete!\n")
cat(sprintf("Processed %d files\n", nrow(df)))
cat(sprintf("Mean WPM: %.1f\n", mean(df$WPM)))
cat(sprintf("Mean pitch: %.1f Hz\n", mean(df$average_pitch)))
```

## Summary

This workflow covers:
1. ✅ Installation with optimizations
2. ✅ Input data preparation
3. ✅ Single and batch processing
4. ✅ Data frame conversion
5. ✅ Descriptive statistics
6. ✅ Visualization
7. ✅ Statistical analysis
8. ✅ Results export
9. ✅ Troubleshooting
10. ✅ Complete example

For more details, see:
- `?lst_voxit` - Function documentation
- `VOXIT_QUICKSTART.md` - Quick start guide
- `VOXIT_INTEGRATION_SUMMARY.md` - Technical details

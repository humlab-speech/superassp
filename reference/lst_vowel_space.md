# Vowel space analysis (F1×F2 area ratio)

Computes vowel space area from formant data using k-means clustering and
reference vowel comparison. The vowel space ratio is computed as the
ratio of the actual vowel space convex hull area to a reference vowel
space area for the same gender.

## Usage

``` r
lst_vowel_space(
  formant_data,
  gender = 1,
  mode = "triangle",
  scaling = FALSE,
  plot_formants = FALSE,
  return_jstf = FALSE
)
```

## Arguments

- formant_data:

  Data frame or matrix with columns F1, F2 (or first two columns). F1
  and F2 should be in Hz.

- gender:

  Speaker gender for reference vowel selection: 0 = female, 1 = male
  (default), 2 = child

- mode:

  Vowel space computation mode: "triangle" = use 3 corner vowels
  (default, more robust) "polygon" = use 4 vowels (may be less stable)

- scaling:

  Apply Bark frequency scaling (default: FALSE)

- plot_formants:

  Show scatter plot of formant points (default: FALSE)

- return_jstf:

  Logical. Return JsonTrackObj instead of data.frame? Default FALSE.
  When both toFile and return_jstf are TRUE, the file is written AND the
  object returned.

## Value

List with elements:

- `vowel_space_ratio`: Numeric, ratio of actual to reference vowel space
  area (0–2+ range)

- `centroids`: Matrix of estimated vowel formant centroids from k-means
  clustering

- `n_frames`: Number of frames used in computation

## Details

**Algorithm**:

1.  Filter formant data by gender-specific frequency ranges to identify
    candidate vowels

2.  Run k-means clustering using reference vowel centroids as
    initialization

3.  Match estimated centroids to reference vowels using Mahalanobis
    distance

4.  Compute convex hull of matched vowel points (actual vowel space)

5.  Compute convex hull of reference vowels

6.  Ratio = actual area / reference area

**Interpretation**:

- Ratio = 1.0: perfect agreement with reference vowel space

- Ratio \< 1.0: speaker has reduced vowel space (dysarthria, age, etc.)

- Ratio \> 1.0: speaker has expanded vowel space (hyperarticulation,
  singing)

- NaN/error: insufficient vowel tokens (\<1000 frames) in valid range

**Reference vowel sets**:

- Triangle: /i/, /a/, /u/ (corner vowels)

- Polygon: /i/, /\\\varepsilon\\/, /a/, /u/ (4-vowel system)

**Minimum frame requirement**: 1000 frames must be within
gender-specific frequency bounds, or function returns 0 ratio
(insufficient data).

## Examples

``` r
if (FALSE) { # \dontrun{
# Compute vowel space from file batch
files <- c("speaker1.wav", "speaker2.wav")
formants <- trk_formant_burg(files, toFile = FALSE)

# Extract F1-F2 and compute vowel space
formant_df <- data.frame(F1 = formants$fm1, F2 = formants$fm2)
vs <- lst_vowel_space(formant_df, gender = 1)
cat("Vowel space ratio:", vs$vowel_space_ratio, "\n")

# For female speaker
vs_female <- lst_vowel_space(formant_df, gender = 0)

# Using Bark scaling
vs_bark <- lst_vowel_space(formant_df, gender = 1, scaling = TRUE)
} # }
```

# Create ggplot with automatic track labels

Creates a ggplot2 plot with automatic axis labels derived from track
names. This is a convenience wrapper around
[`ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html) that
automatically extracts and applies appropriate labels for acoustic track
data.

## Usage

``` r
ggtrack(
  data,
  mapping = ggplot2::aes(),
  ...,
  full_labels = FALSE,
  use_subscripts = TRUE
)
```

## Arguments

- data:

  data.frame or tibble. Data from
  [`as.data.frame.AsspDataObj()`](https://humlab-speech.github.io/superassp/reference/AsspDataObj.md).

- mapping:

  ggplot2::aes() specification. Aesthetic mappings.

- ...:

  Additional ggplot2 layers to add (geoms, scales, themes, etc.).

- full_labels:

  Logical. If TRUE, use full descriptive labels. If FALSE (default), use
  short labels suitable for plot axes.

- use_subscripts:

  Logical. If TRUE (default), use plotmath expressions with subscripts
  (fo → f₀, F1 → F₁). If FALSE, use plain text.

## Value

A ggplot object with automatic axis labels.

## Details

This function simplifies plotting of acoustic track data by
automatically generating appropriate axis labels based on column names.
It extracts the x and y variables from the aesthetic mapping and applies
labels using
[`get_track_label()`](https://humlab-speech.github.io/superassp/reference/get_track_label.md).

**Short labels** (full_labels = FALSE, default):

- "fo \[Hz\]", "F1 \[Hz\]", "H1-H2c \[dB\]"

- Concise, suitable for most plots

**Full labels** (full_labels = TRUE):

- "Frequency of oscillation \[Hz\]"

- "First formant frequency \[Hz\]"

- "H1-H2 corrected for formants \[dB\]"

- Descriptive, suitable for publications

**Additional layers** can be added using the `...` argument or by adding
to the returned ggplot object with `+`.

## See also

- [`as.data.frame.AsspDataObj()`](https://humlab-speech.github.io/superassp/reference/AsspDataObj.md)
  for creating the data frame

- [`get_track_label()`](https://humlab-speech.github.io/superassp/reference/get_track_label.md)
  for label extraction

- [`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
  for the underlying plotting function

## Examples

``` r
if (FALSE) { # \dontrun{
library(ggplot2)

# Get pitch track
pitch <- wrassp::ksvF0("audio.wav", toFile = FALSE)
df <- as.data.frame(pitch)

# Quick plot with automatic labels
ggtrack(df, aes(x = frame_time, y = fo_Hz)) +
  geom_line()
# Y-axis automatically labeled: "fo \[Hz\]"

# With full descriptive labels
ggtrack(df, aes(x = frame_time, y = fo_Hz), full_labels = TRUE) +
  geom_line() +
  theme_minimal()
# Y-axis: "Frequency of oscillation \[Hz\]"

# Adding layers via ...
fms <- wrassp::forest("vowel.wav", toFile = FALSE)
df_fms <- as.data.frame(fms)

ggtrack(df_fms, aes(x = F2_Hz, y = F1_Hz),
        geom_point(alpha = 0.5),
        scale_x_reverse(),
        scale_y_reverse(),
        theme_minimal())
# Automatic labels: "F1 \[Hz\]", "F2 \[Hz\]"

# Can also add layers with +
ggtrack(df, aes(x = frame_time, y = fo_Hz)) +
  geom_line(color = "steelblue") +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  theme_bw()
} # }
```

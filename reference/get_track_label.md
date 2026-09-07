# Get track label for plotting

Extracts the display label for a track column, either from stored
attributes or by intelligent inference.

## Usage

``` r
get_track_label(df, col, full = FALSE)
```

## Arguments

- df:

  data.frame or tibble. Data frame from
  [`as.data.frame.AsspDataObj()`](https://humlab-speech.github.io/superassp/reference/AsspDataObj.md).

- col:

  Character. Column name.

- full:

  Logical. If TRUE, return full descriptive label. If FALSE (default),
  return short label suitable for plot axes.

## Value

Character. Display label for the track.

## Details

This function first checks the `track_labels` or `track_descriptions`
attributes stored by
[`as.data.frame.AsspDataObj()`](https://humlab-speech.github.io/superassp/reference/AsspDataObj.md).
If not found, it intelligently infers the label from the column name.

**Short labels** (full = FALSE):

- "fo \[Hz\]", "F1 \[Hz\]", "H1-H2c \[dB\]"

- Suitable for plot axes

**Full labels** (full = TRUE):

- "Frequency of oscillation \[Hz\]"

- "First formant frequency \[Hz\]"

- "H1-H2 corrected for formants \[dB\]"

- Suitable for documentation, papers

## Examples

``` r
if (FALSE) { # \dontrun{
fms <- wrassp::forest("audio.wav", toFile = FALSE)
df <- as.data.frame(fms)

# Short label
get_track_label(df, "F1_Hz")
# [1] "F1 [Hz]"

# Full descriptive label
get_track_label(df, "F1_Hz", full = TRUE)
# [1] "First formant frequency [Hz]"

# Works even without attributes (inference)
df2 <- data.frame(frame_time = 1:10, fo_Hz = rnorm(10, 120, 10))
get_track_label(df2, "fo_Hz")
# [1] "fo [Hz]"
} # }
```

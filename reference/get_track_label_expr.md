# Get track label as expression for plotting

Extended version of
[`get_track_label()`](https://humlab-speech.github.io/superassp/reference/get_track_label.md)
that returns plotmath expressions with subscripts when
`use_subscripts = TRUE`.

## Usage

``` r
get_track_label_expr(df, col, full = FALSE, use_subscripts = TRUE)
```

## Arguments

- df:

  data.frame or tibble. Data frame from
  [`as.data.frame.AsspDataObj()`](https://humlab-speech.github.io/superassp/reference/AsspDataObj.md).

- col:

  Character. Column name.

- full:

  Logical. If TRUE, return full descriptive label. Default: FALSE.

- use_subscripts:

  Logical. If TRUE, return plotmath expression with subscripts. If
  FALSE, return plain text. Default: TRUE.

## Value

Expression object (if use_subscripts = TRUE) or character string.

## Details

This function extends
[`get_track_label()`](https://humlab-speech.github.io/superassp/reference/get_track_label.md)
with the option to generate plotmath expressions for scientific
subscripts in ggplot2.

**With subscripts** (use_subscripts = TRUE):

- fo → f\[o\] (f subscript o)

- F1 → F\[1\] (F subscript 1)

- H1-H2c → H\[1\]-H\[2c\] (proper subscripts)

**Without subscripts** (use_subscripts = FALSE):

- Same as
  [`get_track_label()`](https://humlab-speech.github.io/superassp/reference/get_track_label.md)

- Returns: "fo \[Hz\]", "F1 \[Hz\]", "H1-H2c \[dB\]"

For full descriptive labels (full = TRUE), subscripts are not used as
the full text doesn't need them (e.g., "Frequency of oscillation
\[Hz\]").

## Examples

``` r
if (FALSE) { # \dontrun{
fms <- wrassp::forest("audio.wav", toFile = FALSE)
df <- as.data.frame(fms)

# Get expression with subscript
expr <- get_track_label_expr(df, "F1_Hz")
# Returns: expression(F\[1\]~"\[Hz\]")

# Use in ggplot
library(ggplot2)
ggplot(df, aes(x = frame_time, y = F1_Hz)) +
  geom_line() +
  labs(y = get_track_label_expr(df, "F1_Hz"))
# Y-axis shows: F1 \[Hz\] (rendered with a subscript)

# Without subscripts
label <- get_track_label_expr(df, "F1_Hz", use_subscripts = FALSE)
# Returns: "F1 \[Hz\]"
} # }
```

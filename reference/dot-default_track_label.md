# Default y-axis label for a track table

Only a table holding exactly one track is unambiguous; anything else
keeps ggplot2's own default label. The sampled waveform track ("audio")
is amplitude, not a measured quantity with a label of its own.

## Usage

``` r
.default_track_label(data, full = FALSE, use_subscripts = TRUE)
```

## Arguments

- data:

  data.frame or NULL. Wide track table.

- full:

  Logical. Use the full descriptive label.

- use_subscripts:

  Logical. Use plotmath subscripts in the short label.

## Value

The label of the single track, a waiver otherwise.

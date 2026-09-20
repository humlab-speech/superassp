# Long time x frequency table from a spectral track

Long time x frequency table from a spectral track

## Usage

``` r
.assp_spec_data(
  data,
  tracks = NULL,
  time = "frame_time",
  freq_hz_per_bin = NULL,
  na.zeros = FALSE
)
```

## Arguments

- data:

  AsspDataObj, data.frame, or long track table.

- tracks:

  Character or NULL. Spectral track (band) to use.

- time:

  Character. Preferred name of the time column.

- freq_hz_per_bin:

  Numeric or NULL. Spacing of the coefficients in Hz, overriding the
  value derived from `origFreq`.

- na.zeros:

  Logical. Convert stored zeros to `NA`.

## Value

A data.frame with columns `frame_time`, `freq` and `value`.

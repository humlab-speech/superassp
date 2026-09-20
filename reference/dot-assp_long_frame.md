# Assemble the long track table

Assemble the long track table

## Usage

``` r
.assp_long_frame(
  frame_time,
  values,
  track,
  band,
  bin,
  na.zeros,
  origFreq = NULL
)
```

## Arguments

- frame_time:

  Numeric vector. Time of each record, in seconds.

- values:

  List of numeric vectors, one per track column.

- track, band:

  Character vectors. Column name and band name per track.

- bin:

  Integer vector. Coefficient index per track (`NA` if not indexed).

- na.zeros:

  Logical. Convert stored zeros to `NA`.

- origFreq:

  Numeric or NULL. Sample rate the tracks were computed from.

## Value

A data.frame with columns `frame_time`, `value`, `track`, `band` and
`bin`.

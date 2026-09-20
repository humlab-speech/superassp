# Long track table from a wide track table

Long track table from a wide track table

## Usage

``` r
.assp_wide_to_long(data, tracks = NULL, time = "frame_time", na.zeros = FALSE)
```

## Arguments

- data:

  data.frame. Wide table, as written by
  [`as.data.frame.AsspDataObj()`](https://humlab-speech.github.io/superassp/reference/AsspDataObj.md).

- tracks:

  Character or NULL. Tracks to keep.

- time:

  Character. Preferred name of the time column.

- na.zeros:

  Logical. Convert stored zeros to `NA`.

## Value

A data.frame as described in
[`.assp_long_frame()`](https://humlab-speech.github.io/superassp/reference/dot-assp_long_frame.md).

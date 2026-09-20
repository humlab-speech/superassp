# Long track table from any supported input

Long track table from any supported input

## Usage

``` r
.assp_long_data(data, tracks = NULL, time = "frame_time", na.zeros = FALSE)
```

## Arguments

- data:

  AsspDataObj, data.frame, or an already long track table.

- tracks:

  Character or NULL. Tracks to keep.

- time:

  Character. Preferred name of the time column.

- na.zeros:

  Logical. Convert stored zeros to `NA`.

## Value

A data.frame with columns `frame_time`, `value`, `track`, `band` and
`bin`, carrying `origFreq` as an attribute when it is known.

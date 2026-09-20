# Long track table from an AsspDataObj

Long track table from an AsspDataObj

## Usage

``` r
.assp_obj_to_long(obj, tracks = NULL, na.zeros = FALSE)
```

## Arguments

- obj:

  AsspDataObj.

- tracks:

  Character or NULL. Tracks to keep.

- na.zeros:

  Logical. Convert stored zeros to `NA`.

## Value

A data.frame as described in
[`.assp_long_frame()`](https://humlab-speech.github.io/superassp/reference/dot-assp_long_frame.md).

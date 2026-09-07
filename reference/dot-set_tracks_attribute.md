# Set tracks attribute on AsspDataObj result

Internal helper to propagate the tracks attribute from a function to its
returned AsspDataObj object(s). This enables proper template expansion
in
[`as.data.frame.AsspDataObj()`](https://humlab-speech.github.io/superassp/reference/AsspDataObj.md).

## Usage

``` r
.set_tracks_attribute(result, func, n_files = 1)
```

## Arguments

- result:

  AsspDataObj or list of AsspDataObj objects

- func:

  Function object whose tracks attribute should be copied

- n_files:

  Number of files processed (for list handling)

## Value

The result with tracks attribute set

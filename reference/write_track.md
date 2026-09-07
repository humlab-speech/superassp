# Write Track to File

Counterpart to
[`read_track()`](https://humlab-speech.github.io/superassp/reference/read_track.md).
Dispatches to
[`write_ssff()`](https://humlab-speech.github.io/superassp/reference/write_ssff.md)
for AsspDataObj or
[`write_jstf()`](https://humlab-speech.github.io/superassp/reference/write_jstf.md)
for JsonTrackObj based on object class.

## Usage

``` r
write_track(obj, file, ...)
```

## Arguments

- obj:

  AsspDataObj or JsonTrackObj to write

- file:

  Output file path

- ...:

  Additional arguments passed to the underlying writer

## Value

Invisibly returns file path

## See also

[`read_track()`](https://humlab-speech.github.io/superassp/reference/read_track.md),
[`write_ssff()`](https://humlab-speech.github.io/superassp/reference/write_ssff.md),
[`write_jstf()`](https://humlab-speech.github.io/superassp/reference/write_jstf.md)

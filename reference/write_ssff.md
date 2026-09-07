# Write an AsspDataObj to an SSFF file

User-facing wrapper around the ASSP C-level writer. Interface mirrors
the legacy `write.AsspDataObj`.

## Usage

``` r
write_ssff(dobj, file = attr(dobj, "filePath"))
```

## Arguments

- dobj:

  An `AsspDataObj`.

- file:

  Output file path. Defaults to the `filePath` attribute of `dobj`.

## Value

Invisibly, the resolved output file path (after `path.expand`). Use this
for chained pipelines. The file written is an SSFF container with the
tracks stored in `dobj` encoded according to
`attr(dobj, "trackFormats")`.

## See also

[`write_jstf`](https://humlab-speech.github.io/superassp/reference/write_jstf.md)
for JSTF (JSON) output.

## Examples

``` r
if (FALSE) { # \dontrun{
wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")
f0  <- trk_pitch_rapt(wav, toFile = FALSE)

out <- tempfile(fileext = ".f0")
write_ssff(f0, file = out)

f0_back <- read_ssff(out)
identical(names(f0), names(f0_back))
} # }
```

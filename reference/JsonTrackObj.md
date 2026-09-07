# JsonTrackObj — JSON Track Format Object

A list-based S3 class representing a JSTF (JSON Speech Track Format)
file in memory: a self-describing container for summary measures
(jitter, shimmer, voice-quality scores, …) with their field schema and
provenance. Produced by `lst_*` functions with `return_jstf = TRUE` or
`toFile = FALSE`, and read back by
[`read_jstf()`](https://humlab-speech.github.io/superassp/reference/read_jstf.md).
Unlike `AsspDataObj` (equally-spaced signal tracks), JSTF holds a small
number of time-bounded *slices*, each a set of named scalar or vector
values.

## Usage

``` r
# S3 method for class 'JsonTrackObj'
print(x, ...)

# S3 method for class 'JsonTrackObj'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'JsonTrackObj'
as_tibble(x, ...)

# S3 method for class 'JsonTrackObj'
summary(object, ...)
```

## Arguments

- x:

  JsonTrackObj

- ...:

  Additional arguments (ignored; present for S3 method signature
  compatibility).

- row.names:

  NULL or character vector of row names

- optional:

  Logical, if TRUE column names are checked for syntactic validity

- object:

  JsonTrackObj

## Methods (by generic)

- `print(JsonTrackObj)`: Print a compact summary of a JsonTrackObj.

- `as.data.frame(JsonTrackObj)`: Convert to a data.frame; each slice
  becomes a row.

- `as_tibble(JsonTrackObj)`: Convert to a tibble with typed columns.
  Requires the tibble package.

- `summary(JsonTrackObj)`: Print a summary of a JsonTrackObj to the
  console.

## Structure

A named list with class `c("JsonTrackObj", "list")`:

- `format`, `version` — always `"JSTF"` and the schema version.

- `function_name` — the `lst_*` function that produced it.

- `file_path`, `sample_rate`, `audio_duration` — source provenance.

- `metadata` — `function_version` and the `parameters` used.

- `field_schema` — named list mapping each value field to its type
  (`"numeric"`, `"numeric_vector"`, `"character"`, …).

- `slices` — list of `{begin_time, end_time, values}` records; `values`
  is a named list keyed by the `field_schema`.

## Inspecting an object

[`print()`](https://rdrr.io/r/base/print.html) shows a compact summary
(format, function, slice count). Read fields directly:
`obj$field_schema` lists available measures, and
`obj$slices[[1]]$values` holds the values for the first slice.
[`sample_rate()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md)
and
[`file_path()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md)
accessors also work. Use
[`write_jstf()`](https://humlab-speech.github.io/superassp/reference/write_jstf.md)
to serialize and
[`read_jstf()`](https://humlab-speech.github.io/superassp/reference/read_jstf.md)
to load.

## See also

[assp_accessors](https://humlab-speech.github.io/superassp/reference/assp_accessors.md)
for accessor generics that work on this class;
[`read_jstf()`](https://humlab-speech.github.io/superassp/reference/read_jstf.md),
[`write_jstf()`](https://humlab-speech.github.io/superassp/reference/write_jstf.md)
for I/O.

## Examples

``` r
if (FALSE) { # \dontrun{
# Produce a JSTF object from a summary function
vq <- lst_voice_report(
  system.file("samples", "sustained", "a1.wav", package = "superassp"),
  toFile = FALSE
)

vq                      # compact summary (print method)
names(vq$field_schema)  # available measures
vq$slices[[1]]$values   # values for the first slice

# Round-trip to disk
f <- tempfile(fileext = ".json")
write_jstf(vq, f)
identical_obj <- read_jstf(f)
} # }
```

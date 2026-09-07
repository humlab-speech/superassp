# Write JSTF Object to File

Writes a JsonTrackObj to a JSTF file using jsonlite.

## Usage

``` r
write_jstf(obj, file, pretty = FALSE, digits = 6, auto_unbox = TRUE)
```

## Arguments

- obj:

  JsonTrackObj to write

- file:

  Output file path

- pretty:

  Logical, pretty-print JSON (default: FALSE for efficiency)

- digits:

  Number of decimal digits for numbers (default: 6)

- auto_unbox:

  Logical, automatically unbox single-element arrays (default: TRUE)

## Value

Invisibly returns file path

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- create_json_track_obj(
  results = list(f0_mean = 150, f0_sd = 20),
  function_name = "lst_example",
  file_path = "audio.wav",
  sample_rate = 16000,
  audio_duration = 5.0
)
write_jstf(obj, "output.jstf")
} # }
```

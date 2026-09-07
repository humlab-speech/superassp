# AsspDataObj — ASSP Data Object

S3 class for in-memory ASSP/SSFF signal data: a set of equally-spaced,
time-aligned tracks (audio samples, or analysis frames such as F0,
formants, RMS). Produced by
[`read_ssff()`](https://humlab-speech.github.io/superassp/reference/read_ssff.md),
[`read_audio()`](https://humlab-speech.github.io/superassp/reference/read_audio.md),
and every `trk_*` function called with `toFile = FALSE`. The layout is
compatible with `emuR` and can be written to disk with
[`write_ssff()`](https://humlab-speech.github.io/superassp/reference/write_ssff.md).

## Usage

``` r
# S3 method for class 'AsspDataObj'
as.data.frame(
  x,
  ...,
  convert_units = TRUE,
  clean_names = TRUE,
  na.zeros = FALSE
)

# S3 method for class 'AsspDataObj'
print(x, ...)

# S3 method for class 'AsspDataObj'
as_tibble(
  x,
  field = NULL,
  beginTime = NULL,
  endTime = NULL,
  na.zeros = TRUE,
  convert_units = TRUE,
  clean_names = TRUE
)

# S3 method for class 'AsspDataObj'
cut(obj, where, n_preceeding, n_following)
```

## Arguments

- x:

  AsspDataObj

- ...:

  additional arguments (ignored)

- convert_units:

  Convert columns with unit labels (e.g., "fo\[Hz\]") to units objects.
  Default: TRUE.

- clean_names:

  Logical. Convert bracket notation to underscore notation. Default:
  TRUE.

- na.zeros:

  Replace zero values with NA. Default: TRUE.

- field:

  Optional field name or index to extract. If NULL, all fields are
  extracted.

- beginTime:

  Start time for the extracted portion in seconds. Default: NULL
  (beginning of signal). Note: uses `beginTime`/`endTime` (seconds)
  matching DSP function conventions, unlike
  [`read_audio()`](https://humlab-speech.github.io/superassp/reference/read_audio.md)
  which uses `begin`/`end`.

- endTime:

  End time for the extracted portion in seconds. Default: NULL (end of
  signal).

- obj:

  AsspDataObj to cut.

- where:

  Relative time 0.0–1.0 where the cutout is centred.

- n_preceeding:

  Max samples to include before the centre.

- n_following:

  Max samples to include after the centre.

## Methods (by generic)

- `as.data.frame(AsspDataObj)`: Convert to a data.frame with template
  expansion, optional clean names and unit assignment.

- `print(AsspDataObj)`: Print a summary of the AsspDataObj (also aliased
  as summary.AsspDataObj).

- `as_tibble(AsspDataObj)`: Convert to a tibble with `times_orig`,
  `times_rel`, `times_norm` and track columns; compatible with emuR
  workflows.

- `cut(AsspDataObj)`: Cut out a portion centred at a relative time point
  (0.0–1.0).

## Structure

An `AsspDataObj` is a named list — one element per track, each a numeric
matrix of `n_records` rows by (channels or coefficients) columns —
carrying metadata as attributes:

- `sampleRate` — frame rate in Hz (sample rate for audio,
  `1000 / windowShift` for analysis tracks).

- `startTime` — time of the first record, in seconds.

- `startRecord`, `endRecord` — record index bounds.

- `trackFormats` — storage type per track (`"INT16"`, `"REAL32"`,
  `"REAL64"`, …).

- `origFreq` — original audio sample rate (analysis tracks only).

- `filePath` — source path.

## Inspecting an object

Use the accessor generics (see
[assp_accessors](https://humlab-speech.github.io/superassp/reference/assp_accessors.md))
rather than reaching into attributes directly:
[`track_names()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md),
[`sample_rate()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md),
[`n_records()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md),
[`signal_duration()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md),
[`start_time()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md),
[`track_formats()`](https://humlab-speech.github.io/superassp/reference/assp_accessors.md).
Track matrices are reached by name, e.g. `obj[["F0"]]`.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) /
`as_tibble()` flatten all tracks into one time-indexed table (a
`frame_time` column plus one column per track/coefficient), which is the
usual bridge to `dplyr`/plotting.

## See also

[assp_accessors](https://humlab-speech.github.io/superassp/reference/assp_accessors.md)
for accessor generics that work on this class;
[`read_audio()`](https://humlab-speech.github.io/superassp/reference/read_audio.md),
[`read_ssff()`](https://humlab-speech.github.io/superassp/reference/read_ssff.md),
[`write_ssff()`](https://humlab-speech.github.io/superassp/reference/write_ssff.md)
for I/O.

## Examples

``` r
wav <- system.file("samples", "sustained", "a1.wav", package = "superassp")

# read -> analyze -> tabulate
rms <- trk_rms(wav, toFile = FALSE, verbose = FALSE)

track_names(rms)      # "RMS[dB]"
#> [1] "RMS[dB]"
n_records(rms)        # number of analysis frames
#> [1] 805
sample_rate(rms)      # frame rate in Hz
#> [1] 199.5475

df <- as.data.frame(rms)
#> Warning: Package 'units' not available. Skipping unit assignment.
head(df)              # frame_time + one column per track
#>    frame_time   RMS_dB
#> 1 0.002505669 26.59310
#> 2 0.007517007 30.23827
#> 3 0.012528345 31.27690
#> 4 0.017539683 30.04947
#> 5 0.022551020 26.95653
#> 6 0.027562358 23.09678
```

# Unified Track Reading Interface

Reads either SSFF or JSTF files transparently based on file extension.

## Usage

``` r
read_track(
  file,
  begin = 0,
  end = 0,
  samples = FALSE,
  validate = TRUE,
  zero_to_na = TRUE,
  tracks = NULL,
  threads = 1L
)
```

## Arguments

- file:

  Path to track file (.f0, .fms, .vat, .vsj, etc.)

- begin:

  Start of region to read (seconds, or samples if `samples = TRUE`).
  Default 0 = file start. For JSTF files, slices whose `begin_time` is
  before this value are excluded.

- end:

  End of region to read (seconds, or samples if `samples = TRUE`).
  Default 0 = file end. For JSTF files, slices whose `begin_time` is
  after this value are excluded.

- samples:

  Logical. If `TRUE`, `begin`/`end` are in samples (divided by
  `sample_rate` in the JSTF header); otherwise in seconds. Default
  `FALSE`.

- validate:

  Logical, validate after reading (default: TRUE, JSTF only)

- zero_to_na:

  Logical, SSFF only. If `TRUE` (default) stored values that are exactly
  `0` are returned as `NA` for every track that is not sampled audio.
  SSFF has no NULL/NA encoding and `0` is its substitute, so this is
  what makes "no value" distinguishable from a measured zero. Use
  [`read_ssff`](https://humlab-speech.github.io/superassp/reference/read_ssff.md)
  for the values exactly as stored.

- tracks:

  Optional character vector of track names to read (SSFF only). `NULL`
  (default) reads every track in the file; unselected tracks are skipped
  without being converted.

- threads:

  Number of threads used to convert large SSFF files (default 1,
  serial). Results are identical regardless of the value.

## Value

AsspDataObj (for SSFF) or JsonTrackObj (for JSTF)

## Examples

``` r
if (FALSE) { # \dontrun{
# Read SSFF pitch track
pitch <- read_track("audio.f0")

# Read SSFF with time windowing
pitch <- read_track("audio.f0", begin = 1.0, end = 2.5)

# Read JSTF voice quality track
vq <- read_track("audio.vat")

# Both can be converted to data.frame
df1 <- as.data.frame(pitch)
df2 <- as.data.frame(vq)
} # }
```

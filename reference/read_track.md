# Unified Track Reading Interface

Reads either SSFF or JSTF files transparently based on file extension.

## Usage

``` r
read_track(file, begin = 0, end = 0, samples = FALSE, validate = TRUE)
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

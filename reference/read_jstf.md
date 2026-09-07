# Read JSTF File

Reads a JSTF file using RcppSimdJson for high-performance parsing,
falling back to jsonlite if RcppSimdJson is unavailable.

## Usage

``` r
read_jstf(file, begin = 0, end = 0, samples = FALSE, validate = TRUE)
```

## Arguments

- file:

  Path to JSTF file.

- begin:

  Start of region to read in seconds (or samples if `samples = TRUE`).
  Default 0 = file start. Slices whose `begin_time` falls within
  `[begin, end]` are retained.

- end:

  End of region to read in seconds (or samples if `samples = TRUE`).
  Default 0 = no upper limit (entire file).

- samples:

  Logical. If `TRUE`, `begin`/`end` are interpreted as sample indices
  and converted to seconds using the file's `sample_rate` field.

- validate:

  Logical, validate after reading (default: TRUE).

## Value

JsonTrackObj (possibly subset to the requested time window)

## Examples

``` r
if (FALSE) { # \dontrun{
obj <- read_jstf("output.jstf")
# Read only slices between 1 s and 3 s
obj_sub <- read_jstf("output.jstf", begin = 1, end = 3)
df  <- as.data.frame(obj)
} # }
```

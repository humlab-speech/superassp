# Assign units to data frame columns

Assigns R units package units to columns based on their unit suffix.

## Usage

``` r
.assign_track_units(df)
```

## Arguments

- df:

  Data frame. Data frame with cleaned column names.

## Value

Data frame with units assigned to appropriate columns.

## Details

Maps common unit suffixes to R units package units:

- `_Hz` → Hz

- `_dB` → dB (dimensionless, but labeled)

- `_pct` → percent (dimensionless)

- `_us` → microseconds

- `_ms` → milliseconds

- `_s` → seconds

- `_Bark` → Bark (dimensionless)

- `_mel` → mel (dimensionless)

- `_ERB` → ERB (dimensionless)

Requires the 'units' package to be installed.

## Examples

``` r
if (FALSE) { # \dontrun{
df <- data.frame(frame_time = 1:10, fo_Hz = rnorm(10, 120, 10))
df <- .assign_track_units(df)
class(df$fo_Hz)  # "units"
} # }
```

# Clean track names for R data frames

Converts SSFF-style track names (with brackets and hyphens) to
R-friendly column names (with underscores).

## Usage

``` r
.clean_track_names(names)
```

## Arguments

- names:

  Character vector. Track names to clean.

## Value

Character vector of cleaned names.

## Details

Applies the following transformations:

1.  `\[Hz\]` → `_Hz` (brackets to underscores)

2.  `H1-H2` → `H1_H2` (hyphens to underscores)

3.  `(local)` → `_local` (parentheses to underscores)

4.  Multiple spaces → single underscore

This follows the hybrid naming strategy where SSFF files use scientific
notation with brackets, but R data frames use clean underscore notation.

## Examples

``` r
if (FALSE) { # \dontrun{
.clean_track_names(c("fo[Hz]", "F1[Hz]", "H1-H2c[dB]"))
# [1] "fo_Hz" "F1_Hz" "H1_H2c_dB"

.clean_track_names(c("Jitter (local)[%]", "Shimmer (dda)[%]"))
# [1] "Jitter_local_pct" "Shimmer_dda_pct"
} # }
```

# Expand track template to column names

Expands a track name template containing placeholder 'i' to a vector of
numbered column names.

## Usage

``` r
.expand_track_template(template, n_cols)
```

## Arguments

- template:

  Character. Track name template (e.g., "Fi\[Hz\]", "LPCi").

- n_cols:

  Integer. Number of columns to generate.

## Value

Character vector of expanded column names.

## Details

Uses the uniform placeholder pattern where the last 'i' after an
uppercase letter and before '\[' or end-of-string is replaced with
sequential numbers.

For templates without placeholders, falls back to appending "\_1",
"\_2", etc.

## Examples

``` r
if (FALSE) { # \dontrun{
.expand_track_template("Fi[Hz]", 4)
# [1] "F1[Hz]" "F2[Hz]" "F3[Hz]" "F4[Hz]"

.expand_track_template("LPCi", 12)
# [1] "LPC1" "LPC2" "LPC3" ... "LPC12"

.expand_track_template("Hi[dB]", 3)
# [1] "H1[dB]" "H2[dB]" "H3[dB]"

# Non-template (fallback)
.expand_track_template("unknown", 3)
# [1] "unknown_1" "unknown_2" "unknown_3"
} # }
```

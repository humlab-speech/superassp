# Convert Frequency to ERB-rate Scale

Converts frequency values from Hz to the ERB-rate scale (also called
ERBS or Cams). The ERB-rate scale is based on the Equivalent Rectangular
Bandwidth of auditory filters and provides a perceptually uniform
frequency scale.

## Usage

``` r
ucnv_hz_to_erb(freq, method = c("glasberg1990", "moore1983"), as_units = NULL)
```

## Arguments

- freq:

  Numeric vector or units object with frequency values in Hz. If a units
  object is provided, it will be converted to Hz.

- method:

  Character string specifying the conversion formula to use. Options
  are:

  - `"glasberg1990"` (default): Glasberg & Moore (1990) linear
    approximation

  - `"moore1983"`: Moore & Glasberg (1983) polynomial approximation

- as_units:

  Logical. If TRUE, returns a units object with "ERB" units. If FALSE,
  returns a plain numeric vector. Default is TRUE if the units package
  is available and input has units.

## Value

If `as_units = TRUE` and the units package is available: a units object
with "ERB" units. Otherwise, a numeric vector of ERB-rate values.

## Details

Two different approximation formulas are available:

**Glasberg & Moore (1990)** (default): \$\$ERBS(f) = 21.4 \cdot
\log\_{10}(1 + 0.00437 \cdot f)\$\$

Valid range: 100-10,000 Hz. This is the most commonly used formula.

**Moore & Glasberg (1983)**: \$\$ERBS(f) = 11.17 \cdot
\ln\left(\frac{f + 312}{f + 14675}\right) + 43.0\$\$

where f is in Hz. Valid range: 100-6500 Hz.

## References

(Glasberg and Moore 1990)

(Moore and Glasberg 1983)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
ucnv_hz_to_erb(1000)
ucnv_hz_to_erb(c(100, 500, 1000, 2000, 4000))

# With units package
library(units)
freq <- set_units(1000, Hz)
ucnv_hz_to_erb(freq)

# Different methods
ucnv_hz_to_erb(1000, method = "moore1983")
} # }
```

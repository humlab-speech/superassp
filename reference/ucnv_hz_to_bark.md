# Convert Frequency to Bark Scale

Converts frequency values from Hz to the Bark scale, a psychoacoustic
scale that approximates the human ear's critical bands. The Bark scale
is linear up to 1 Bark per critical band (24 critical bands total).

## Usage

``` r
ucnv_hz_to_bark(
  freq,
  method = c("traunmuller", "zwicker", "wang"),
  as_units = NULL
)
```

## Arguments

- freq:

  Numeric vector or units object with frequency values in Hz. If a units
  object is provided, it will be converted to Hz.

- method:

  Character string specifying the conversion formula to use. Options
  are:

  - `"traunmuller"` (default): Traunmüller (1990) approximation

  - `"zwicker"`: Zwicker's formula

  - `"wang"`: Wang, Sekey & Gersho (1992) formula

- as_units:

  Logical. If TRUE, returns a units object with "Bark" units. If FALSE,
  returns a plain numeric vector. Default is TRUE if the units package
  is available and input has units.

## Value

If `as_units = TRUE` and the units package is available: a units object
with "Bark" units. Otherwise, a numeric vector of Bark values.

## Details

Three different approximation formulas are available:

**Traunmüller (1990)** (default): \$\$Bark = \frac{26.81 \cdot f}{1960 +
f} - 0.53\$\$

Valid range: approximately 20-15500 Hz. Values below 2 Bark are adjusted
by \\Bark\_{adjusted} = Bark + 0.15 \cdot (2 - Bark)\\.

**Zwicker**: \$\$Bark = 13 \cdot \arctan(0.00076 \cdot f) + 3.5 \cdot
\arctan\left(\frac{f}{7500}\right)^2\$\$

**Wang, Sekey & Gersho (1992)**: \$\$Bark = 6 \cdot
\sinh^{-1}\left(\frac{f}{600}\right)\$\$

## References

(Traunmüller 1990)

(Zwicker 1961)

(Wang et al. 1992)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage with numeric input
ucnv_hz_to_bark(1000)
ucnv_hz_to_bark(c(100, 500, 1000, 2000, 4000))

# With units package
library(units)
freq <- set_units(c(100, 500, 1000, 2000), Hz)
ucnv_hz_to_bark(freq)

# Different methods
ucnv_hz_to_bark(1000, method = "zwicker")
ucnv_hz_to_bark(1000, method = "wang")

# Without units in output
ucnv_hz_to_bark(1000, as_units = FALSE)
} # }
```

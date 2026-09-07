# Convert Bark Scale to Frequency

Converts values from the Bark scale to frequency in Hz. This is the
inverse of
[`ucnv_hz_to_bark`](https://humlab-speech.github.io/superassp/reference/ucnv_hz_to_bark.md).

## Usage

``` r
ucnv_bark_to_hz(
  bark,
  method = c("traunmuller", "zwicker", "wang"),
  as_units = NULL
)
```

## Arguments

- bark:

  Numeric vector or units object with Bark values. If a units object is
  provided with "Bark" units, it will be converted to numeric.

- method:

  Character string specifying which inverse formula to use. Must match
  the forward conversion method. Options are:

  - `"traunmuller"` (default): Inverse of Traunmüller (1990)

  - `"zwicker"`: Inverse of Zwicker's formula (numerical)

  - `"wang"`: Inverse of Wang, Sekey & Gersho (1992)

- as_units:

  Logical. If TRUE, returns a units object with "Hz" units. If FALSE,
  returns a plain numeric vector. Default is TRUE if the units package
  is available and input has units.

## Value

If `as_units = TRUE` and the units package is available: a units object
with "Hz" units. Otherwise, a numeric vector of frequencies in Hz.

## Details

The inverse formulas are derived analytically where possible:

**Traunmüller (1990)** (default): \$\$f = \frac{1960 \cdot (Bark +
0.53)}{26.81 - (Bark + 0.53)}\$\$

Note: For Bark \< 2, the input is adjusted to account for the
low-frequency correction in the forward transform.

**Wang, Sekey & Gersho (1992)**: \$\$f = 600 \cdot
\sinh\left(\frac{Bark}{6}\right)\$\$

**Zwicker**: Uses numerical root-finding (uniroot) since the inverse is
not analytically solvable.

## See also

[`ucnv_hz_to_bark`](https://humlab-speech.github.io/superassp/reference/ucnv_hz_to_bark.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Basic usage
ucnv_bark_to_hz(10)
ucnv_bark_to_hz(c(5, 10, 15, 20))

# Round trip conversion
freq <- 1000
bark <- ucnv_hz_to_bark(freq)
ucnv_bark_to_hz(bark)  # Should return ~1000

# With units package
library(units)
b <- set_units(10, Bark)
ucnv_bark_to_hz(b)

# Different methods
ucnv_bark_to_hz(10, method = "wang")
} # }
```

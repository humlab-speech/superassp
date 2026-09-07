# Convert track name to plotmath expression

Converts a cleaned track name to a plotmath expression with proper
subscripts for scientific notation.

## Usage

``` r
.col_to_expression(col, full = FALSE)
```

## Arguments

- col:

  Character. Column name (cleaned format, e.g., "fo_Hz", "F1_Hz").

- full:

  Logical. If TRUE, use full descriptive text. Default: FALSE.

## Value

Expression object suitable for ggplot2 labels.

## Details

Creates plotmath expressions with subscripts following Titze 2015
notation:

- `fo_Hz` → expression(f\[o\]~"\[Hz\]") \# f subscript o

- `F1_Hz` → expression(F\[1\]~"\[Hz\]") \# F subscript 1

- `H1_H2c_dB` → expression(H\[1\]-H\["2c"\]~"\[dB\]") \# H subscript 1
  minus H subscript 2c

- `LPC12` → expression(LPC\[12\]) \# LPC subscript 12

The `~` operator adds a small space between the parameter and unit in
plotmath.

For full descriptive labels, returns a character string instead of
expression, as full text like "Frequency of oscillation \[Hz\]" doesn't
need subscripts.

## Examples

``` r
if (FALSE) { # \dontrun{
# Short labels (expressions with subscripts)
.col_to_expression("fo_Hz")       # f\[o\] \[Hz\]
.col_to_expression("F1_Hz")       # F\[1\] \[Hz\]
.col_to_expression("H1_H2c_dB")   # H\[1\]-H\[2c\] \[dB\]

# Full labels (text strings)
.col_to_expression("fo_Hz", full = TRUE)  # "Frequency of oscillation \[Hz\]"
} # }
```

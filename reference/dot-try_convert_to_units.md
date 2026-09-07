# Helper function to try converting column to units

Attempts to convert a numeric column to units. If successful, returns
the units object. If it fails, returns the original column and issues a
warning.

## Usage

``` r
.try_convert_to_units(col, unit_str, col_name)
```

## Arguments

- col:

  Numeric vector; column data

- unit_str:

  Character; unit string to convert to

- col_name:

  Character; column name for warning messages

## Value

Numeric vector or units object

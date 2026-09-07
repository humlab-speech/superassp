# Generate expression labels for data frame columns

Creates a named list of plotmath expressions for all columns in a data
frame.

## Usage

``` r
.generate_track_expressions(col_names, full = FALSE)
```

## Arguments

- col_names:

  Character vector. Column names.

- full:

  Logical. If TRUE, use full descriptive labels. Default: FALSE.

## Value

Named list mapping column names to expressions.

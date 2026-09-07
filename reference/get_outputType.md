# Get the return format of a wrassp/superassp speech signal processing function

Speech signal processing functions in superassp may return an SSFF track
with one or more columns, or a simple, not nested list. This function
may be used to learn the output type of a specific function.

## Usage

``` r
get_outputType(x, package = "superassp")
```

## Arguments

- x:

  The name of a speech signal processing function that are defined in
  the superassp or wrassp packages.

- package:

  The name of the package where the function is defined.

## Value

Either "SSFF" or "list".

## Examples

``` r
superassp:::get_outputType("trk_formant_forest")
#> [1] "SSFF"
superassp:::get_outputType("lst_avqi")
#> [1] "JSTF"
superassp:::get_outputType("trk_formant_burg")
#> [1] "SSFF"
```

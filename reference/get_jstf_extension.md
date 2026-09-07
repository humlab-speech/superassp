# Get JSTF extension for function

Get JSTF extension for function

## Usage

``` r
get_jstf_extension(function_name)
```

## Arguments

- function_name:

  Name of lst\_\* function

## Value

File extension (without dot)

## Examples

``` r
superassp:::get_jstf_extension("lst_avqi")  # "avq"
#> [1] "avq"
superassp:::get_jstf_extension("lst_dsi")  # "dsi"
#> [1] "dsi"
```

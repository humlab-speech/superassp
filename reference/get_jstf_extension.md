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
get_jstf_extension("lst_vat")  # "vat"
#> Error in get_jstf_extension("lst_vat"): could not find function "get_jstf_extension"
get_jstf_extension("lst_voice_sauce")  # "vsj"
#> Error in get_jstf_extension("lst_voice_sauce"): could not find function "get_jstf_extension"
```

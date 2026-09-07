# Get the (default) extension for an SSFF producing function or a signal file

Get the (default) extension for an SSFF producing function or a signal
file

## Usage

``` r
get_extension(x)
```

## Arguments

- x:

  The name of function defined to output trackdata, or the full path to
  a file that can be read using
  [wrassp::read.AsspDataObj](https://rdrr.io/pkg/wrassp/man/read.AsspDataObj.html).

## Value

A string indicating the default file extension of the SSFF generating
function, or the file extension of the signal file.

## Examples

``` r
superassp:::get_extension("trk_formant_forest")
#> [1] "fms"
superassp:::get_extension("trk_formant_burg")
#> [1] "pfm"
```

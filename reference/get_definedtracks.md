# Get available tracks for a function or an SSFF file

Get available tracks for a function or an SSFF file

## Usage

``` r
get_definedtracks(x)
```

## Arguments

- x:

  The name of function defined to output trackdata, or the full path to
  a file that can be read using
  [wrassp::read.AsspDataObj](https://rdrr.io/pkg/wrassp/man/read.AsspDataObj.html).

## Value

A vector of tracks that the function is defined to return, or are
contained within the file.

## Examples

``` r
superassp:::get_definedtracks("trk_formant_forest")
#> [1] "F[Hz]" "B[Hz]"
superassp:::get_definedtracks("trk_formant_burg")
#>  [1] "fm1" "fm2" "fm3" "fm4" "fm5" "bw1" "bw2" "bw3" "bw4" "bw5" "L1"  "L2" 
#> [13] "L3"  "L4"  "L5" 
```

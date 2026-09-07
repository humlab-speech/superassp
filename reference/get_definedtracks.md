# Get available tracks for a function or an SSFF file

Get available tracks for a function or an SSFF file

## Usage

``` r
get_definedtracks(x)
```

## Arguments

- x:

  The name of function defined to output trackdata, or the full path to
  a file that can be read using wrassp::read.AsspDataObj.

## Value

A vector of tracks that the function is defined to return, or are
contained within the file.

## Examples

``` r
get_definedtracks("trk_formant_forest")
#> Error in get_definedtracks("trk_formant_forest"): could not find function "get_definedtracks"
get_definedtracks("trk_formant_burg")
#> Error in get_definedtracks("trk_formant_burg"): could not find function "get_definedtracks"
```

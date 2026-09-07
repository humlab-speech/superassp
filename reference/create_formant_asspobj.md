# Convert Snack formant result to AsspDataObj

Convert Snack formant result to AsspDataObj

## Usage

``` r
create_formant_asspobj(res, windowShift, numFormants)
```

## Arguments

- res:

  List from snackf_cpp (fm, bw, sample_rate, n_frames)

- windowShift:

  Frame shift in milliseconds

- numFormants:

  Number of formants

## Value

AsspDataObj with fm and bw tracks

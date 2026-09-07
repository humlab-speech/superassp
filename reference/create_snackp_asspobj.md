# Convert Snack pitch result (4 tracks) to AsspDataObj

Convert Snack pitch result (4 tracks) to AsspDataObj

## Usage

``` r
create_snackp_asspobj(res, windowShift)
```

## Arguments

- res:

  List from snackp_cpp (f0, voicing, rms, acpeak, sample_rate, n_frames)

- windowShift:

  Frame shift in milliseconds

## Value

AsspDataObj with f0, voicing, rms, acpeak tracks

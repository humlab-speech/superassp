# Convert CheapTrick C++ result to AsspDataObj

Convert CheapTrick C++ result to AsspDataObj

## Usage

``` r
create_spectrogram_asspobj(ct_result, windowShift)
```

## Arguments

- ct_result:

  List returned from cheap_trick_cpp

- windowShift:

  Frame shift in milliseconds

## Value

AsspDataObj with 'sp' track (n_frames x fft_size/2+1 matrix)

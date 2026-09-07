# Convert SPTK REAPER epoch times to pitch mark AsspDataObj

Convert SPTK REAPER epoch times to pitch mark AsspDataObj

## Usage

``` r
create_pitchmark_asspobj(epoch_times, sample_rate, windowShift)
```

## Arguments

- epoch_times:

  Vector of epoch times in seconds (irregular intervals)

- sample_rate:

  Original audio sample rate

- windowShift:

  Frame shift in milliseconds

## Value

AsspDataObj with pm (pitch mark) track as binary indicator

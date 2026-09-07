# Internal: D4C aperiodicity on raw PCM for voxit pipeline

Internal: D4C aperiodicity on raw PCM for voxit pipeline

## Usage

``` r
d4c_r(wave, fs, temporal_positions, f0, threshold = 0.85)
```

## Arguments

- wave:

  Numeric vector; PCM samples

- fs:

  Integer; sample rate

- temporal_positions:

  Numeric vector of frame times

- f0:

  Numeric vector of F0 values

- threshold:

  D4C VUV threshold

## Value

list(aperiodicity, temporal_positions, fs)

# Internal: CheapTrick spectral envelope on raw PCM for voxit pipeline

Internal: CheapTrick spectral envelope on raw PCM for voxit pipeline

## Usage

``` r
cheap_trick_r(wave, fs, temporal_positions, f0, q1 = -0.15, f0_floor = 71)
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

- q1:

  Regularization parameter

- f0_floor:

  Lower F0 bound for FFT size

## Value

list(spectrogram, temporal_positions, fs, fft_size)

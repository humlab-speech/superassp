# Internal: Harvest F0 on raw PCM for voxit pipeline

Internal: Harvest F0 on raw PCM for voxit pipeline

## Usage

``` r
harvest_r(wave, fs, f0_floor = 71, f0_ceil = 800, frame_period = 5)
```

## Arguments

- wave:

  Numeric vector; PCM samples

- fs:

  Integer; sample rate

- f0_floor:

  Numeric; lower F0 bound (Hz)

- f0_ceil:

  Numeric; upper F0 bound (Hz)

- frame_period:

  Numeric; frame shift in ms

## Value

list(f0, temporal_positions, fs)

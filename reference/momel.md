# Run MOMEL algorithm on F0 values

Run MOMEL algorithm on F0 values

## Usage

``` r
momel(
  f0_values,
  window_length = 30L,
  min_f0 = 60,
  max_f0 = 750,
  max_error = 1.04,
  reduced_window_length = 20L,
  minimal_distance = 5,
  minimal_frequency_ratio = 0.05
)
```

## Arguments

- f0_values:

  numeric vector of F0 in Hz (10ms frames), 0 = unvoiced

- window_length:

  window in samples for cible (default 30ms / 10ms = 3)

- min_f0:

  minimum F0 in Hz

- max_f0:

  maximum F0 in Hz

- max_error:

  maximum error ratio

- reduced_window_length:

  window for reduction

- minimal_distance:

  minimal distance in frames

- minimal_frequency_ratio:

  minimal frequency ratio for merging

## Value

data.frame with columns time (frames, 0-based) and frequency (Hz)

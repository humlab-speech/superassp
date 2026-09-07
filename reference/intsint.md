# Run INTSINT algorithm on MOMEL targets

Run INTSINT algorithm on MOMEL targets

## Usage

``` r
intsint(targets)
```

## Arguments

- targets:

  data.frame with time (frames) and frequency (Hz) from momel()

## Value

list(targets_df, range, key) where targets_df has columns time_sec,
tone, target_hz, estimate_hz

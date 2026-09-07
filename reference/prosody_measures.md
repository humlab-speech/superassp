# Compute prosodic measures from audio file or Sound object

Compute prosodic measures from audio file or Sound object

## Usage

``` r
prosody_measures(
  soundPath = NULL,
  sound = NULL,
  minF = 60,
  maxF = 750,
  windowShift = 10
)
```

## Arguments

- soundPath:

  path to WAV file (optional if sound provided)

- sound:

  pladdrr Sound object (optional if soundPath provided)

- minF:

  minimum F0 Hz (default 60)

- maxF:

  maximum F0 Hz (default 750)

- windowShift:

  window shift in ms (default 10)

## Value

named list of ~180 features, or NULL if file \< 1s
